#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import re
import sys
from pathlib import Path

# Uniform logging shared with the other pangen steps (see inst/utils/comb_logging.py)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils"))
from comb_logging import Logger, add_log_args

LOCUS_RE = re.compile(
    r"^locus_(\d+)(?:_aligned)?\.fasta$",
    re.IGNORECASE
)

def parse_input_paths(input_txt: Path) -> list[Path]:
    """
    input.txt contains paths, one per line. Genome name = basename of the path.
    Order is preserved. Duplicates (by basename) are removed, first one kept.
    """
    paths: list[Path] = []
    seen = set()
    with input_txt.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            base = os.path.basename(line.rstrip("/\\"))
            if base and base not in seen:
                paths.append(Path(line))
                seen.add(base)
    return paths

def count_lines(path: Path) -> int:
    """
    Number of lines in a file, counted over raw bytes (these files hold one
    long sequence per line, so decoding them just to count would be wasteful).
    A missing trailing newline still counts as a line.
    """
    n = 0
    last = b""
    with path.open("rb") as f:
        while True:
            chunk = f.read(1 << 20)
            if not chunk:
                break
            n += chunk.count(b"\n")
            last = chunk[-1:]
    if last and last != b"\n":
        n += 1
    return n

def expected_locus_count(paths: list[Path]) -> int:
    """
    The per-genome files fed to comb_07 have exactly one line per locus, so
    their line count is the authoritative number of loci for this chromosome.
    Returns 0 if none of them can be read.
    """
    readable = [p for p in paths if p.is_file()]
    if not readable:
        return 0
    # All of them have the same length (comb_07 refuses to run otherwise),
    # so the smallest one on disk is the cheapest to scan.
    return count_lines(min(readable, key=lambda p: p.stat().st_size))

def index_locus_files(folder: Path) -> dict[int, Path]:
    """
    Returns mapping: locus_number -> filepath.

    Both locus_<n>.fasta (comb_07) and locus_<n>_aligned.fasta (comb_08) map to
    the same locus. In a clean run only one of the two exists for a given <n>;
    if both are there, the comb_08 output wins, so that the choice never
    depends on the order the directory happens to be listed in.
    """
    idx: dict[int, Path] = {}
    aligned: set[int] = set()
    for p in sorted(folder.iterdir()):
        if not p.is_file():
            continue
        m = LOCUS_RE.match(p.name)
        if not m:
            continue
        n = int(m.group(1))
        is_aligned = m.group(0).lower().endswith("_aligned.fasta")
        if n in idx and not is_aligned and n in aligned:
            continue
        idx[n] = p
        if is_aligned:
            aligned.add(n)
    return idx

def read_one_line_fasta(path: Path) -> dict[str, str]:
    """
    FASTA where each header >name is followed by exactly ONE sequence line.
    Only the first sequence line after the header is taken
    (extra lines, if any, are ignored).
    """
    data: dict[str, str] = {}
    current = None
    expecting_seq = False

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].strip() or None
                expecting_seq = True
                continue
            if expecting_seq and current is not None:
                data[current] = line
                expecting_seq = False

    return data

def alignment_length(locus_map: dict[str, str]) -> int:
    for seq in locus_map.values():
        return len(seq)
    return 0

def safe_filename(name: str) -> str:
    return re.sub(r'[<>:"/\\|?*\x00-\x1F]', "_", name)

def main():
    ap = argparse.ArgumentParser(
        description=(
            "Transpose locus_*.fasta files into genome-based files; "
            "missing loci -> empty line, "
            "missing genomes in a locus -> gap string of required length. "
            "Loci are processed in batches."
        )
    )
    ap.add_argument("-i", "--input", default="input.txt", help="File containing path list (basename = genome name).")
    ap.add_argument("-d", "--dir", default=".", help="Directory containing locus_*.fasta files.")
    ap.add_argument("-o", "--out", default="out", help="Output directory for genome files.")
    ap.add_argument("-b", "--batch", type=int, default=1000, help="Locus batch size (default 1000).")
    add_log_args(ap)  # logging (uniform with the R steps)
    args = ap.parse_args()

    log = Logger(args.path_log, args.log_level, script=os.path.basename(__file__))

    folder = Path(args.dir).resolve()
    input_txt = Path(args.input).resolve()
    out_dir = Path(args.out).resolve()

    if not input_txt.exists():
        log.log(f"ERROR: input.txt not found: {input_txt}", echo=True)
        print(f"ERROR: input.txt not found: {input_txt}", file=sys.stderr)
        sys.exit(1)

    input_paths = parse_input_paths(input_txt)
    genomes = [p.name for p in input_paths]
    if not genomes:
        log.log("ERROR: input.txt does not contain valid paths/genome names.", echo=True)
        print("ERROR: input.txt does not contain valid paths/genome names.", file=sys.stderr)
        sys.exit(1)

    locus_idx = index_locus_files(folder)
    if not locus_idx:
        log.log(f"ERROR: no files matching locus_*.fasta found in folder {folder}", echo=True)
        print(f"ERROR: no files matching locus_*.fasta found in folder {folder}", file=sys.stderr)
        sys.exit(1)

    # The number of loci comes from the input files, not from whatever happens
    # to lie in the folder. Locus numbering is derived from the breaks and is
    # recomputed on every run, so files left by an earlier run are not
    # overwritten and would otherwise be pulled into the alignment as real loci.
    max_locus = expected_locus_count(input_paths)
    if max_locus <= 0:
        max_locus = max(locus_idx.keys())
        log.log("WARNING: could not read the input files, falling back to the "
                "folder content to get the locus range", echo=True)
    else:
        stale = [n for n in locus_idx if n > max_locus]
        if stale:
            for n in stale:
                del locus_idx[n]
            msg = ("WARNING: %d locus file(s) beyond locus %d ignored in %s "
                   "(leftovers of a previous run; rerun with -cleanup)"
                   % (len(stale), max_locus, folder))
            log.log(msg, echo=True)
            print(msg, file=sys.stderr)

    out_dir.mkdir(parents=True, exist_ok=True)

    log.log("input:   %s" % input_txt)
    log.log("dir:     %s" % folder)
    log.log("out:     %s" % out_dir)
    log.log("genomes: %d  locus files: %d  locus range: 1..%d"
            % (len(genomes), len(locus_idx), max_locus))

    # Prepare output files: create/clear them so we only append later.
    out_paths: dict[str, Path] = {}
    for genome in genomes:
        out_path = out_dir / f"{safe_filename(genome)}.txt"
        out_paths[genome] = out_path
        with out_path.open("w", encoding="utf-8", newline="\n") as out:
            pass

    batch_size = max(1, args.batch)
    missing_files_total = 0

    # Process loci in batches
    for start in range(1, max_locus + 1, batch_size):
        end = min(max_locus, start + batch_size - 1)

        # 1) read locus batch into memory
        locus_maps: dict[int, dict[str, str] | None] = {}
        locus_lens: dict[int, int] = {}

        for i in range(start, end + 1):
            fp = locus_idx.get(i)
            if fp is None:
                locus_maps[i] = None
                locus_lens[i] = 0
                missing_files_total += 1
            else:
                mp = read_one_line_fasta(fp)
                locus_maps[i] = mp
                locus_lens[i] = alignment_length(mp)

        # 2) append lines for each genome for this batch
        for genome in genomes:
            out_path = out_paths[genome]
            with out_path.open("a", encoding="utf-8", newline="\n") as out:
                for i in range(start, end + 1):
                    mp = locus_maps[i]
                    if mp is None:
                        out.write("\n")  # empty line
                    else:
                        if genome in mp:
                            out.write(mp[genome] + "\n")
                        else:
                            L = locus_lens[i]
                            out.write(("-" * L if L > 0 else "") + "\n")

        # optional progress
        print(f"OK: processed locus batch {start}..{end}")
        log.log("progress: locus batch %d..%d" % (start, end), echo=False)

    print(f"OK: files created: {len(genomes)}")
    print(f"Locus range: 1..{max_locus}")
    if missing_files_total:
        print(f"Missing locus_*.fasta files in range: {missing_files_total} (empty lines inserted)")
    print(f"Output directory: {out_dir}")

    log.log("DONE: genome files=%d  locus range=1..%d  missing locus files=%d  elapsed=%.1fs"
            % (len(genomes), max_locus, missing_files_total, log.elapsed()), echo=True)

if __name__ == "__main__":
    main()