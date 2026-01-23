#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Per-locus MSA pipeline (genome-major inputs) — CLEAN DEFAULTS.

ADDED (as requested):
- After each MAFFT or FAMSA alignment, we apply the "bad alignment" criterion
  equivalent to your R code:
    * if alignment length < 1000 => do NOT flag as bad
    * compute per-column variation (1 if NOT exactly one of A/C/G/T is present)
    * for each sequence, find non-gap blocks (continuous runs where char != '-')
    * for each block compute pi = (#variation cols in block) / block_len
    * if ANY block in ANY sequence has pi > 0.2 => alignment is BAD
  If BAD => output timeout_mark ('*' by default) for ALL genomes at this locus.

Notes:
- Gap char is fixed to '-'
- Variation ignores gaps and ignores non-ACGT letters (e.g. N) the same way
  mx2profile typically profiles A/C/G/T; columns with 0 or >=2 A/C/G/T are "variation".
"""

import argparse
import os
import sys
import subprocess
import shlex
from concurrent.futures import ProcessPoolExecutor
from typing import Dict, List, Tuple, Optional


GAP = "-"  # fixed, no CLI

# --- Bad-alignment criterion params (from your R) ---
SIM_CUTOFF = 0.2
MIN_LEN_CHECK = 1000  # if first aligned seq length < 1000 => do not mark bad


# -------------------------
# CLI (minimal)
# -------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Per-locus MSA (clean defaults)")

    # Required
    p.add_argument("--inputs-list", required=True,
                   help="File with input genome paths, one per line (order preserved).")
    p.add_argument("--outdir", required=True,
                   help="Output directory.")

    # Outputs
    p.add_argument("--out-suffix", default=".aln.txt",
                   help="Suffix for output files (default: .aln.txt).")

    # Parallelism / streaming
    p.add_argument("--block-size", type=int, default=1000,
                   help="Number of loci processed per block.")
    p.add_argument("--threads", type=int, default=8,
                   help="Number of worker processes (parallel loci).")
    p.add_argument("--chunksize", type=int, default=64,
                   help="ProcessPool map chunksize.")
    p.add_argument("--out-stripe", type=int, default=16,
                   help="How many output files to keep open at once.")
    p.add_argument("--progress-every", type=int, default=2000,
                   help="Print progress every N loci (stderr).")
    p.add_argument("--expected-lines", type=int, default=0,
                   help="If >0, enforce exact loci count; else read to EOF.")

    # Normalization
    p.add_argument("--uppercase", action="store_true",
                   help="Uppercase sequences before alignment.")
    p.add_argument("--strip-spaces", action="store_true",
                   help="Strip leading/trailing whitespace per line.")

    # Aligner selection
    p.add_argument("--aligner", choices=["abpoa", "mafft", "famsa", "none"], default="abpoa",
                   help="Aligner for non-trivial loci (default: abpoa).")

    # Timeout / failure marker
    p.add_argument("--timeout-sec", type=int, default=600,
                   help="Per-locus timeout for aligner call in seconds (default: 600).")
    p.add_argument("--timeout-mark", default="*",
                   help="String to output for a timed-out/failed/bad locus (default: '*').")

    # Binaries + optional extras
    p.add_argument("--abpoa-bin", default="abpoa", help="abPOA executable name/path.")
    p.add_argument("--mafft-bin", default="mafft", help="MAFFT executable name/path.")
    p.add_argument("--famsa-bin", default="famsa", help="FAMSA executable name/path.")
    p.add_argument("--famsa-extra", default="", help="Extra FAMSA args (string).")

    # dump problematic loci
    p.add_argument("--dump-fasta-dir", default="",
                   help="If set (non-empty), write locus_X.fasta for loci that output timeout_mark ('*').")


    return p.parse_args()


def read_inputs_list(path: str) -> List[str]:
    with open(path, "r", encoding="utf-8") as f:
        files = [line.strip() for line in f if line.strip()]
    if not files:
        raise ValueError("inputs-list is empty.")
    return files


def ensure_outdir(outdir: str) -> None:
    os.makedirs(outdir, exist_ok=True)


def out_path(outdir: str, in_path: str, suffix: str) -> str:
    return os.path.join(outdir, os.path.basename(in_path) + suffix)


# -------------------------
# FASTA / MSA parsing helpers
# -------------------------

def parse_fasta(text: str) -> Tuple[List[str], List[str]]:
    """Parse FASTA. Returns (ids, seqs) preserving order."""
    ids: List[str] = []
    seqs: List[str] = []
    cur_id: Optional[str] = None
    cur_chunks: List[str] = []

    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if cur_id is not None:
                ids.append(cur_id)
                seqs.append("".join(cur_chunks))
            cur_id = line[1:].strip()
            cur_chunks = []
        else:
            cur_chunks.append(line)

    if cur_id is not None:
        ids.append(cur_id)
        seqs.append("".join(cur_chunks))

    if not seqs:
        raise RuntimeError("Empty FASTA output.")
    return ids, seqs


def parse_abpoa_msa(text: str) -> Tuple[List[str], List[str]]:
    """
    Parse abPOA MSA output (often PIR-like). Supports PIR-ish and FASTA-ish.
    Returns (ids, aligned_seqs).
    """
    lines = [ln.rstrip("\n") for ln in text.splitlines() if ln.strip() != ""]
    ids: List[str] = []
    seqs: List[str] = []

    i = 0
    while i < len(lines):
        ln = lines[i]
        if not ln.startswith(">"):
            i += 1
            continue

        header = ln[1:]
        name = header.strip()

        # PIR headers: P1;name or F1;name, followed by a description line
        if header.startswith("P1;") or header.startswith("F1;"):
            name = header.split(";", 1)[1].strip()
            i += 1
            if i < len(lines):
                i += 1
        else:
            i += 1

        seq_chunks: List[str] = []
        while i < len(lines) and not lines[i].startswith(">"):
            s = lines[i].strip()
            if "*" in s:
                s = s.replace("*", "")
                if s:
                    seq_chunks.append(s)
                i += 1
                break
            seq_chunks.append(s)
            i += 1

        ids.append(name)
        seqs.append("".join(seq_chunks))

    if not seqs:
        raise RuntimeError("Failed to parse abPOA output (empty MSA).")
    return ids, seqs


def _map_back_by_numeric_ids(ids: List[str], aln: List[str], n: int) -> List[str]:
    """
    We build FASTA with numeric ids 0..k-1.
    Map alignments back into that order; fallback to first n if parsing fails.
    """
    out: List[Optional[str]] = [None] * n
    for name, aseq in zip(ids, aln):
        try:
            j = int(name)
        except ValueError:
            j = None
        if j is not None and 0 <= j < n:
            out[j] = aseq
    if any(x is None for x in out):
        return aln[:n]
    return [x for x in out if x is not None]


# -------------------------
# BAD alignment criterion (R -> Python)
# -------------------------

def _find_ones_blocks(bits: List[int]) -> List[Tuple[int, int]]:
    """Equivalent of findOnes() for a 0/1 list. Returns [(beg,end), ...] inclusive indices."""
    blocks: List[Tuple[int, int]] = []
    i = 0
    n = len(bits)
    while i < n:
        if bits[i] == 1:
            beg = i
            i += 1
            while i < n and bits[i] == 1:
                i += 1
            end = i - 1
            blocks.append((beg, end))
        else:
            i += 1
    return blocks


def is_bad_alignment_like_r(aligned: List[str],
                            sim_cutoff: float = SIM_CUTOFF,
                            min_len_check: int = MIN_LEN_CHECK) -> bool:
    """
    Replicates the R criterion:

    - if nchar(seqs[1]) < 1000: return FALSE (not bad)
    - pos_variation[col] = 1 if NOT exactly one of A/C/G/T is present in that column
      (len(present_nucs) != 1)
    - for each sequence: blocks = contiguous positions where s != '-'
      pi(block) = sum(pos_variation[block]) / block_len
    - if any block has pi > sim_cutoff => BAD
    """
    if not aligned:
        return False
    L = len(aligned[0])
    if L < min_len_check:
        return False

    # safety: ensure all same length
    for s in aligned:
        if len(s) != L:
            # If aligner produced inconsistent lengths, treat as bad (safer)
            return True

    # uppercase (R code: toupper(mx))
    aligned_u = [s.upper() for s in aligned]

    # Compute per-column variation
    # R logic corresponds to: variation = 0 only when exactly 1 of A/C/G/T is present
    pos_variation = [0] * L
    acgt = {"A", "C", "G", "T"}
    for j in range(L):
        present = set()
        for s in aligned_u:
            ch = s[j]
            if ch in acgt:
                present.add(ch)
                if len(present) > 1:
                    break
        pos_variation[j] = 0 if len(present) == 1 else 1

    # For each sequence, find non-gap blocks and compute pi
    for s in aligned_u:
        bits = [1 if ch != GAP else 0 for ch in s]
        blocks = _find_ones_blocks(bits)
        if not blocks:
            continue
        for beg, end in blocks:
            blen = end - beg + 1
            # pi = fraction of variable columns in that block
            v = 0
            # tight loop
            for k in range(beg, end + 1):
                v += pos_variation[k]
            pi = v / blen
            if pi > sim_cutoff:
                return True

    return False


# -------------------------
# Aligner runners (stdin + timeout)
# -------------------------

def _build_fasta(unique_seqs: List[str]) -> str:
    lines: List[str] = []
    for idx, s in enumerate(unique_seqs):
        lines.append(f">{idx}")
        lines.append(s)
    return "\n".join(lines) + "\n"


def run_abpoa_msa_stdin(unique_seqs: List[str], abpoa_bin: str, timeout_sec: int) -> List[str]:
    """
    Fixed abPOA defaults (as requested):
      - global mode: -m 0
      - match/mismatch: -M 2 -X 4
      - gap penalties: -O 20,60  -E 1,0
    """
    fasta_text = _build_fasta(unique_seqs)
    cmd = [
        abpoa_bin,
        "-r", "1",
        "-m", "0",
        "-M", "2",
        "-X", "4",
        "-O", "20,60",
        "-E", "1,0",
        "-"
    ]
    proc = subprocess.run(cmd, input=fasta_text, capture_output=True, text=True, timeout=timeout_sec)
    if proc.returncode != 0:
        raise RuntimeError(f"abPOA failed (code {proc.returncode}). stderr:\n{proc.stderr[:2000]}")
    ids, aln = parse_abpoa_msa(proc.stdout)
    return _map_back_by_numeric_ids(ids, aln, len(unique_seqs))


def run_mafft_msa_stdin(unique_seqs: List[str], mafft_bin: str, timeout_sec: int) -> List[str]:
    """
    Fixed MAFFT defaults (as requested):
      --quiet --op 3 --ep 0.1 --thread 1
    """
    fasta_text = _build_fasta(unique_seqs)
    cmd = [
        mafft_bin,
        "--quiet",
        "--thread", "1",
        "--op", "3",
        "--ep", "0.1",
        "-"
    ]
    proc = subprocess.run(cmd, input=fasta_text, capture_output=True, text=True, timeout=timeout_sec)
    if proc.returncode != 0:
        raise RuntimeError(f"MAFFT failed (code {proc.returncode}). stderr:\n{proc.stderr[:2000]}")
    ids, aln = parse_fasta(proc.stdout)
    return _map_back_by_numeric_ids(ids, aln, len(unique_seqs))


def run_famsa_msa_stdin(unique_seqs: List[str], famsa_bin: str, extra: str, timeout_sec: int) -> List[str]:
    """
    FAMSA reads FASTA from stdin and writes aligned FASTA to stdout.
    Threading: keep single-thread inside each locus; parallelize loci via --threads.
    """
    fasta_text = _build_fasta(unique_seqs)
    cmd = [famsa_bin]
    if extra.strip():
        cmd += shlex.split(extra)
    proc = subprocess.run(cmd, input=fasta_text, capture_output=True, text=True, timeout=timeout_sec)
    if proc.returncode != 0:
        raise RuntimeError(f"FAMSA failed (code {proc.returncode}). stderr:\n{proc.stderr[:2000]}")
    ids, aln = parse_fasta(proc.stdout)
    return _map_back_by_numeric_ids(ids, aln, len(unique_seqs))


# -------------------------
# Worker config + per-locus logic
# -------------------------

class _Cfg:
    __slots__ = (
        "uppercase",
        "aligner",
        "timeout_sec",
        "timeout_mark",
        "abpoa_bin",
        "mafft_bin",
        "famsa_bin",
        "famsa_extra",
    )

    def __init__(self, args: argparse.Namespace):
        self.uppercase = bool(args.uppercase)
        self.aligner = args.aligner
        self.timeout_sec = int(args.timeout_sec)
        self.timeout_mark = args.timeout_mark

        self.abpoa_bin = args.abpoa_bin
        self.mafft_bin = args.mafft_bin
        self.famsa_bin = args.famsa_bin
        self.famsa_extra = args.famsa_extra


_CFG: Optional[_Cfg] = None


def _init_worker(cfg: _Cfg) -> None:
    global _CFG
    _CFG = cfg


def pad_to_maxlen(seqs: List[str]) -> List[str]:
    if not seqs:
        return []
    max_len = max(len(s) for s in seqs)
    return [s + (GAP * (max_len - len(s))) for s in seqs]


def align_one_locus(locus_lines: List[str]) -> List[str]:
    """
    locus_lines: list[str] length = n_genomes; empty string = missing.
    Output: list[str] length = n_genomes; missing genomes -> gaps of alignment length (or "" if length=0).

    Trivial loci (skip aligner):
      - no present sequences, OR
      - all present sequences are unique
    Non-trivial -> chosen aligner, with timeout. On timeout/failure -> timeout_mark for all genomes.

      - For MAFFT and FAMSA results, apply bad-alignment criterion; if bad => timeout_mark for all genomes.
    """
    global _CFG
    if _CFG is None:
        raise RuntimeError("Worker config not initialized.")

    n = len(locus_lines)

    seq_to_idxs: Dict[str, List[int]] = {}
    present_count = 0

    for i, s in enumerate(locus_lines):
        if s == "":
            continue
        present_count += 1
        if _CFG.uppercase:
            s = s.upper()
        seq_to_idxs.setdefault(s, []).append(i)

    # Trivial A: no present sequences => empty line
    if present_count == 0:
        return [""] * n

    unique_seqs = list(seq_to_idxs.keys())

    # Trivial B: all present are unique (no duplicates among present)
    all_unique_present = (present_count == len(unique_seqs))
    if all_unique_present:
        padded_unique = pad_to_maxlen(unique_seqs)
        aln_len = len(padded_unique[0]) if padded_unique else 0
        gaps = GAP * aln_len
        out = [gaps] * n
        for seq, aseq in zip(unique_seqs, padded_unique):
            for gi in seq_to_idxs[seq]:
                out[gi] = aseq
        return out

    # Non-trivial: run aligner with timeout
    try:
        used_aligner = _CFG.aligner
        if used_aligner == "abpoa":
            aligned_unique = run_abpoa_msa_stdin(
                unique_seqs,
                abpoa_bin=_CFG.abpoa_bin,
                timeout_sec=_CFG.timeout_sec,
            )
        elif used_aligner == "mafft":
            aligned_unique = run_mafft_msa_stdin(
                unique_seqs,
                mafft_bin=_CFG.mafft_bin,
                timeout_sec=_CFG.timeout_sec,
            )
        elif used_aligner == "famsa":
            aligned_unique = run_famsa_msa_stdin(
                unique_seqs,
                famsa_bin=_CFG.famsa_bin,
                extra=_CFG.famsa_extra,
                timeout_sec=_CFG.timeout_sec,
            )
        else:  # none
            aligned_unique = pad_to_maxlen(unique_seqs)

    except subprocess.TimeoutExpired:
        return [_CFG.timeout_mark] * n
    except Exception as e:
        print(f"[align-error] {type(e).__name__}: {e}", file=sys.stderr)
        raise

    # apply bad-alignment criterion only after MAFFT or FAMSA
    if _CFG.aligner in ("mafft", "famsa"):
        if is_bad_alignment_like_r(aligned_unique, sim_cutoff=SIM_CUTOFF, min_len_check=MIN_LEN_CHECK):
            return [_CFG.timeout_mark] * n

    aln_len = len(aligned_unique[0]) if aligned_unique else 0
    gaps = GAP * aln_len

    out = [gaps] * n
    for seq, aseq in zip(unique_seqs, aligned_unique):
        for gi in seq_to_idxs[seq]:
            out[gi] = aseq
    return out


# -------------------------
# IO: striped write
# -------------------------

def write_block_striped(out_paths: List[str], per_genome_lines: List[List[str]], stripe: int) -> None:
    n = len(out_paths)
    for start in range(0, n, stripe):
        end = min(n, start + stripe)
        fhs = []
        try:
            for i in range(start, end):
                fhs.append(open(out_paths[i], "a", encoding="utf-8", newline=""))
            for j, fh in enumerate(fhs):
                gi = start + j
                fh.writelines(per_genome_lines[gi])
        finally:
            for fh in fhs:
                try:
                    fh.close()
                except Exception:
                    pass


# -------------------------
# dump problematic loci as FASTA
# -------------------------

def dump_locus_fasta(dump_dir: str,
                     locus_number_1based: int,
                     locus_lines: List[str],
                     genome_names: List[str],
                     uppercase: bool) -> None:
    """
    Write locus_X.fasta with original sequences.
    - headers are genome_names (basenames)
    - empty strings are skipped
    """
    if not dump_dir:
        return

    out_path = os.path.join(dump_dir, f"locus_{locus_number_1based}.fasta")
    tmp_path = out_path + ".tmp"

    with open(tmp_path, "w", encoding="utf-8", newline="\n") as f:
        wrote = False
        for name, seq in zip(genome_names, locus_lines):
            if seq == "":
                continue
            s = seq.upper() if uppercase else seq
            f.write(f">{name}\n{s}\n")
            wrote = True

    # If everything was empty, remove file (per your requirement: empty lines not represented)
    if not os.path.getsize(tmp_path):
        try:
            os.remove(tmp_path)
        except Exception:
            pass
        return

    os.replace(tmp_path, out_path)


# -------------------------
# Main
# -------------------------

def main() -> int:
    args = parse_args()
    in_files = read_inputs_list(args.inputs_list)
    ensure_outdir(args.outdir)

    if args.dump_fasta_dir:
        ensure_outdir(args.dump_fasta_dir)

    if args.threads < 1:
        raise ValueError("--threads must be >= 1")
    if args.block_size < 1:
        raise ValueError("--block-size must be >= 1")
    if args.out_stripe < 1:
        raise ValueError("--out-stripe must be >= 1")
    if args.timeout_sec < 1:
        raise ValueError("--timeout-sec must be >= 1")
    if args.timeout_mark == "":
        raise ValueError("--timeout-mark must be non-empty (e.g. '*').")

    n_genomes = len(in_files)
    genome_names = [os.path.basename(p) for p in in_files]

    # Open inputs
    infhs = [open(p, "r", encoding="utf-8", newline="") for p in in_files]
    out_paths = [out_path(args.outdir, p, args.out_suffix) for p in in_files]

    # Truncate/create outputs
    for op in out_paths:
        with open(op, "w", encoding="utf-8", newline=""):
            pass

    cfg = _Cfg(args)
    locus_idx = 0  # number of loci read so far (1-based will be locus_idx)
    block: List[List[str]] = []
    block_nums: List[int] = []  # 1-based locus numbers aligned with block entries

    pool = ProcessPoolExecutor(
        max_workers=args.threads,
        initializer=_init_worker,
        initargs=(cfg,),
    )

    try:
        while True:
            locus_lines: List[Optional[str]] = []
            any_eof = False

            for fh in infhs:
                s = fh.readline()
                if s == "":
                    any_eof = True
                    locus_lines.append(None)
                else:
                    s = s.rstrip("\n")
                    if args.strip_spaces:
                        s = s.strip()
                    locus_lines.append(s)

            if any_eof:
                if not all(x is None for x in locus_lines):
                    raise RuntimeError("Input files have different number of lines (EOF mismatch).")
                break

            # Now it's a real locus
            locus_idx += 1
            locus_line_strs = [x for x in locus_lines if x is not None]  # type: ignore

            if args.expected_lines > 0 and locus_idx > args.expected_lines:
                raise RuntimeError(f"Read more than expected-lines={args.expected_lines}.")

            block.append(locus_line_strs)
            block_nums.append(locus_idx)

            if len(block) >= args.block_size:
                aligned_block = list(pool.map(align_one_locus, block, chunksize=args.chunksize))

                # dump problematic loci (all outputs are timeout_mark)
                if args.dump_fasta_dir:
                    for loc_no, orig_lines, locus_out in zip(block_nums, block, aligned_block):
                        if locus_out and all(x == args.timeout_mark for x in locus_out):
                            dump_locus_fasta(
                                dump_dir=args.dump_fasta_dir,
                                locus_number_1based=loc_no,
                                locus_lines=orig_lines,
                                genome_names=genome_names,
                                uppercase=args.uppercase,
                            )

                per_genome_buf: List[List[str]] = [[] for _ in range(n_genomes)]
                for locus_out in aligned_block:
                    for gi, s2 in enumerate(locus_out):
                        per_genome_buf[gi].append(s2 + "\n")

                write_block_striped(out_paths, per_genome_buf, args.out_stripe)
                block.clear()
                block_nums.clear()

                if args.progress_every and locus_idx % args.progress_every == 0:
                    print(f"[progress] loci processed: {locus_idx}", file=sys.stderr)

        # flush tail
        if block:
            aligned_block = list(pool.map(align_one_locus, block, chunksize=args.chunksize))

            if args.dump_fasta_dir:
                for loc_no, orig_lines, locus_out in zip(block_nums, block, aligned_block):
                    if locus_out and all(x == args.timeout_mark for x in locus_out):
                        dump_locus_fasta(
                            dump_dir=args.dump_fasta_dir,
                            locus_number_1based=loc_no,
                            locus_lines=orig_lines,
                            genome_names=genome_names,
                            uppercase=args.uppercase,
                        )

            per_genome_buf: List[List[str]] = [[] for _ in range(n_genomes)]
            for locus_out in aligned_block:
                for gi, s2 in enumerate(locus_out):
                    per_genome_buf[gi].append(s2 + "\n")
            write_block_striped(out_paths, per_genome_buf, args.out_stripe)

        if args.expected_lines > 0 and locus_idx != args.expected_lines:
            raise RuntimeError(f"Expected exactly {args.expected_lines} loci, but read {locus_idx}.")

        print(f"[done] genomes={n_genomes}, loci={locus_idx}, outdir={args.outdir}", file=sys.stderr)
        return 0

    finally:
        try:
            pool.shutdown(wait=True, cancel_futures=False)
        except Exception:
            pass
        for fh in infhs:
            try:
                fh.close()
            except Exception:
                pass


if __name__ == "__main__":
    raise SystemExit(main())
