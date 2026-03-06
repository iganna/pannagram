#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Tuple, Optional, Dict

GAP = "-"
SIM_CUTOFF = 0.2
MIN_LEN_CHECK = 1000  # if alignment length < 1000 => do NOT mark bad


# -------------------------
# FASTA helpers
# -------------------------

def build_fasta(headers: List[str], seqs: List[str]) -> str:
    out = []
    for h, s in zip(headers, seqs):
        out.append(f">{h}")
        out.append(s)
    return "\n".join(out) + "\n"

def parse_fasta(text: str) -> Tuple[List[str], List[str]]:
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
    return ids, seqs


# -------------------------
# BAD alignment criterion
# -------------------------

def _find_ones_blocks(bits: List[int]) -> List[Tuple[int, int]]:
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
    if not aligned:
        return False
    L = len(aligned[0])
    if L < min_len_check:
        return False

    for s in aligned:
        if len(s) != L:
            return True

    aligned_u = [s.upper() for s in aligned]
    acgt = {"A", "C", "G", "T"}

    pos_variation = [0] * L
    for j in range(L):
        present = set()
        for s in aligned_u:
            ch = s[j]
            if ch in acgt:
                present.add(ch)
                if len(present) > 1:
                    break
        pos_variation[j] = 0 if len(present) == 1 else 1

    for s in aligned_u:
        bits = [1 if ch != GAP else 0 for ch in s]
        blocks = _find_ones_blocks(bits)
        for beg, end in blocks:
            blen = end - beg + 1
            v = sum(pos_variation[beg:end + 1])
            if (v / blen) > sim_cutoff:
                return True
    return False


# -------------------------
# Aligner runners
# -------------------------

def run_mafft(headers: List[str], seqs: List[str],
              mafft_bin: str, timeout_sec: int) -> Tuple[List[str], List[str]]:
    # MAFFT defaults you wanted:
    # --quiet --op 3 --ep 0.1 --thread 1
    fasta_in = build_fasta(headers, seqs)
    cmd = [mafft_bin, "--quiet", "--thread", "1", "--op", "3", "--ep", "0.1", "-"]
    proc = subprocess.run(cmd, input=fasta_in, text=True, capture_output=True, timeout=timeout_sec)
    if proc.returncode != 0:
        raise RuntimeError(f"MAFFT failed (code {proc.returncode}). stderr:\n{proc.stderr[:2000]}")
    return parse_fasta(proc.stdout)

def run_famsa(headers: List[str], seqs: List[str],
              famsa_bin: str, timeout_sec: int) -> Tuple[List[str], List[str]]:
    """
    FAMSA CLI differs between builds/packages.
    We'll try several stdin/stdout conventions until one works.
    """
    fasta_in = build_fasta(headers, seqs)

    # Try multiple invocation styles (most robust in practice)
    cmd_variants = [
        [famsa_bin, "-t", "1"],                            # some builds: stdin->stdout by default
        [famsa_bin, "-t", "1", "-", "-"],                  # dash convention
        [famsa_bin, "-t", "1", "STDIN", "STDOUT"],         # per README (may fail on your build)
        [famsa_bin, "-t", "1", "stdin", "stdout"],         # lowercase variant
        [famsa_bin, "-t", "1", "/dev/stdin", "/dev/stdout"]# explicit devices (Linux)
    ]

    last_err = None
    for cmd in cmd_variants:
        try:
            proc = subprocess.run(
                cmd,
                input=fasta_in,
                text=True,
                capture_output=True,
                timeout=timeout_sec
            )
        except subprocess.TimeoutExpired:
            # propagate timeout to be handled by caller same as MAFFT
            raise

        if proc.returncode == 0 and proc.stdout.strip():
            return parse_fasta(proc.stdout)

        # keep last stderr for debugging if all variants fail
        last_err = (cmd, proc.returncode, proc.stderr)

        # If it printed "Unable to open input file ..." or similar, try next variant
        continue

    cmd, code, err = last_err if last_err is not None else (["<unknown>"], -1, "")
    raise RuntimeError(
        f"FAMSA failed for all stdin/stdout conventions. "
        f"Last cmd={cmd} code={code}. stderr:\n{(err or '')[:2000]}"
    )

# -------------------------
# IO utils
# -------------------------

def read_inputs_list(path: str) -> List[str]:
    with open(path, "r", encoding="utf-8") as f:
        files = [ln.strip() for ln in f if ln.strip()]
    if not files:
        raise ValueError("inputs-list is empty")
    return files

def ensure_dir(p: str) -> None:
    os.makedirs(p, exist_ok=True)

def write_locus_fasta(out_path: str, headers: List[str], seqs: List[str]) -> None:
    tmp = out_path + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="\n") as f:
        f.write(build_fasta(headers, seqs))
    os.replace(tmp, out_path)


# -------------------------
# Worker
# -------------------------

def process_one_locus(locus_no: int,
                      genome_names: List[str],
                      locus_lines: List[str],
                      outdir: str,
                      baddir: str,
                      aligner: str,
                      mafft_bin: str,
                      famsa_bin: str,
                      timeout_sec: int,
                      uppercase: bool,
                      strip_spaces: bool) -> Tuple[int, str]:

    # Collect only non-empty ones.
    headers: List[str] = []
    seqs: List[str] = []
    for name, s in zip(genome_names, locus_lines):
        if strip_spaces:
            s = s.strip()
        if not s:
            continue
        if uppercase:
            s = s.upper()
        headers.append(name)
        seqs.append(s)

    out_path = os.path.join(outdir, f"locus_{locus_no}.fasta")
    bad_path = os.path.join(baddir, f"locus_{locus_no}.fasta")

    # If all are empty, we write an empty locus and OK.
    if not seqs:
        write_locus_fasta(out_path, [], [])
        return locus_no, "ok"

    # SAFE IDS (that aligner don't change headers)
    safe_headers = [f"s{i+1:06d}" for i in range(len(headers))]
    safe2orig = dict(zip(safe_headers, headers))
    orig2safe = dict(zip(headers, safe_headers))

    try:
        if aligner == "mafft":
            ids, aln = run_mafft(safe_headers, seqs, mafft_bin=mafft_bin, timeout_sec=timeout_sec)
        else:
            ids, aln = run_famsa(safe_headers, seqs, famsa_bin=famsa_bin, timeout_sec=timeout_sec)

        m_safe: Dict[str, str] = {i: a for i, a in zip(ids, aln)}

        # We strictly require that all non-empty ones are returned.
        missing_safe = [sh for sh in safe_headers if sh not in m_safe]
        if missing_safe:
            ex = [safe2orig[x] for x in missing_safe[:5]]
            raise RuntimeError(f"Aligner output missing {len(missing_safe)} sequences, e.g. {ex}")

        # restore the order and the original headers
        aligned_in_input_order = [m_safe[orig2safe[h]] for h in headers]

        # write_locus_fasta(out_path, headers, aligned_in_input_order)

        if is_bad_alignment_like_r(aligned_in_input_order):
            # write into bad initial sequences
            write_locus_fasta(bad_path, headers, seqs)
            return locus_no, "bad"

        write_locus_fasta(out_path, headers, aligned_in_input_order)
        return locus_no, "ok"

    except subprocess.TimeoutExpired:
        write_locus_fasta(bad_path, headers, seqs)
        return locus_no, "bad"
    except Exception as e:
        print(f"[locus {locus_no}] error: {type(e).__name__}: {e}", file=sys.stderr)
        write_locus_fasta(bad_path, headers, seqs)
        return locus_no, "bad"


# -------------------------
# Main
# -------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser("Per-locus MSA (mafft|famsa) -> per-locus FASTA outputs")
    p.add_argument("--inputs-list", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--baddir", required=True)

    p.add_argument("--aligner", choices=["mafft", "famsa"], default="mafft",
                   help="Which aligner to run per locus (default: mafft).")

    p.add_argument("--mafft-bin", default="mafft")
    p.add_argument("--famsa-bin", default="famsa")

    p.add_argument("--timeout-sec", type=int, default=180)
    p.add_argument("--threads", type=int, default=8,
                   help="Parallel loci (aligner itself is forced single-thread).")
    p.add_argument("--uppercase", action="store_true")
    p.add_argument("--strip-spaces", action="store_true")
    p.add_argument("--expected-lines", type=int, default=0)
    p.add_argument("--progress-every", type=int, default=2000)
    return p.parse_args()

def main() -> int:
    args = parse_args()
    in_files = read_inputs_list(args.inputs_list)
    genome_names = [os.path.basename(p) for p in in_files]

    ensure_dir(args.outdir)
    ensure_dir(args.baddir)

    infhs = [open(p, "r", encoding="utf-8", newline="") for p in in_files]

    locus_no = 0
    futures = []

    with ProcessPoolExecutor(max_workers=max(1, args.threads)) as pool:
        try:
            while True:
                lines: List[Optional[str]] = []
                any_eof = False
                for fh in infhs:
                    s = fh.readline()
                    if s == "":
                        any_eof = True
                        lines.append(None)
                    else:
                        lines.append(s.rstrip("\n"))

                if any_eof:
                    if not all(x is None for x in lines):
                        raise RuntimeError("Input files have different number of lines (EOF mismatch).")
                    break

                locus_no += 1
                if args.expected_lines and locus_no > args.expected_lines:
                    raise RuntimeError(f"Read more than expected-lines={args.expected_lines}")

                locus_lines = [x for x in lines if x is not None]  # type: ignore

                futures.append(pool.submit(
                    process_one_locus,
                    locus_no,
                    genome_names,
                    locus_lines,
                    args.outdir,
                    args.baddir,
                    args.aligner,
                    args.mafft_bin,
                    args.famsa_bin,
                    args.timeout_sec,
                    args.uppercase,
                    args.strip_spaces
                ))

                if args.progress_every and locus_no % args.progress_every == 0:
                    print(f"[queued] loci: {locus_no}", file=sys.stderr)

            ok = bad = 0
            for fut in as_completed(futures):
                _, status = fut.result()
                if status == "ok":
                    ok += 1
                else:
                    bad += 1

            if args.expected_lines and locus_no != args.expected_lines:
                raise RuntimeError(f"Expected {args.expected_lines} loci, got {locus_no}")

            print(f"[done] aligner={args.aligner} genomes={len(in_files)} loci={locus_no} ok={ok} bad={bad}", file=sys.stderr)
            return 0

        finally:
            for fh in infhs:
                try:
                    fh.close()
                except Exception:
                    pass

if __name__ == "__main__":
    raise SystemExit(main())
