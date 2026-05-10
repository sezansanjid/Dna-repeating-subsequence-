"""
Maximum-Oriented Detection of Repeating Subsequences in DNA Strings
====================================================================
Python prototype implementing Algorithm 1 (MaxOrientedRepSubseq) from:

  "Maximum-Oriented Detection of Repeating Subsequences in DNA Strings"

The ExtendkRepSubseq call in Algorithm 1 is replaced here by a greedy
extension heuristic, exactly as described in Section V-B and Section VII-B
of the paper.  The theoretical guarantees of Section VI apply to the exact
subroutine; this prototype is the research vehicle used for the experiments
in Section VIII.

Usage
-----
  # Inline sequence
  python dna_repeating_subseq.py --seq ATGCATGCATGC --k 2

  # FASTA file (first record is used)
  python dna_repeating_subseq.py --fasta sequence.fa --k 2

  # FASTA file, all records
  python dna_repeating_subseq.py --fasta sequence.fa --k 2 --all

Requirements: Python 3.8+  (no third-party libraries)
"""

import argparse
import sys
import time

# ---------------------------------------------------------------------------
# FASTA I/O
# ---------------------------------------------------------------------------

def parse_fasta(path: str) -> list[tuple[str, str]]:
    """Parse a FASTA file; return list of (header, sequence) pairs."""
    records: list[tuple[str, str]] = []
    header = ""
    parts: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    records.append((header, "".join(parts).upper()))
                header = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if header:
        records.append((header, "".join(parts).upper()))
    return records


# ---------------------------------------------------------------------------
# Core: k-repeating subsequence check
# ---------------------------------------------------------------------------

def is_k_repeating(S: str, X: str, k: int) -> bool:
    """
    Return True iff X^k is a subsequence of S.

    Greedily match k consecutive copies of X inside S.
    Time: O(k * |X|) per call (linear in the combined length).
    """
    pos = 0
    n = len(S)
    for _ in range(k):
        for ch in X:
            while pos < n and S[pos] != ch:
                pos += 1
            if pos >= n:
                return False
            pos += 1
    return True


# ---------------------------------------------------------------------------
# Greedy extension heuristic  (replaces ExtendkRepSubseq in the prototype)
# ---------------------------------------------------------------------------

def greedy_extend(S: str, X: str, k: int) -> str:
    """
    Greedy extension heuristic described in Section VII-B of the paper.

    Iterates over positions in S from left to right; at each position
    attempts to insert S[pos] into the current candidate at every possible
    insertion point.  The first insertion point that maintains X^k ⊆ S is
    accepted, and the scan continues from the next position.

    This is a heuristic — it does not implement the exact ExtendkRepSubseq
    algorithm of Gong et al. (Algorithm 2 in the paper), but it is the
    procedure used in the reported experiments.
    """
    for pos in range(len(S)):
        ch = S[pos]
        # Try inserting ch at every position in X
        for insert_at in range(len(X) + 1):
            candidate = X[:insert_at] + ch + X[insert_at:]
            if is_k_repeating(S, candidate, k):
                X = candidate
                break
    return X


# ---------------------------------------------------------------------------
# Algorithm 1: MaxOrientedRepSubseq
# ---------------------------------------------------------------------------

DNA_ALPHABET = ("A", "T", "G", "C")


def max_oriented_rep_subseq(S: str, k: int = 2) -> dict:
    """
    Algorithm 1 — MaxOrientedRepSubseq(S, k).

    For each nucleotide σ ∈ {A, T, G, C}:
      1. Count occurrences r = floor(count(σ) / k).
         If r = 0, skip (σ appears fewer than k times).
      2. Form initial candidate X = σ^r  (trivially satisfies X^k ⊆ S).
      3. Run greedy_extend(S, X, k).
    Return the longest result together with the winning nucleotide.

    Returns a dict with keys:
      best_seq       – the best subsequence found (str)
      best_len       – its length (int)
      best_sigma     – the nucleotide that yielded the best result (str)
      per_sigma      – {σ: length} for every σ that was tried (dict)
    """
    best_seq = ""
    best_sigma = ""
    per_sigma: dict[str, int] = {}

    for sigma in DNA_ALPHABET:
        r = S.count(sigma) // k
        if r == 0:
            continue                          # Line 4 of Algorithm 1: skip
        X_init = sigma * r                    # initial k-repeating candidate
        X = greedy_extend(S, X_init, k)      # Algorithm 2 placeholder
        per_sigma[sigma] = len(X)
        if len(X) > len(best_seq):
            best_seq = X
            best_sigma = sigma

    return {
        "best_seq": best_seq,
        "best_len": len(best_seq),
        "best_sigma": best_sigma,
        "per_sigma": per_sigma,
    }


# ---------------------------------------------------------------------------
# Single-start baseline  (initialises from A only — Section VII-B)
# ---------------------------------------------------------------------------

def single_start_baseline(S: str, k: int = 2, sigma: str = "A") -> dict:
    """
    Single-start baseline: runs one extension call from the given nucleotide.
    Defaults to σ = A, matching the baseline described in the paper.
    """
    r = S.count(sigma) // k
    if r == 0:
        return {"best_seq": "", "best_len": 0, "best_sigma": sigma}
    X_init = sigma * r
    X = greedy_extend(S, X_init, k)
    return {"best_seq": X, "best_len": len(X), "best_sigma": sigma}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def validate_dna(seq: str, label: str = "") -> str:
    """Strip whitespace, upper-case, and verify only {A,T,G,C}."""
    seq = seq.replace(" ", "").replace("\n", "").upper()
    bad = set(seq) - set("ATGC")
    if bad:
        sys.exit(
            f"Error{' (' + label + ')' if label else ''}: "
            f"Non-DNA characters found: {bad}"
        )
    return seq


def run_on_sequence(header: str, S: str, k: int, show_seq: bool) -> None:
    """Run both baseline and multi-start; print a formatted report."""
    n = len(S)
    print(f"\n{'=' * 60}")
    print(f"Sequence : {header}")
    print(f"Length   : {n}")
    print(f"k        : {k}")
    print("-" * 60)

    t0 = time.perf_counter()
    baseline = single_start_baseline(S, k, sigma="A")
    t_base = (time.perf_counter() - t0) * 1000

    t0 = time.perf_counter()
    result = max_oriented_rep_subseq(S, k)
    t_multi = (time.perf_counter() - t0) * 1000

    b_len = baseline["best_len"]
    m_len = result["best_len"]
    delta = (
        100.0 * (m_len - b_len) / b_len if b_len > 0 else float("nan")
    )

    print(f"Baseline (σ=A)     : length = {b_len}  ({t_base:.2f} ms)")
    print(f"Multi-start        : length = {m_len}  ({t_multi:.2f} ms)")
    print(f"Best nucleotide    : {result['best_sigma']}")
    print(f"Improvement        : {delta:+.2f}%")
    print(f"Per-nucleotide     : {result['per_sigma']}")
    if show_seq and m_len <= 200:
        print(f"Best subsequence   : {result['best_seq']}")
    elif show_seq:
        print(f"Best subsequence   : {result['best_seq'][:100]}...  (truncated)")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Maximum-Oriented Repeating Subsequence Detection in DNA "
            "(prototype implementation of Algorithm 1)"
        )
    )
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument(
        "--seq", metavar="DNA",
        help="DNA sequence string (letters A/T/G/C only)"
    )
    src.add_argument(
        "--fasta", metavar="FILE",
        help="Path to a FASTA file"
    )
    parser.add_argument(
        "--k", type=int, default=2,
        help="Repetition factor k ≥ 2 (default: 2 = square subsequence)"
    )
    parser.add_argument(
        "--all", action="store_true",
        help="Process all records in the FASTA file (default: first record only)"
    )
    parser.add_argument(
        "--show-seq", action="store_true",
        help="Print the best subsequence found (truncated at 200 chars)"
    )
    args = parser.parse_args()

    if args.k < 2:
        sys.exit("Error: k must be ≥ 2.")

    if args.seq:
        S = validate_dna(args.seq, label="--seq")
        run_on_sequence("(command-line input)", S, args.k, args.show_seq)
    else:
        records = parse_fasta(args.fasta)
        if not records:
            sys.exit(f"Error: no FASTA records found in {args.fasta}")
        subset = records if args.all else records[:1]
        for header, seq in subset:
            S = validate_dna(seq, label=header)
            run_on_sequence(header, S, args.k, args.show_seq)


if __name__ == "__main__":
    main()
