# Dna-repeating-subsequence-
Python prototype for maximum-oriented detection of repeating subsequences in DNA strings
# Maximum-Oriented Detection of Repeating Subsequences in DNA Strings

Python prototype for the multi-start framework described in:

> "Maximum-Oriented Detection of Repeating Subsequences in DNA Strings"  

Implements **Algorithm 1** (MaxOrientedRepSubseq) with a greedy extension
heuristic in place of the exact ExtendkRepSubseq subroutine, as described
in Sections V-B and VII-B of the paper. Supports FASTA input.

## Requirements

Python 3.8 or later. No third-party libraries required.

## Usage

**Inline sequence:**
```bash
python dna_repeating_subseq.py --seq ATGCATGCATGC --k 2
```

**FASTA file (first record):**
```bash
python dna_repeating_subseq.py --fasta example.fa --k 2
```

**FASTA file (all records):**
```bash
python dna_repeating_subseq.py --fasta example.fa --k 2 --all
```

**Show the best subsequence found:**
```bash
python dna_repeating_subseq.py --fasta example.fa --k 2 --show-seq
```

## Options

| Flag | Description |
|------|-------------|
| `--seq DNA` | Inline DNA string (A/T/G/C only) |
| `--fasta FILE` | Path to FASTA file |
| `--k INT` | Repetition factor k ≥ 2 (default: 2) |
| `--all` | Process all records in FASTA (default: first only) |
| `--show-seq` | Print the best subsequence (truncated at 200 chars) |

## Example output
