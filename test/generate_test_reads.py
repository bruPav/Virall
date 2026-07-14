#!/usr/bin/env python3
"""Generate synthetic paired-end reads from a reference FASTA for pipeline testing.
No external dependencies — stdlib only. Deterministic output (fixed random seed)."""

import gzip
import random
import sys
from pathlib import Path

SEED = 42
READ_LEN = 150
INSERT_SIZE = 300
NUM_PAIRS = 200
PHRED = 40  # perfect quality


def rc(seq: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def main():
    ref_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("test/data/phiX174.fasta")
    out_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("test/data")

    # read reference
    with open(ref_path) as f:
        lines = [l.strip() for l in f if not l.startswith(">")]
    genome = "".join(lines)

    random.seed(SEED)
    glen = len(genome)

    r1_lines = []
    r2_lines = []

    for i in range(NUM_PAIRS):
        start = random.randint(0, glen - INSERT_SIZE - 1)
        end = start + INSERT_SIZE
        frag = genome[start:end]
        read1 = frag[:READ_LEN]
        read2 = rc(frag[-READ_LEN:])
        qual = chr(PHRED + 33) * READ_LEN
        r1_lines.append(f"@read_{i}/1\n{read1}\n+\n{qual}\n")
        r2_lines.append(f"@read_{i}/2\n{read2}\n+\n{qual}\n")

    out_dir.mkdir(parents=True, exist_ok=True)
    r1_path = out_dir / "test_R1.fastq.gz"
    r2_path = out_dir / "test_R2.fastq.gz"

    with gzip.open(r1_path, "wt") as f:
        f.writelines(r1_lines)
    with gzip.open(r2_path, "wt") as f:
        f.writelines(r2_lines)

    print(f"Generated {NUM_PAIRS} paired-end reads ({READ_LEN} bp each) from {ref_path.name}")
    print(f"  {r1_path} ({r1_path.stat().st_size} bytes)")
    print(f"  {r2_path} ({r2_path.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
