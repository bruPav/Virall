#!/usr/bin/env python3
"""
Extract FASTA sequences by ID list. Used by Nextflow FILTER_VIRAL.
Reads IDs from a file (one per line) or from Kaiju output (C lines, column 2).
Writes only those sequences from the input FASTA.
"""
import argparse
import sys
from pathlib import Path

try:
    from Bio import SeqIO
except ImportError:
    sys.exit("extract_fasta_by_ids requires biopython: pip install biopython")


def main():
    parser = argparse.ArgumentParser(description="Extract FASTA sequences by ID list")
    parser.add_argument("--ids", type=Path, help="File with one ID per line")
    parser.add_argument("--kaiju", type=Path, help="Kaiju results TSV (use C lines, column 2 as IDs)")
    parser.add_argument("--fasta", type=Path, required=True, help="Input FASTA")
    parser.add_argument("--out", type=Path, required=True, help="Output FASTA")
    parser.add_argument("--min-len", type=int, default=0, help="Minimum contig length (bp) to keep")
    args = parser.parse_args()

    if args.ids:
        with open(args.ids) as f:
            wanted = {line.strip().split()[0] for line in f if line.strip()}
    elif args.kaiju:
        wanted = set()
        with open(args.kaiju) as f:
            for line in f:
                if line.startswith("C"):  # Classified
                    parts = line.strip().split("\t")
                    if len(parts) >= 2:
                        wanted.add(parts[1])
    else:
        sys.exit("Provide --ids or --kaiju")

    count = 0
    skipped_short = 0
    with open(args.out, "w") as out_f:
        for rec in SeqIO.parse(args.fasta, "fasta"):
            # match by full id or id before first space
            seq_id = rec.id.split()[0]
            if seq_id in wanted or rec.id in wanted:
                if args.min_len and len(rec.seq) < args.min_len:
                    skipped_short += 1
                    continue
                SeqIO.write(rec, out_f, "fasta")
                count += 1
    if skipped_short:
        print(f"Skipped {skipped_short} contigs shorter than {args.min_len} bp", file=sys.stderr)
    print(f"Wrote {count} sequences to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
