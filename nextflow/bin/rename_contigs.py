#!/usr/bin/env python3
"""
Rename viral contigs based on their lowest Kaiju taxonomic classification.

Replaces assembler-generated names (e.g. NODE_1_length_5000_cov_100 or
contig_1) with informative taxonomy-based names (e.g. Escherichia_phage_T4_1).

Outputs:
  - Renamed FASTA file with original name preserved in header description
  - Updated Kaiju results TSV files with new contig IDs
  - Name mapping file (original -> new) for traceability
"""

import argparse
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path


def sanitize_name(name):
    """Sanitize a taxonomy name for use as a contig identifier."""
    name = re.sub(r'[<>:"/\\|?*\[\]{}()\'`~!@#$%^&+=,;]', '', name)
    name = re.sub(r'\s+', '_', name)
    name = re.sub(r'_+', '_', name)
    name = name.strip('_.')
    return name if name else "unknown"


def get_lowest_classification(lineage_str):
    """
    Extract the lowest (most specific) non-empty classification level
    from a Kaiju semicolon-separated lineage string.

    Example: "Viruses; Uroviricota; Caudoviricetes; ...; Escherichia phage T4;"
    Returns: "Escherichia_phage_T4"
    """
    if not lineage_str or not lineage_str.strip():
        return "unclassified_virus"

    levels = [l.strip().rstrip(';') for l in lineage_str.split(';')]
    levels = [l for l in levels if l and l != 'NA' and l != 'Viruses']

    if not levels:
        return "unclassified_virus"

    return sanitize_name(levels[-1])


def parse_kaiju_taxonomy(kaiju_file):
    """Parse Kaiju results-with-names to get contig_id -> lineage mapping."""
    taxonomy = {}

    if not kaiju_file.exists() or kaiju_file.stat().st_size == 0:
        return taxonomy

    with open(kaiju_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 4 and parts[0] == 'C':
                contig_id = parts[1].strip()
                lineage = parts[3].strip()
                taxonomy[contig_id] = lineage

    return taxonomy


def parse_fasta_ids(fasta_file):
    """Get ordered list of contig IDs from a FASTA file."""
    ids = []
    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                cid = line[1:].strip().split()[0]
                ids.append(cid)
    return ids


def build_name_mapping(contig_ids, taxonomy):
    """
    Build old_name -> new_name mapping based on Kaiju taxonomy.

    Each contig gets a name based on its lowest taxonomic classification,
    with a numeric suffix for uniqueness (e.g. Siphoviridae_1, Siphoviridae_2).
    """
    base_names = {}
    for cid in contig_ids:
        lineage = taxonomy.get(cid, '')
        base_names[cid] = get_lowest_classification(lineage)

    name_counts = defaultdict(int)
    mapping = {}

    for cid in contig_ids:
        base = base_names[cid]
        name_counts[base] += 1
        mapping[cid] = f"{base}_{name_counts[base]}"

    return mapping


def rename_fasta(input_fasta, output_fasta, mapping):
    """Write FASTA file with renamed contig headers."""
    with open(input_fasta) as fin, open(output_fasta, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                header = line[1:].strip()
                old_id = header.split()[0]
                rest = header[len(old_id):]
                new_id = mapping.get(old_id, old_id)
                fout.write(f">{new_id} original_name={old_id}{rest}\n")
            else:
                fout.write(line)


def rename_kaiju_tsv(input_tsv, output_tsv, mapping):
    """Write Kaiju TSV with renamed contig IDs (column 2)."""
    with open(input_tsv) as fin, open(output_tsv, 'w') as fout:
        for line in fin:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2 and parts[1].strip() in mapping:
                parts[1] = mapping[parts[1].strip()]
            fout.write('\t'.join(parts) + '\n')


def write_mapping(mapping, output_file):
    """Write name mapping TSV for traceability."""
    with open(output_file, 'w') as f:
        f.write("original_name\tnew_name\n")
        for old, new in mapping.items():
            f.write(f"{old}\t{new}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Rename viral contigs using Kaiju taxonomy classification"
    )
    parser.add_argument("--fasta", type=Path, required=True,
                        help="Input viral contigs FASTA")
    parser.add_argument("--kaiju-dir", type=Path, required=True,
                        help="Input Kaiju results directory")
    parser.add_argument("--out-fasta", type=Path, required=True,
                        help="Output renamed FASTA")
    parser.add_argument("--out-kaiju-dir", type=Path, required=True,
                        help="Output directory for updated Kaiju results")
    parser.add_argument("--out-mapping", type=Path, required=True,
                        help="Output name mapping TSV")

    args = parser.parse_args()

    kaiju_names_file = args.kaiju_dir / "kaiju_results_with_names.tsv"
    kaiju_raw_file = args.kaiju_dir / "kaiju_results.tsv"

    # Handle empty input gracefully
    if not args.fasta.exists() or args.fasta.stat().st_size == 0:
        print("[rename_contigs] No contigs to rename (empty FASTA)", file=sys.stderr)
        args.out_fasta.touch()
        args.out_kaiju_dir.mkdir(parents=True, exist_ok=True)
        if kaiju_names_file.exists():
            shutil.copy2(kaiju_names_file,
                         args.out_kaiju_dir / "kaiju_results_with_names.tsv")
        else:
            (args.out_kaiju_dir / "kaiju_results_with_names.tsv").touch()
        if kaiju_raw_file.exists():
            shutil.copy2(kaiju_raw_file,
                         args.out_kaiju_dir / "kaiju_results.tsv")
        with open(args.out_mapping, 'w') as f:
            f.write("original_name\tnew_name\n")
        return

    # Parse inputs
    print(f"[rename_contigs] Reading taxonomy from {kaiju_names_file}",
          file=sys.stderr)
    taxonomy = parse_kaiju_taxonomy(kaiju_names_file)
    contig_ids = parse_fasta_ids(args.fasta)

    print(f"[rename_contigs] Found {len(contig_ids)} contigs in FASTA",
          file=sys.stderr)
    print(f"[rename_contigs] Found taxonomy for {len(taxonomy)} contigs "
          f"in Kaiju results", file=sys.stderr)

    # Build name mapping
    mapping = build_name_mapping(contig_ids, taxonomy)

    # Write renamed FASTA
    rename_fasta(args.fasta, args.out_fasta, mapping)

    # Write updated Kaiju results
    args.out_kaiju_dir.mkdir(parents=True, exist_ok=True)
    if kaiju_names_file.exists():
        rename_kaiju_tsv(
            kaiju_names_file,
            args.out_kaiju_dir / "kaiju_results_with_names.tsv",
            mapping
        )
    else:
        (args.out_kaiju_dir / "kaiju_results_with_names.tsv").touch()

    if kaiju_raw_file.exists():
        rename_kaiju_tsv(
            kaiju_raw_file,
            args.out_kaiju_dir / "kaiju_results.tsv",
            mapping
        )

    # Write mapping file
    write_mapping(mapping, args.out_mapping)

    # Print summary
    print(f"[rename_contigs] Renamed {len(mapping)} contigs:", file=sys.stderr)
    for old, new in list(mapping.items())[:10]:
        print(f"  {old} -> {new}", file=sys.stderr)
    if len(mapping) > 10:
        print(f"  ... and {len(mapping) - 10} more", file=sys.stderr)


if __name__ == "__main__":
    main()
