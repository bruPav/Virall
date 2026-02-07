#!/usr/bin/env python3
"""
Build 10x-compatible sparse matrix from single-cell viral UMI counts.

Outputs:
  - matrix.mtx: Sparse matrix in Matrix Market format
  - barcodes.tsv: Cell barcodes (one per line)
  - features.tsv: Viral contig features (contig_id, contig_name, feature_type)
  - sc_viral_summary.tsv: Per-cell viral load statistics
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path

try:
    from scipy.io import mmwrite
    from scipy.sparse import csr_matrix
    import numpy as np
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


def parse_kaiju_taxonomy(kaiju_dir):
    """Parse Kaiju results to get contig taxonomy."""
    taxonomy = {}
    kaiju_file = Path(kaiju_dir) / "kaiju_results_with_names.tsv"
    
    if not kaiju_file.exists():
        return taxonomy
    
    with open(kaiju_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4 and parts[0] == 'C':
                contig_id = parts[1]
                lineage = parts[3] if len(parts) > 3 else ''
                # Extract species/genus name from lineage
                levels = [l.strip() for l in lineage.split(';') if l.strip() and l.strip() != 'NA']
                name = levels[-1] if levels else contig_id
                taxonomy[contig_id] = name
    
    return taxonomy


def parse_counts(counts_file):
    """
    Parse UMI counts from umi_tools count output or simple count table.
    
    Expected format (tab-separated):
    cell    gene    count
    ACGT... contig1 5
    """
    counts = defaultdict(lambda: defaultdict(int))
    
    with open(counts_file) as f:
        header = f.readline().strip().split('\t')
        
        # Handle different column names
        cell_col = 0
        gene_col = 1
        count_col = 2
        
        for i, col in enumerate(header):
            col_lower = col.lower()
            if col_lower in ('cell', 'cell_barcode', 'barcode', 'cb'):
                cell_col = i
            elif col_lower in ('gene', 'contig', 'feature', 'gene_id'):
                gene_col = i
            elif col_lower in ('count', 'counts', 'umi_count', 'n'):
                count_col = i
        
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) > max(cell_col, gene_col, count_col):
                cell = parts[cell_col]
                gene = parts[gene_col]
                try:
                    count = int(parts[count_col])
                except ValueError:
                    continue
                counts[cell][gene] += count
    
    return counts


def parse_contig_lengths(contigs_file):
    """Parse contig lengths from FASTA file."""
    lengths = {}
    current_id = None
    current_len = 0
    
    with open(contigs_file) as f:
        for line in f:
            if line.startswith('>'):
                if current_id:
                    lengths[current_id] = current_len
                current_id = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line.strip())
        
        if current_id:
            lengths[current_id] = current_len
    
    return lengths


def write_mtx_manual(matrix, row_indices, col_indices, filepath):
    """Write Matrix Market format manually if scipy not available."""
    n_rows = max(row_indices) + 1 if row_indices else 0
    n_cols = max(col_indices) + 1 if col_indices else 0
    n_entries = len(matrix)
    
    with open(filepath, 'w') as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n")
        f.write(f"{n_rows} {n_cols} {n_entries}\n")
        for i, (row, col) in enumerate(zip(row_indices, col_indices)):
            # MTX format is 1-indexed
            f.write(f"{row + 1} {col + 1} {matrix[i]}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Build 10x-compatible MTX matrix from viral UMI counts"
    )
    parser.add_argument("--counts", required=True,
                        help="UMI counts TSV from umi_tools count")
    parser.add_argument("--contigs", required=True,
                        help="Viral contigs FASTA")
    parser.add_argument("--kaiju-dir", required=True,
                        help="Kaiju results directory")
    parser.add_argument("--min-reads", type=int, default=100,
                        help="Minimum reads per cell to include (default: 100)")
    parser.add_argument("--min-viral-umis", type=int, default=1,
                        help="Minimum viral UMIs to call cell infected (default: 1)")
    parser.add_argument("--output-dir", default=".",
                        help="Output directory")
    
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Parse input files
    print("[sc_build_matrix] Parsing counts...", file=sys.stderr)
    counts = parse_counts(args.counts)
    print(f"[sc_build_matrix] Found {len(counts)} cells with counts", file=sys.stderr)
    
    print("[sc_build_matrix] Parsing Kaiju taxonomy...", file=sys.stderr)
    taxonomy = parse_kaiju_taxonomy(args.kaiju_dir)
    print(f"[sc_build_matrix] Found taxonomy for {len(taxonomy)} contigs", file=sys.stderr)
    
    print("[sc_build_matrix] Parsing contig lengths...", file=sys.stderr)
    contig_lengths = parse_contig_lengths(args.contigs)
    print(f"[sc_build_matrix] Found {len(contig_lengths)} contigs", file=sys.stderr)
    
    # Get all unique features (viral contigs)
    all_features = set()
    for cell_counts in counts.values():
        all_features.update(cell_counts.keys())
    
    # Filter to only contigs we have (in case of mapping artifacts)
    features = sorted([f for f in all_features if f in contig_lengths])
    feature_to_idx = {f: i for i, f in enumerate(features)}
    
    # Filter cells by minimum reads
    filtered_cells = {}
    for cell, cell_counts in counts.items():
        total_umis = sum(cell_counts.values())
        if total_umis >= args.min_reads:
            filtered_cells[cell] = cell_counts
    
    print(f"[sc_build_matrix] {len(filtered_cells)} cells pass min_reads filter ({args.min_reads})", 
          file=sys.stderr)
    
    if not filtered_cells or not features:
        print("[sc_build_matrix] No data to output!", file=sys.stderr)
        # Write empty files
        with open(output_dir / "matrix.mtx", 'w') as f:
            f.write("%%MatrixMarket matrix coordinate integer general\n")
            f.write("0 0 0\n")
        with open(output_dir / "barcodes.tsv", 'w') as f:
            pass
        with open(output_dir / "features.tsv", 'w') as f:
            pass
        with open(output_dir / "sc_viral_summary.tsv", 'w') as f:
            f.write("cell_barcode\ttotal_viral_umis\tn_viral_contigs\tinfected\ttop_virus\ttop_virus_umis\n")
        return
    
    # Build sorted cell list
    barcodes = sorted(filtered_cells.keys())
    barcode_to_idx = {b: i for i, b in enumerate(barcodes)}
    
    # Build sparse matrix data
    row_indices = []
    col_indices = []
    data = []
    
    for cell, cell_counts in filtered_cells.items():
        cell_idx = barcode_to_idx[cell]
        for feature, count in cell_counts.items():
            if feature in feature_to_idx:
                row_indices.append(feature_to_idx[feature])
                col_indices.append(cell_idx)
                data.append(count)
    
    # Write matrix.mtx
    print("[sc_build_matrix] Writing matrix.mtx...", file=sys.stderr)
    if HAS_SCIPY:
        sparse_matrix = csr_matrix(
            (data, (row_indices, col_indices)),
            shape=(len(features), len(barcodes))
        )
        mmwrite(str(output_dir / "matrix.mtx"), sparse_matrix)
    else:
        write_mtx_manual(data, row_indices, col_indices, output_dir / "matrix.mtx")
    
    # Write barcodes.tsv
    print("[sc_build_matrix] Writing barcodes.tsv...", file=sys.stderr)
    with open(output_dir / "barcodes.tsv", 'w') as f:
        for barcode in barcodes:
            f.write(f"{barcode}\n")
    
    # Write features.tsv (10x format: gene_id, gene_name, feature_type)
    print("[sc_build_matrix] Writing features.tsv...", file=sys.stderr)
    with open(output_dir / "features.tsv", 'w') as f:
        for feature in features:
            name = taxonomy.get(feature, feature)
            f.write(f"{feature}\t{name}\tViral\n")
    
    # Write summary statistics
    print("[sc_build_matrix] Writing sc_viral_summary.tsv...", file=sys.stderr)
    with open(output_dir / "sc_viral_summary.tsv", 'w') as f:
        f.write("cell_barcode\ttotal_viral_umis\tn_viral_contigs\tinfected\ttop_virus\ttop_virus_umis\n")
        
        for cell in barcodes:
            cell_counts = filtered_cells[cell]
            total_umis = sum(cell_counts.values())
            n_contigs = len(cell_counts)
            infected = "yes" if total_umis >= args.min_viral_umis else "no"
            
            # Find top virus
            if cell_counts:
                top_virus = max(cell_counts.items(), key=lambda x: x[1])
                top_name = taxonomy.get(top_virus[0], top_virus[0])
                top_umis = top_virus[1]
            else:
                top_name = "none"
                top_umis = 0
            
            f.write(f"{cell}\t{total_umis}\t{n_contigs}\t{infected}\t{top_name}\t{top_umis}\n")
    
    # Print summary
    n_infected = sum(1 for c in barcodes if sum(filtered_cells[c].values()) >= args.min_viral_umis)
    print(f"[sc_build_matrix] Output summary:", file=sys.stderr)
    print(f"  - {len(barcodes)} cells", file=sys.stderr)
    print(f"  - {len(features)} viral features", file=sys.stderr)
    print(f"  - {len(data)} non-zero entries", file=sys.stderr)
    print(f"  - {n_infected} infected cells (>= {args.min_viral_umis} viral UMIs)", file=sys.stderr)


if __name__ == "__main__":
    main()
