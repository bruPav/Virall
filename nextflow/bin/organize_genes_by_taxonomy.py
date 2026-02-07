#!/usr/bin/env python3
"""
Organize predicted genes by Kaiju taxonomy.

Creates a directory structure:
  by_taxonomy/
  ├── gene_taxonomy_mapping.tsv
  ├── <Class1>/
  │   ├── <Family>_<Genus>.faa/.fna
  │   ├── <Family>_unclassified.faa/.fna
  │   └── unclassified.faa/.fna
  ├── <Class2>/
  │   └── ...
  └── unclassified/
      └── all.faa/.fna
"""

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

# Try to import BioPython for FASTA parsing
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False


def parse_prodigal_header(header):
    """
    Parse Prodigal FASTA header to extract contig_id and coordinates.
    
    Example header:
    >NODE_193_length_5055_cov_8896.076325_1 # 947 # 1807 # -1 # ID=1_1;partial=00;...
    
    Returns: (gene_id, contig_id, start, end, strand)
    """
    parts = header.split(" # ")
    gene_id = parts[0].lstrip(">")
    
    # Extract contig_id by removing the last _N (gene number)
    # Gene IDs are like: NODE_193_length_5055_cov_8896.076325_1
    match = re.match(r"(.+)_(\d+)$", gene_id)
    if match:
        contig_id = match.group(1)
    else:
        contig_id = gene_id
    
    start = int(parts[1]) if len(parts) > 1 else None
    end = int(parts[2]) if len(parts) > 2 else None
    strand = int(parts[3]) if len(parts) > 3 else 1
    
    return gene_id, contig_id, start, end, strand


def parse_kaiju_results(kaiju_file):
    """
    Parse Kaiju results file to get contig_id -> taxonomy mapping.
    
    Returns dict: contig_id -> {
        'taxon_id': str,
        'lineage': str,
        'parsed_lineage': list of taxonomy levels
    }
    """
    contig_taxonomy = {}
    
    with open(kaiju_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split("\t")
            status = parts[0]
            contig_id = parts[1]
            
            if status == "C" and len(parts) >= 4:
                taxon_id = parts[2]
                lineage = parts[3]
                
                # Parse lineage: "Viruses; Nucleocytoviricota; Megaviricetes; ..."
                levels = [p.strip() for p in lineage.split(";") if p.strip()]
                
                contig_taxonomy[contig_id] = {
                    'taxon_id': taxon_id,
                    'lineage': lineage,
                    'parsed_lineage': levels
                }
            else:
                # Unclassified
                contig_taxonomy[contig_id] = {
                    'taxon_id': '0',
                    'lineage': '',
                    'parsed_lineage': []
                }
    
    return contig_taxonomy


def get_taxonomy_key(lineage_levels):
    """
    Generate a taxonomy key for file/directory naming.
    
    Lineage levels typically: [Viruses, Phylum, Class, Order, Family, Genus, Species]
    We use Class (index 2) for directory, and Family_Genus for filename.
    
    Returns: (class_dir, file_basename)
    """
    if not lineage_levels:
        return ("unclassified", "all")
    
    # Skip "Viruses" if present at position 0
    levels = lineage_levels[1:] if lineage_levels[0] == "Viruses" else lineage_levels
    
    # Standard viral taxonomy: Phylum, Class, Order, Family, Genus, Species
    # We want Class for directory (index 1 after removing Viruses)
    # and Family_Genus for filename (indices 3, 4)
    
    # Find class (typically index 1 or 2 depending on structure)
    # For viruses: Phylum (1), Class (2), Order (3), Family (4), Genus (5), Species (6)
    # After removing "Viruses": Phylum (0), Class (1), Order (2), Family (3), Genus (4), Species (5)
    
    class_dir = "unclassified"
    file_parts = []
    
    # Determine class directory (use first non-NA after Viruses, typically Phylum or Class)
    # For better grouping, use index 1 (Class level) if available
    if len(levels) >= 2:
        # Index 1 is typically Class for viruses
        class_candidate = levels[1] if levels[1] != "NA" else levels[0]
        if class_candidate != "NA":
            class_dir = sanitize_filename(class_candidate)
    elif len(levels) >= 1 and levels[0] != "NA":
        class_dir = sanitize_filename(levels[0])
    
    # Build filename from Family + Genus (indices 3, 4 after Viruses removal)
    # Or use lowest non-NA level with _unclassified suffix
    family_idx = 3  # Family position
    genus_idx = 4   # Genus position
    species_idx = 5  # Species position
    
    if len(levels) > family_idx and levels[family_idx] != "NA":
        file_parts.append(levels[family_idx])
        
        if len(levels) > genus_idx and levels[genus_idx] != "NA":
            file_parts.append(levels[genus_idx])
            
            if len(levels) > species_idx and levels[species_idx] != "NA":
                # Full classification - add species
                file_parts.append(levels[species_idx])
            else:
                # Genus known but species = NA
                file_parts.append("unclassified")
        else:
            # Family known but genus = NA
            file_parts.append("unclassified")
    elif len(levels) > 0:
        # No family, find lowest non-NA level
        lowest_non_na = None
        for level in reversed(levels):
            if level != "NA":
                lowest_non_na = level
                break
        
        if lowest_non_na:
            file_parts.append(lowest_non_na)
            file_parts.append("unclassified")
        else:
            file_parts.append("unclassified")
    else:
        file_parts.append("all")
    
    file_basename = "_".join(sanitize_filename(p) for p in file_parts)
    
    return (class_dir, file_basename)


def sanitize_filename(name):
    """Sanitize a string to be safe for use as filename."""
    # Replace problematic characters
    name = re.sub(r'[<>:"/\\|?*]', '_', name)
    name = re.sub(r'\s+', '_', name)
    name = re.sub(r'_+', '_', name)
    name = name.strip('_')
    return name if name else "unknown"


def reverse_complement(seq):
    """Return reverse complement of DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                  'N': 'N', 'n': 'n'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def extract_gene_nucleotide(contig_seq, start, end, strand):
    """Extract gene nucleotide sequence from contig (1-based coordinates)."""
    # Prodigal uses 1-based coordinates
    gene_seq = contig_seq[start-1:end]
    
    if strand == -1:
        gene_seq = reverse_complement(gene_seq)
    
    return gene_seq


def simple_fasta_parser(fasta_file):
    """Simple FASTA parser when BioPython not available."""
    records = {}
    current_id = None
    current_seq = []
    
    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    records[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_id:
            records[current_id] = ''.join(current_seq)
    
    return records


def main():
    parser = argparse.ArgumentParser(
        description="Organize predicted genes by Kaiju taxonomy"
    )
    parser.add_argument("--proteins", required=True, 
                        help="Input proteins FASTA (from Prodigal)")
    parser.add_argument("--contigs", required=True,
                        help="Input viral contigs FASTA (nucleotide)")
    parser.add_argument("--kaiju", required=True,
                        help="Kaiju results with names TSV")
    parser.add_argument("--outdir", required=True,
                        help="Output directory for organized files")
    
    args = parser.parse_args()
    
    proteins_file = Path(args.proteins)
    contigs_file = Path(args.contigs)
    kaiju_file = Path(args.kaiju)
    out_dir = Path(args.outdir)
    
    # Validate inputs
    for f in [proteins_file, contigs_file, kaiju_file]:
        if not f.exists():
            print(f"Error: {f} not found", file=sys.stderr)
            sys.exit(1)
    
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"[organize_genes] Loading Kaiju taxonomy from {kaiju_file}", file=sys.stderr)
    contig_taxonomy = parse_kaiju_results(kaiju_file)
    print(f"[organize_genes] Loaded taxonomy for {len(contig_taxonomy)} contigs", file=sys.stderr)
    
    # Load contig sequences
    print(f"[organize_genes] Loading contig sequences from {contigs_file}", file=sys.stderr)
    if HAS_BIOPYTHON:
        contig_seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(contigs_file, "fasta")}
    else:
        contig_seqs = simple_fasta_parser(contigs_file)
    print(f"[organize_genes] Loaded {len(contig_seqs)} contig sequences", file=sys.stderr)
    
    # Parse proteins and organize by taxonomy
    print(f"[organize_genes] Processing proteins from {proteins_file}", file=sys.stderr)
    
    # Structure: {(class_dir, file_basename): [(gene_id, aa_seq, nt_seq, full_header, taxonomy_info), ...]}
    genes_by_taxon = defaultdict(list)
    
    # Also collect mapping data for TSV output
    mapping_data = []
    
    # Read proteins file
    if HAS_BIOPYTHON:
        for record in SeqIO.parse(proteins_file, "fasta"):
            header = record.description
            aa_seq = str(record.seq)
            
            gene_id, contig_id, start, end, strand = parse_prodigal_header(header)
            
            # Get taxonomy for this contig
            tax_info = contig_taxonomy.get(contig_id, {
                'taxon_id': '0',
                'lineage': '',
                'parsed_lineage': []
            })
            
            # Get taxonomy key for organizing
            class_dir, file_basename = get_taxonomy_key(tax_info['parsed_lineage'])
            
            # Extract nucleotide sequence
            nt_seq = ""
            if contig_id in contig_seqs and start and end:
                nt_seq = extract_gene_nucleotide(contig_seqs[contig_id], start, end, strand)
            
            genes_by_taxon[(class_dir, file_basename)].append({
                'gene_id': gene_id,
                'header': header,
                'aa_seq': aa_seq,
                'nt_seq': nt_seq,
                'contig_id': contig_id,
                'start': start,
                'end': end,
                'strand': strand,
                'lineage': tax_info['lineage']
            })
            
            mapping_data.append({
                'gene_id': gene_id,
                'contig_id': contig_id,
                'start': start,
                'end': end,
                'strand': strand,
                'class_dir': class_dir,
                'file_basename': file_basename,
                'lineage': tax_info['lineage']
            })
    else:
        # Simple parsing without BioPython
        current_header = None
        current_seq = []
        
        with open(proteins_file) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_header:
                        aa_seq = ''.join(current_seq)
                        gene_id, contig_id, start, end, strand = parse_prodigal_header(current_header)
                        
                        tax_info = contig_taxonomy.get(contig_id, {
                            'taxon_id': '0',
                            'lineage': '',
                            'parsed_lineage': []
                        })
                        
                        class_dir, file_basename = get_taxonomy_key(tax_info['parsed_lineage'])
                        
                        nt_seq = ""
                        if contig_id in contig_seqs and start and end:
                            nt_seq = extract_gene_nucleotide(contig_seqs[contig_id], start, end, strand)
                        
                        genes_by_taxon[(class_dir, file_basename)].append({
                            'gene_id': gene_id,
                            'header': current_header,
                            'aa_seq': aa_seq,
                            'nt_seq': nt_seq,
                            'contig_id': contig_id,
                            'start': start,
                            'end': end,
                            'strand': strand,
                            'lineage': tax_info['lineage']
                        })
                        
                        mapping_data.append({
                            'gene_id': gene_id,
                            'contig_id': contig_id,
                            'start': start,
                            'end': end,
                            'strand': strand,
                            'class_dir': class_dir,
                            'file_basename': file_basename,
                            'lineage': tax_info['lineage']
                        })
                    
                    current_header = line[1:]
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # Don't forget last record
            if current_header:
                aa_seq = ''.join(current_seq)
                gene_id, contig_id, start, end, strand = parse_prodigal_header(current_header)
                
                tax_info = contig_taxonomy.get(contig_id, {
                    'taxon_id': '0',
                    'lineage': '',
                    'parsed_lineage': []
                })
                
                class_dir, file_basename = get_taxonomy_key(tax_info['parsed_lineage'])
                
                nt_seq = ""
                if contig_id in contig_seqs and start and end:
                    nt_seq = extract_gene_nucleotide(contig_seqs[contig_id], start, end, strand)
                
                genes_by_taxon[(class_dir, file_basename)].append({
                    'gene_id': gene_id,
                    'header': current_header,
                    'aa_seq': aa_seq,
                    'nt_seq': nt_seq,
                    'contig_id': contig_id,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'lineage': tax_info['lineage']
                })
                
                mapping_data.append({
                    'gene_id': gene_id,
                    'contig_id': contig_id,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'class_dir': class_dir,
                    'file_basename': file_basename,
                    'lineage': tax_info['lineage']
                })
    
    total_genes = len(mapping_data)
    print(f"[organize_genes] Processed {total_genes} genes", file=sys.stderr)
    
    # Write organized files
    files_written = 0
    for (class_dir, file_basename), genes in genes_by_taxon.items():
        # Create class directory
        class_path = out_dir / class_dir
        class_path.mkdir(parents=True, exist_ok=True)
        
        # Write protein FASTA
        faa_file = class_path / f"{file_basename}.faa"
        with open(faa_file, 'w') as f:
            for gene in genes:
                f.write(f">{gene['header']}\n")
                # Wrap sequence at 60 characters
                seq = gene['aa_seq']
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i+60] + "\n")
        
        # Write nucleotide FASTA (only if we have sequences)
        genes_with_nt = [g for g in genes if g['nt_seq']]
        if genes_with_nt:
            fna_file = class_path / f"{file_basename}.fna"
            with open(fna_file, 'w') as f:
                for gene in genes_with_nt:
                    f.write(f">{gene['gene_id']} {gene['contig_id']}:{gene['start']}-{gene['end']}({'+' if gene['strand'] == 1 else '-'})\n")
                    seq = gene['nt_seq']
                    for i in range(0, len(seq), 60):
                        f.write(seq[i:i+60] + "\n")
            files_written += 2
        else:
            files_written += 1
    
    # Write complete genes.fna (all genes nucleotide)
    all_fna_file = out_dir / "genes.fna"
    with open(all_fna_file, 'w') as f:
        for (class_dir, file_basename), genes in genes_by_taxon.items():
            for gene in genes:
                if gene['nt_seq']:
                    f.write(f">{gene['gene_id']} {gene['contig_id']}:{gene['start']}-{gene['end']}({'+' if gene['strand'] == 1 else '-'})\n")
                    seq = gene['nt_seq']
                    for i in range(0, len(seq), 60):
                        f.write(seq[i:i+60] + "\n")
    
    # Write gene taxonomy mapping TSV
    mapping_file = out_dir / "gene_taxonomy_mapping.tsv"
    with open(mapping_file, 'w') as f:
        f.write("gene_id\tcontig_id\tstart\tend\tstrand\tclass_dir\tfile_basename\tlineage\n")
        for row in mapping_data:
            strand_str = "+" if row['strand'] == 1 else "-"
            f.write(f"{row['gene_id']}\t{row['contig_id']}\t{row['start']}\t{row['end']}\t{strand_str}\t{row['class_dir']}\t{row['file_basename']}\t{row['lineage']}\n")
    
    # Print summary
    print(f"[organize_genes] Created {len(genes_by_taxon)} taxonomy groups", file=sys.stderr)
    print(f"[organize_genes] Output directory: {out_dir}", file=sys.stderr)
    
    # Summary by class
    class_summary = defaultdict(int)
    for (class_dir, _), genes in genes_by_taxon.items():
        class_summary[class_dir] += len(genes)
    
    print(f"[organize_genes] Genes by class:", file=sys.stderr)
    for class_name, count in sorted(class_summary.items(), key=lambda x: -x[1]):
        print(f"  {class_name}: {count} genes", file=sys.stderr)


if __name__ == "__main__":
    main()
