#!/usr/bin/env python3
"""
Organize predicted genes by Kaiju taxonomy and annotate with VOG functions.

Creates a directory structure:
  by_taxonomy/
  ├── gene_annotation.tsv          (comprehensive per-gene table)
  ├── gene_taxonomy_mapping.tsv    (symlink, backwards compat)
  ├── functional_summary.tsv       (counts per VOG functional category)
  ├── <Class1>/
  │   ├── <Family>_<Genus>.faa/.fna
  │   ├── <Family>_unclassified.faa/.fna
  │   └── unclassified.faa/.fna
  ├── <Class2>/
  │   └── ...
  └── unclassified/
      └── all.faa/.fna

When VOG metadata files are provided (--vog-domains, --vog-annotations, etc.),
each gene is annotated with its best VOG hit, functional category, description,
and virus-specificity flags.
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


# ---------------------------------------------------------------------------
# VOG functional annotation helpers
# ---------------------------------------------------------------------------

def parse_hmmscan_domtblout(domtblout_file, evalue_threshold=1e-5):
    """
    Parse HMMER domtblout file and return the best VOG hit per query gene.

    domtblout columns (space-delimited, comment lines start with #):
      0  target name (VOG ID, e.g. VOG00123)
      1  target accession
      2  tlen
      3  query name (gene_id from Prodigal)
      4  query accession
      5  qlen
      6  full-seq E-value
      7  full-seq score
      8  full-seq bias
      ...  (domain-level fields follow)

    Returns: dict gene_id -> {'vog_id': str, 'evalue': float, 'score': float}
    """
    best_hits = {}  # gene_id -> {vog_id, evalue, score}

    with open(domtblout_file) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 8:
                continue

            vog_id = fields[0]
            gene_id = fields[3]
            try:
                evalue = float(fields[6])
                score = float(fields[7])
            except ValueError:
                continue

            if evalue > evalue_threshold:
                continue

            prev = best_hits.get(gene_id)
            if prev is None or evalue < prev['evalue']:
                best_hits[gene_id] = {
                    'vog_id': vog_id,
                    'evalue': evalue,
                    'score': score,
                }

    return best_hits


def load_vog_annotations(annotations_file):
    """
    Load vog.annotations.tsv → dict VOG_ID -> (functional_category, description).

    File format (tab-separated, no header):
      GroupName  ProteinCount  SpeciesCount  FunctionalCategory  ConsensusFunctionalDescription
    """
    annotations = {}
    with open(annotations_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            vog_id = parts[0]
            func_cat = parts[3]
            description = parts[4]
            annotations[vog_id] = {
                'functional_category': func_cat,
                'description': description,
            }
    return annotations


def load_vog_virusonly(virusonly_file):
    """
    Load vog.virusonly.tsv → dict VOG_ID -> (high, medium, low) booleans.

    File format (tab-separated, no header):
      GroupName  High  Medium  Low   (1=True, 0=False)
    """
    virusonly = {}
    with open(virusonly_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            vog_id = parts[0]
            try:
                virusonly[vog_id] = {
                    'high': bool(int(parts[1])),
                    'medium': bool(int(parts[2])),
                    'low': bool(int(parts[3])),
                }
            except (ValueError, IndexError):
                continue
    return virusonly


def load_functional_categories(categories_file):
    """
    Load vogdb.functional_categories.txt → dict letter_code -> category_name.

    File format (one per line):
      Xb beneficial for host
      Xh harmful to host
      Xp involved in viral DNA processing
      ...
    """
    categories = {}
    with open(categories_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # First token is the code, rest is the name
            parts = line.split(None, 1)
            if len(parts) == 2:
                categories[parts[0]] = parts[1]
    return categories


def resolve_category_name(func_cat_str, categories_lookup):
    """
    Resolve a functional category string (may contain multiple letters)
    into a human-readable name.  E.g. 'Xr' -> 'DNA, RNA and nucleotide metabolism'.
    Multiple categories are joined with '; '.
    """
    if not func_cat_str or not categories_lookup:
        return ""
    names = []
    for code in func_cat_str.split(";"):
        code = code.strip()
        if code in categories_lookup:
            names.append(categories_lookup[code])
    return "; ".join(names) if names else func_cat_str


def main():
    parser = argparse.ArgumentParser(
        description="Organize predicted genes by Kaiju taxonomy and annotate with VOG functions"
    )
    parser.add_argument("--proteins", required=True, 
                        help="Input proteins FASTA (from Prodigal)")
    parser.add_argument("--contigs", required=True,
                        help="Input viral contigs FASTA (nucleotide)")
    parser.add_argument("--kaiju", required=True,
                        help="Kaiju results with names TSV")
    parser.add_argument("--outdir", required=True,
                        help="Output directory for organized files")
    # Optional VOG functional annotation inputs
    parser.add_argument("--vog-domains", default=None,
                        help="hmmscan domtblout file (vog_domains.txt)")
    parser.add_argument("--vog-annotations", default=None,
                        help="VOG annotations TSV (vog.annotations.tsv)")
    parser.add_argument("--vog-virusonly", default=None,
                        help="VOG virus-only TSV (vog.virusonly.tsv)")
    parser.add_argument("--vog-categories", default=None,
                        help="VOG functional categories file (vogdb.functional_categories.txt)")
    parser.add_argument("--vog-evalue", type=float, default=1e-5,
                        help="E-value threshold for VOG hits (default: 1e-5)")
    
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
    
    # -----------------------------------------------------------------------
    # Load VOG functional annotation data (if provided)
    # -----------------------------------------------------------------------
    vog_hits = {}         # gene_id -> {vog_id, evalue, score}
    vog_annotations = {}  # vog_id  -> {functional_category, description}
    vog_virusonly = {}    # vog_id  -> {high, medium, low}
    vog_categories = {}   # letter_code -> category_name
    has_vog = False

    if args.vog_domains and Path(args.vog_domains).exists() and Path(args.vog_domains).stat().st_size > 0:
        print(f"[organize_genes] Loading VOG domain hits from {args.vog_domains}", file=sys.stderr)
        vog_hits = parse_hmmscan_domtblout(args.vog_domains, evalue_threshold=args.vog_evalue)
        print(f"[organize_genes] Loaded best VOG hit for {len(vog_hits)} genes (e-value <= {args.vog_evalue})", file=sys.stderr)
        has_vog = True

    if args.vog_annotations and Path(args.vog_annotations).exists():
        print(f"[organize_genes] Loading VOG annotations from {args.vog_annotations}", file=sys.stderr)
        vog_annotations = load_vog_annotations(args.vog_annotations)
        print(f"[organize_genes] Loaded annotations for {len(vog_annotations)} VOGs", file=sys.stderr)

    if args.vog_virusonly and Path(args.vog_virusonly).exists():
        print(f"[organize_genes] Loading VOG virus-only flags from {args.vog_virusonly}", file=sys.stderr)
        vog_virusonly = load_vog_virusonly(args.vog_virusonly)
        print(f"[organize_genes] Loaded virus-only flags for {len(vog_virusonly)} VOGs", file=sys.stderr)

    if args.vog_categories and Path(args.vog_categories).exists():
        print(f"[organize_genes] Loading functional categories from {args.vog_categories}", file=sys.stderr)
        vog_categories = load_functional_categories(args.vog_categories)
        print(f"[organize_genes] Loaded {len(vog_categories)} functional category definitions", file=sys.stderr)
    
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
    
    def _build_vog_fields(gene_id):
        """Build VOG annotation fields for a gene (returns dict)."""
        hit = vog_hits.get(gene_id)
        if hit is None:
            return {
                'vog_id': '',
                'vog_evalue': '',
                'vog_score': '',
                'functional_category': '',
                'functional_category_name': '',
                'functional_description': '',
                'virus_only_high': '',
                'virus_only_medium': '',
                'virus_only_low': '',
            }
        vid = hit['vog_id']
        ann = vog_annotations.get(vid, {})
        vo = vog_virusonly.get(vid, {})
        func_cat = ann.get('functional_category', '')
        return {
            'vog_id': vid,
            'vog_evalue': f"{hit['evalue']:.2e}",
            'vog_score': f"{hit['score']:.1f}",
            'functional_category': func_cat,
            'functional_category_name': resolve_category_name(func_cat, vog_categories),
            'functional_description': ann.get('description', ''),
            'virus_only_high': '1' if vo.get('high') else ('0' if vo else ''),
            'virus_only_medium': '1' if vo.get('medium') else ('0' if vo else ''),
            'virus_only_low': '1' if vo.get('low') else ('0' if vo else ''),
        }

    def _process_gene(gene_id, contig_id, start, end, strand, header, aa_seq, nt_seq, tax_info):
        """Shared logic for processing a single gene record."""
        class_dir, file_basename = get_taxonomy_key(tax_info['parsed_lineage'])

        genes_by_taxon[(class_dir, file_basename)].append({
            'gene_id': gene_id,
            'header': header,
            'aa_seq': aa_seq,
            'nt_seq': nt_seq,
            'contig_id': contig_id,
            'start': start,
            'end': end,
            'strand': strand,
            'lineage': tax_info['lineage'],
        })

        row = {
            'gene_id': gene_id,
            'contig_id': contig_id,
            'start': start,
            'end': end,
            'strand': strand,
            'class_dir': class_dir,
            'file_basename': file_basename,
            'lineage': tax_info['lineage'],
        }
        if has_vog:
            row.update(_build_vog_fields(gene_id))
        mapping_data.append(row)

    # Read proteins file
    if HAS_BIOPYTHON:
        for record in SeqIO.parse(proteins_file, "fasta"):
            header = record.description
            aa_seq = str(record.seq)
            gene_id, contig_id, start, end, strand = parse_prodigal_header(header)
            tax_info = contig_taxonomy.get(contig_id, {
                'taxon_id': '0', 'lineage': '', 'parsed_lineage': []
            })
            nt_seq = ""
            if contig_id in contig_seqs and start and end:
                nt_seq = extract_gene_nucleotide(contig_seqs[contig_id], start, end, strand)
            _process_gene(gene_id, contig_id, start, end, strand, header, aa_seq, nt_seq, tax_info)
    else:
        # Simple parsing without BioPython
        current_header = None
        current_seq = []

        def _flush_record():
            if current_header is None:
                return
            aa_seq = ''.join(current_seq)
            gene_id, contig_id, start, end, strand = parse_prodigal_header(current_header)
            tax_info = contig_taxonomy.get(contig_id, {
                'taxon_id': '0', 'lineage': '', 'parsed_lineage': []
            })
            nt_seq = ""
            if contig_id in contig_seqs and start and end:
                nt_seq = extract_gene_nucleotide(contig_seqs[contig_id], start, end, strand)
            _process_gene(gene_id, contig_id, start, end, strand, current_header, aa_seq, nt_seq, tax_info)

        with open(proteins_file) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    _flush_record()
                    current_header = line[1:]
                    current_seq = []
                else:
                    current_seq.append(line)
            _flush_record()
    
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
    
    # Write gene annotation TSV (enriched with VOG data when available)
    base_columns = ["gene_id", "contig_id", "start", "end", "strand", "class_dir", "file_basename", "lineage"]
    vog_columns = ["vog_id", "vog_evalue", "vog_score", "functional_category",
                    "functional_category_name", "functional_description",
                    "virus_only_high", "virus_only_medium", "virus_only_low"]
    columns = base_columns + (vog_columns if has_vog else [])

    annotation_file = out_dir / "gene_annotation.tsv"
    with open(annotation_file, 'w') as f:
        f.write("\t".join(columns) + "\n")
        for row in mapping_data:
            strand_str = "+" if row['strand'] == 1 else "-"
            vals = [
                row['gene_id'], row['contig_id'],
                str(row['start'] or ''), str(row['end'] or ''),
                strand_str, row['class_dir'], row['file_basename'], row['lineage'],
            ]
            if has_vog:
                vals.extend([
                    row.get('vog_id', ''),
                    row.get('vog_evalue', ''),
                    row.get('vog_score', ''),
                    row.get('functional_category', ''),
                    row.get('functional_category_name', ''),
                    row.get('functional_description', ''),
                    row.get('virus_only_high', ''),
                    row.get('virus_only_medium', ''),
                    row.get('virus_only_low', ''),
                ])
            f.write("\t".join(vals) + "\n")

    # Backwards-compatible symlink
    compat_link = out_dir / "gene_taxonomy_mapping.tsv"
    if compat_link.exists() or compat_link.is_symlink():
        compat_link.unlink()
    try:
        compat_link.symlink_to("gene_annotation.tsv")
    except OSError:
        # Fallback: hard copy on filesystems that don't support symlinks
        import shutil
        shutil.copy2(annotation_file, compat_link)

    # Write functional summary (only when VOG data is available)
    if has_vog:
        func_counts = defaultdict(lambda: {'total': 0, 'virus_only': 0})
        genes_with_hit = 0
        genes_without_hit = 0

        for row in mapping_data:
            vid = row.get('vog_id', '')
            if vid:
                genes_with_hit += 1
                cat = row.get('functional_category', '') or 'Unknown'
                cat_name = row.get('functional_category_name', '') or cat
                key = f"{cat} ({cat_name})" if cat_name != cat else cat
                func_counts[key]['total'] += 1
                if row.get('virus_only_medium') == '1':
                    func_counts[key]['virus_only'] += 1
            else:
                genes_without_hit += 1

        summary_file = out_dir / "functional_summary.tsv"
        with open(summary_file, 'w') as f:
            f.write("functional_category\tgene_count\tvirus_only_count\n")
            for cat, counts in sorted(func_counts.items(), key=lambda x: -x[1]['total']):
                f.write(f"{cat}\t{counts['total']}\t{counts['virus_only']}\n")
            f.write(f"No VOG hit\t{genes_without_hit}\t\n")

        print(f"[organize_genes] VOG hits: {genes_with_hit} genes annotated, "
              f"{genes_without_hit} without hit", file=sys.stderr)
        print(f"[organize_genes] Functional summary written to {summary_file}", file=sys.stderr)

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
