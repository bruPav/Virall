"""
Gene prediction and annotation module for viral sequences.
"""

import os
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
from loguru import logger

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from .vog_annotator import VOGAnnotator


class ViralGenePredictor:
    """
    Comprehensive gene prediction and annotation for viral sequences.
    """
    
    def __init__(self, threads: int = 8, vog_db_path: Optional[str] = None, config: Optional[Dict] = None):
        """Initialize the gene predictor."""
        self.threads = threads
        self.config = config or {}
        self.temp_dir = Path(tempfile.mkdtemp(prefix="viral_genes_"))
        self.vog_annotator = VOGAnnotator(vog_db_path=vog_db_path, config=config, threads=threads)
        logger.info("ViralGenePredictor initialized")
    
    def predict_genes_comprehensive(
        self,
        contigs_file: str,
        viral_classifications: Dict[str, Dict],
        output_dir: Optional[str] = None
    ) -> Dict[str, Dict]:
        """
        Comprehensive gene prediction based on viral type.
        
        Args:
            contigs_file: Path to viral contigs file
            viral_classifications: DIAMOND classification results
            output_dir: Output directory for results
            
        Returns:
            Dictionary with gene prediction results
        """
        if output_dir is None:
            output_dir = self.temp_dir / "gene_predictions"
        
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)
        
        # Store contigs file path for nucleotide extraction
        self.contigs_file = contigs_file
        logger.info(f"Stored contigs file path: {self.contigs_file}")
        
        logger.info(f"Running comprehensive gene prediction on {contigs_file}")
        
        # Check if we're in RNA mode or Kaiju mode
        # Look for Kaiju classification results in the viral_classifications
        is_rna_mode = False
        is_kaiju_mode = False
        
        # Check for Kaiju classification results
        for classification in viral_classifications.values():
            if isinstance(classification, dict):
                # Check if this is a Kaiju classification result
                if classification.get('method') == 'kaiju':
                    is_kaiju_mode = True
                    # Check if it's RNA mode by looking at the classification method
                    if 'rna' in str(classification.get('classification_method', '')).lower():
                        is_rna_mode = True
                # Also check for the old classification_method field
                elif classification.get('classification_method') == 'kaiju_rna_mode':
                    is_rna_mode = True
                    is_kaiju_mode = True
                elif classification.get('classification_method') in ['kaiju_dna_mode', 'kaiju_rna_mode']:
                    is_kaiju_mode = True
        
        # If we have any Kaiju results, we're in Kaiju mode
        if not is_kaiju_mode and viral_classifications:
            # Check if any classification has Kaiju-related content
            for classification in viral_classifications.values():
                if isinstance(classification, dict) and 'kaiju' in str(classification).lower():
                    is_kaiju_mode = True
                    break
        
        # Note: RNA mode detection removed - it was incorrectly detecting short viral contigs as RNA
        # RNA mode should only be set explicitly in the configuration, not inferred from sequence length
        
        if is_rna_mode or is_kaiju_mode:
            # Both RNA mode and Kaiju mode use Prodigal for gene prediction
            mode_name = "RNA" if is_rna_mode else "Kaiju"
            logger.info(f"{mode_name} mode detected - running Prodigal gene prediction directly on contigs")
            kaiju_genes = self._run_prodigal_gene_prediction(contigs_file, output_dir)
        else:
            # Legacy mode (should not be used anymore)
            logger.info("Legacy mode detected - this should not happen with Kaiju approach")
            kaiju_genes = {}
        
        # Run additional gene prediction methods
        additional_predictions = self._run_additional_gene_prediction(
            contigs_file, viral_classifications, output_dir
        )
        
        # Combine and annotate all gene predictions
        combined_results = self._combine_gene_predictions(
            kaiju_genes, additional_predictions, output_dir, viral_classifications
        )
        
        # Annotate proteins
        protein_annotations = self._annotate_proteins(combined_results, output_dir)
        
        # VOG annotation for functional classification
        vog_annotations = self._run_vog_annotation(combined_results, output_dir, contigs_file)

        # Additionally, write per-species protein FASTA files using Kaiju classifications
        # Build a simple contig_id -> { 'classification': species_name } map from provided Kaiju results
        try:
            kaiju_species_map = {}
            for key, value in viral_classifications.items():
                if isinstance(value, dict) and value.get('method') == 'kaiju' and 'classifications' in value:
                    for contig_id, cls in value['classifications'].items():
                        # Prefer explicit classification/taxon_name fields
                        species_label = cls.get('classification') or cls.get('taxon_name') or 'Unknown'
                        kaiju_species_map[contig_id] = {'classification': species_label}

            if kaiju_species_map:
                self._write_species_protein_files(combined_results, kaiju_species_map, output_dir, contigs_file)
            else:
                logger.info("No Kaiju species map found in provided classifications - keeping VOG-based grouping only")
        except Exception as e:
            logger.warning(f"Failed to write species protein files from Kaiju classifications: {e}")
        
        # Generate summary
        summary = self._generate_gene_summary(combined_results, protein_annotations, vog_annotations)
        
        results = {
            "status": "completed",
            "kaiju_genes": kaiju_genes,
            "additional_predictions": additional_predictions,
            "combined_results": combined_results,
            "protein_annotations": protein_annotations,
            "vog_annotations": vog_annotations,
            "summary": summary,
            "output_directory": str(output_dir)
        }
        
        logger.info("Comprehensive gene prediction completed")
        return results
    
    def _extract_kaiju_genes(self, viral_classifications: Dict[str, Dict]) -> Dict[str, List[Dict]]:
        """Extract gene predictions from classification results."""
        logger.info("Extracting gene predictions from classification results")
        
        kaiju_genes = {}
        
        for contig_id, classification in viral_classifications.items():
            # Look for classification output files
            classification_output_dir = classification.get('output_dir')
            if not classification_output_dir:
                continue
            
            output_path = Path(classification_output_dir)
            
            # Find gene prediction files
            gff_file = output_path / "iter-0" / "all.pdg.gff"
            faa_file = output_path / "iter-0" / "all.pdg.faa"
            
            if gff_file.exists() and faa_file.exists():
                genes = self._parse_gff_file(gff_file, contig_id)
                proteins = self._parse_faa_file(faa_file, contig_id)
                
                # Combine gene and protein information
                for gene in genes:
                    gene_id = gene.get('id', gene.get('gene_id', 'unknown'))
                    # Try multiple matching strategies
                    protein_seq = None
                    if gene_id in proteins:
                        protein_seq = proteins[gene_id]
                    else:
                        # Try with full contig ID prefix
                        full_gene_id = f"{contig_id}||{gene_id}"
                        if full_gene_id in proteins:
                            protein_seq = proteins[full_gene_id]
                        else:
                            # Try partial matching
                            for protein_id, seq in proteins.items():
                                if gene_id in protein_id or protein_id.endswith(gene_id):
                                    protein_seq = seq
                                    break
                    
                    if protein_seq:
                        gene['protein_sequence'] = protein_seq
                
                kaiju_genes[contig_id] = genes
                logger.info(f"Extracted {len(genes)} genes for {contig_id}")
        
        return kaiju_genes
    
    def _parse_gff_file(self, gff_file: Path, contig_id: str) -> List[Dict]:
        """Parse GFF file to extract gene information."""
        genes = []
        
        try:
            with open(gff_file, 'r') as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    
                    fields = line.strip().split('\t')
                    if len(fields) < 9:
                        continue
                    
                    # Only process genes from the target contig
                    if contig_id not in fields[0]:
                        continue
                    
                    if fields[2] == 'CDS':  # Coding sequence
                        gene_info = {
                            'contig_id': fields[0],
                            'source': fields[1],
                            'type': fields[2],
                            'start': int(fields[3]),
                            'end': int(fields[4]),
                            'score': float(fields[5]) if fields[5] != '.' else 0.0,
                            'strand': fields[6],
                            'phase': int(fields[7]) if fields[7] != '.' else 0,
                            'attributes': fields[8]
                        }
                        
                        # Parse attributes
                        attrs = {}
                        for attr in fields[8].split(';'):
                            if '=' in attr:
                                key, value = attr.split('=', 1)
                                attrs[key] = value
                        
                        gene_info['id'] = attrs.get('ID', '')
                        gene_info['partial'] = attrs.get('partial', '00')
                        gene_info['start_type'] = attrs.get('start_type', '')
                        gene_info['rbs_motif'] = attrs.get('rbs_motif', '')
                        gene_info['gc_content'] = float(attrs.get('gc_cont', 0))
                        gene_info['confidence'] = float(attrs.get('conf', 0))
                        
                        genes.append(gene_info)
        
        except Exception as e:
            logger.warning(f"Failed to parse GFF file {gff_file}: {e}")
        
        return genes
    
    def _parse_faa_file(self, faa_file: Path, contig_id: str) -> Dict[str, str]:
        """Parse FASTA file to extract protein sequences."""
        proteins = {}
        
        try:
            for record in SeqIO.parse(faa_file, "fasta"):
                if contig_id in record.id:
                    # Extract gene ID from ID= in description
                    gene_id = None
                    if 'ID=' in record.description:
                        import re
                        match = re.search(r'ID=([^;]+)', record.description)
                        if match:
                            gene_id = match.group(1)
                    
                    # Fallback patterns if ID= not found
                    if not gene_id:
                        if '||' in record.id:
                            # Format: contig_id||gene_id
                            gene_id = record.id.split('||')[1]
                        elif '_' in record.id:
                            # Alternative format: contig_id_gene_id
                            parts = record.id.split('_')
                            if len(parts) > 1:
                                gene_id = parts[-1]
                    
                    # Use extracted gene_id or fallback to full record ID
                    if gene_id:
                        proteins[gene_id] = str(record.seq)
                    else:
                        proteins[record.id] = str(record.seq)
                    
                    # Also store with full ID for flexibility
                    proteins[record.id] = str(record.seq)
        
        except Exception as e:
            logger.warning(f"Failed to parse FASTA file {faa_file}: {e}")
        
        return proteins

    def _run_prodigal_gene_prediction(self, contigs_file: str, output_dir: Path) -> Dict[str, List[Dict]]:
        """Run Prodigal gene prediction directly on contigs."""
        logger.info("Running Prodigal gene prediction on viral contigs")
        
        # Run Prodigal on the contigs file
        prodigal_output = output_dir / "prodigal"
        prodigal_output.mkdir(exist_ok=True)
        
        gff_file = prodigal_output / "prodigal.gff"
        faa_file = prodigal_output / "prodigal.faa"
        
        # Run Prodigal with viral-specific parameters
        cmd = [
            "prodigal",
            "-i", contigs_file,
            "-o", str(gff_file),
            "-a", str(faa_file),
            "-f", "gff",
            "-p", "meta"  # Use meta mode for better gene prediction
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                logger.error(f"Prodigal failed: {result.stderr}")
                return {}
            
            # Parse the GFF file to extract gene information
            genes = self._parse_prodigal_gff(str(gff_file))
            
            # Group genes by contig
            contig_genes = {}
            for gene in genes:
                contig_id = gene['contig_id']
                if contig_id not in contig_genes:
                    contig_genes[contig_id] = []
                contig_genes[contig_id].append(gene)
            
            logger.info(f"Prodigal predicted {len(genes)} genes across {len(contig_genes)} contigs")
            return contig_genes
            
        except Exception as e:
            logger.error(f"Error running Prodigal: {e}")
            return {}

    def _parse_prodigal_gff(self, gff_file: str) -> List[Dict]:
        """Parse Prodigal GFF file to extract gene information."""
        genes = []
        
        try:
            with open(gff_file, 'r') as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    
                    parts = line.strip().split('\t')
                    if len(parts) >= 9 and parts[2] == 'CDS':
                        contig_id = parts[0]
                        start = int(parts[3])
                        end = int(parts[4])
                        strand = parts[6]
                        
                        # Parse attributes
                        attributes = {}
                        for attr in parts[8].split(';'):
                            if '=' in attr:
                                key, value = attr.split('=', 1)
                                attributes[key.strip()] = value.strip()
                        
                        gene_id = attributes.get('ID', f"{contig_id}_{start}_{end}")
                        
                        genes.append({
                            'gene_id': gene_id,
                            'contig_id': contig_id,
                            'start': start,
                            'end': end,
                            'strand': strand,
                            'attributes': attributes
                        })
        
        except Exception as e:
            logger.error(f"Error parsing Prodigal GFF file: {e}")
        
        return genes
    
    def _run_additional_gene_prediction(
        self,
        contigs_file: str,
        viral_classifications: Dict[str, Dict],
        output_dir: Path
    ) -> Dict[str, Dict]:
        """Run additional gene prediction methods based on viral type."""
        logger.info("Running additional gene prediction methods")
        
        additional_predictions = {}
        
        # Group contigs by viral type
        viral_groups = self._group_contigs_by_type(contigs_file, viral_classifications)
        
        for viral_type, contig_ids in viral_groups.items():
            logger.info(f"Predicting genes for {viral_type}: {len(contig_ids)} contigs")
            
            # Sanitize viral type name for file system compatibility
            safe_viral_type = self._sanitize_filename(viral_type)
            
            # Create temporary file with contigs of this type
            type_file = output_dir / f"{safe_viral_type}_contigs.fasta"
            self._extract_contigs_by_ids(contigs_file, contig_ids, type_file)
            
            # Run appropriate gene predictor
            if viral_type in ['dsDNAphage', 'NCLDV']:
                # Use Prodigal for DNA viruses
                genes = self._run_prodigal(type_file, output_dir / f"{safe_viral_type}_prodigal")
            elif viral_type in ['RNA']:
                # Use GeneMarkS for RNA viruses
                genes = self._run_genemark(type_file, output_dir / f"{safe_viral_type}_genemark")
            else:
                # Use Prodigal as default
                genes = self._run_prodigal(type_file, output_dir / f"{safe_viral_type}_prodigal")
            
            additional_predictions[viral_type] = genes
        
        return additional_predictions
    
    def _sanitize_filename(self, filename: str) -> str:
        """Sanitize filename by replacing invalid characters."""
        import re
        # Replace spaces and other problematic characters with underscores
        sanitized = re.sub(r'[/\\:*?"<>|\s]+', '_', filename)
        # Replace multiple underscores with single underscore
        sanitized = re.sub(r'_+', '_', sanitized)
        # Remove leading/trailing underscores
        sanitized = sanitized.strip('_')
        return sanitized
    
    def _group_contigs_by_type(
        self,
        contigs_file: str,
        viral_classifications: Dict[str, Dict]
    ) -> Dict[str, List[str]]:
        """Group contigs by their viral type."""
        viral_groups = {}
        
        # Extract actual contig classifications from Kaiju or DIAMOND results
        for contig_id, contig_classification in viral_classifications.items():
            # Handle both Kaiju and DIAMOND formats
            if 'classification' in contig_classification:
                # Kaiju format
                viral_type = contig_classification.get('classification', 'unknown')
                if viral_type == 'Unclassified':
                    viral_type = 'unknown'
            else:
                # DIAMOND format - use consensus species as viral type
                consensus = contig_classification.get('consensus', {})
                viral_type = consensus.get('species', 'unknown')
                if viral_type == 'unknown':
                    viral_type = consensus.get('family', 'unknown')
            
            if viral_type not in viral_groups:
                viral_groups[viral_type] = []
            viral_groups[viral_type].append(contig_id)
        
        return viral_groups
    
    def _extract_contigs_by_ids(
        self,
        contigs_file: str,
        contig_ids: List[str],
        output_file: Path
    ) -> None:
        """Extract specific contigs by their IDs."""
        contig_ids_set = set(contig_ids)
        extracted_contigs = []
        
        for record in SeqIO.parse(contigs_file, "fasta"):
            # Try multiple matching strategies
            contig_id = record.id.split()[0]  # Get main ID part
            full_id = record.id  # Full ID
            
            # Check if any part of the contig ID matches
            found = False
            for target_id in contig_ids_set:
                if (contig_id == target_id or 
                    full_id == target_id or 
                    contig_id in target_id or 
                    target_id in contig_id):
                    found = True
                    break
            
            if found:
                extracted_contigs.append(record)
        
        SeqIO.write(extracted_contigs, output_file, "fasta")
        logger.info(f"Extracted {len(extracted_contigs)} contigs to {output_file}")
        
        if len(extracted_contigs) == 0:
            logger.warning(f"No contigs found matching IDs: {contig_ids}")
            logger.warning(f"Available contig IDs in file:")
            for record in SeqIO.parse(contigs_file, "fasta"):
                logger.warning(f"  {record.id}")
                break  # Just show first one as example
    
    def _run_prodigal(self, contigs_file: Path, output_dir: Path) -> List[Dict]:
        """Run Prodigal gene prediction."""
        output_dir.mkdir(exist_ok=True)
        
        gff_file = output_dir / "prodigal.gff"
        faa_file = output_dir / "prodigal.faa"
        
        cmd = [
            "prodigal",
            "-i", str(contigs_file),
            "-o", str(gff_file),
            "-a", str(faa_file),
            "-f", "gff",
            "-p", "meta"  # Metagenomic mode
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode == 0:
                genes = self._parse_gff_file(gff_file, "")
                logger.info(f"Prodigal predicted {len(genes)} genes")
                return genes
            else:
                logger.warning(f"Prodigal failed: {result.stderr}")
                return []
        except FileNotFoundError:
            logger.warning("Prodigal not found, skipping gene prediction")
            return []
    
    def _run_genemark(self, contigs_file: Path, output_dir: Path) -> List[Dict]:
        """Run GeneMark gene prediction."""
        output_dir.mkdir(exist_ok=True)
        
        # GeneMark requires specific setup, so we'll use Prodigal as fallback
        logger.info("GeneMark not implemented, using Prodigal as fallback")
        return self._run_prodigal(contigs_file, output_dir)
    
    def _combine_gene_predictions(
        self,
        kaiju_genes: Dict[str, List[Dict]],
        additional_predictions: Dict[str, Dict],
        output_dir: Path,
        viral_classifications: Dict[str, Dict]
    ) -> Dict[str, List[Dict]]:
        """Combine gene predictions from different methods."""
        logger.info("Combining gene predictions")
        
        combined_results = {}
        
        # Start with Kaiju predictions
        logger.info(f"Processing {len(kaiju_genes)} contigs from kaiju_genes")
        for contig_id, genes in kaiju_genes.items():
            logger.info(f"Processing contig {contig_id} with {len(genes)} genes")
            # Get viral type from classification data
            viral_type = 'unknown'
            if contig_id in viral_classifications:
                # Try different possible keys for viral type
                viral_type = (viral_classifications[contig_id].get('viral_group') or 
                             viral_classifications[contig_id].get('taxon_name') or 
                             viral_classifications[contig_id].get('classification') or 
                             'unknown')
            
            logger.info(f"Calling _add_protein_sequences for contig {contig_id} with viral_type {viral_type}")
            # Extract protein sequences for Kaiju genes
            genes_with_proteins = self._add_protein_sequences(genes, contig_id, viral_type, output_dir)
            
            combined_results[contig_id] = {
                'kaiju_genes': genes_with_proteins,
                'total_genes': len(genes_with_proteins),
                'method': 'kaiju',
                'viral_type': viral_type
            }
        
        # Process additional predictions (Prodigal results in RNA mode)
        logger.info(f"Processing {len(additional_predictions)} viral types from additional_predictions")
        for viral_type, genes in additional_predictions.items():
            logger.info(f"Processing viral type {viral_type} with {len(genes)} genes")
            
            # Group genes by contig for processing
            contig_genes = {}
            for gene in genes:
                contig_id = gene.get('contig_id', '').split('||')[0]
                if contig_id not in contig_genes:
                    contig_genes[contig_id] = []
                contig_genes[contig_id].append(gene)
            
            # Process each contig's genes
            for contig_id, contig_genes_list in contig_genes.items():
                logger.info(f"Processing additional prediction contig {contig_id} with {len(contig_genes_list)} genes")
                logger.info(f"Calling _add_protein_sequences for contig {contig_id} with viral_type {viral_type}")
                
                # Extract protein sequences for additional prediction genes
                genes_with_proteins = self._add_protein_sequences(contig_genes_list, contig_id, viral_type, output_dir)
                
                if contig_id in combined_results:
                    # Add to existing results
                    combined_results[contig_id]['additional_genes'] = genes_with_proteins
                    combined_results[contig_id]['total_genes'] += len(genes_with_proteins)
                else:
                    # Create new entry
                    combined_results[contig_id] = {
                        'additional_genes': genes_with_proteins,
                        'total_genes': len(genes_with_proteins),
                        'method': 'additional',
                        'viral_type': viral_type
                    }
        
        
        # Write combined results
        self._write_combined_results(combined_results, output_dir)
        
        return combined_results
    
    def _add_protein_sequences(self, genes: List[Dict], contig_id: str, viral_type: str, output_dir: Path) -> List[Dict]:
        """Add protein sequences to gene predictions by reading from FAA files."""
        genes_with_proteins = []
        
        # Look for FAA files in the output directory
        # Check for generic prodigal directory first (new format)
        faa_files = list(output_dir.glob("prodigal/prodigal.faa"))
        logger.info(f"Looking for generic prodigal files: {faa_files}")
        
        # If not found, try the old RNA mode directory
        if not faa_files:
            faa_files = list(output_dir.glob("prodigal_rna/prodigal.faa"))
            logger.info(f"Looking for RNA mode prodigal files: {faa_files}")
        
        # If not found, try the normal viral type pattern
        if not faa_files:
            # Sanitize the viral type name to match the directory name
            safe_viral_type = self._sanitize_filename(viral_type)
            faa_files = list(output_dir.glob(f"{safe_viral_type}_prodigal/prodigal.faa"))
            logger.info(f"Looking for viral type prodigal files: {safe_viral_type}_prodigal/prodigal.faa -> {faa_files}")
        
        # If still not found, try to find any prodigal.faa file in the output directory
        if not faa_files:
            faa_files = list(output_dir.glob("**/prodigal.faa"))
            logger.info(f"Looking for any prodigal.faa files: {faa_files}")
        
        if not faa_files:
            logger.warning(f"No FAA file found for {viral_type} contigs")
            return genes
        
        faa_file = faa_files[0]
        protein_sequences = self._parse_faa_file(faa_file, contig_id)
        
        for gene in genes:
            gene_id = gene.get('id', gene.get('gene_id', ''))
            if gene_id in protein_sequences:
                gene['protein_sequence'] = protein_sequences[gene_id]
            else:
                # Try alternative ID matching
                for protein_id, sequence in protein_sequences.items():
                    if gene_id in protein_id or protein_id in gene_id:
                        gene['protein_sequence'] = sequence
                        break
                else:
                    logger.debug(f"No protein sequence found for gene {gene_id}")
                    gene['protein_sequence'] = ''
            
            genes_with_proteins.append(gene)
        
        # Log summary only for significant numbers to avoid spam
        if len(genes_with_proteins) > 50:
            logger.info(f"Added protein sequences for {len(genes_with_proteins)} genes from {viral_type}")
        return genes_with_proteins
    
    def _write_combined_results(
        self,
        combined_results: Dict[str, List[Dict]],
        output_dir: Path
    ) -> None:
        """Write combined gene prediction results to files."""
        
        # Write summary table
        summary_file = output_dir / "gene_prediction_summary.tsv"
        with open(summary_file, 'w') as f:
            f.write("contig_id\tviral_type\ttotal_genes\tmethod\tavg_gene_length\tavg_gc_content\n")
            
            for contig_id, data in combined_results.items():
                # Get genes from both kaiju_genes and additional_genes
                kaiju_genes = data.get('kaiju_genes', [])
                additional_genes = data.get('additional_genes', [])
                all_genes = kaiju_genes + additional_genes
                
                if all_genes:
                    avg_length = sum(g['end'] - g['start'] + 1 for g in all_genes) / len(all_genes)
                    avg_gc = sum(g.get('gc_content', 0) for g in all_genes) / len(all_genes)
                    f.write(f"{contig_id}\t{data.get('viral_type', 'unknown')}\t"
                           f"{len(all_genes)}\t{data.get('method', 'unknown')}\t"
                           f"{avg_length:.1f}\t{avg_gc:.3f}\n")
        
        # Write detailed GFF
        gff_file = output_dir / "combined_genes.gff"
        with open(gff_file, 'w') as f:
            f.write("##gff-version 3\n")
            f.write("# Combined gene predictions from Kaiju and additional methods\n")
            
            for contig_id, data in combined_results.items():
                # Get genes from both kaiju_genes and additional_genes
                kaiju_genes = data.get('kaiju_genes', [])
                additional_genes = data.get('additional_genes', [])
                all_genes = kaiju_genes + additional_genes
                
                for gene in all_genes:
                    # Handle different gene formats (Kaiju vs Prodigal)
                    source = gene.get('source', 'Prodigal')
                    gene_type = gene.get('type', 'CDS')
                    score = gene.get('score', '.')
                    phase = gene.get('phase', '.')
                    attributes = gene.get('attributes', {})
                    
                    # Convert attributes dict to string if needed
                    if isinstance(attributes, dict):
                        attr_str = ';'.join([f"{k}={v}" for k, v in attributes.items()])
                    else:
                        attr_str = str(attributes)
                    
                    f.write(f"{gene['contig_id']}\t{source}\t{gene_type}\t"
                           f"{gene['start']}\t{gene['end']}\t{score}\t"
                           f"{gene['strand']}\t{phase}\t{attr_str}\n")
        
        logger.info(f"Combined results written to {output_dir}")
    
    def _annotate_proteins(
        self,
        combined_results: Dict[str, List[Dict]],
        output_dir: Path
    ) -> Dict[str, Dict]:
        """Annotate proteins using multiple databases."""
        logger.info("Annotating proteins")
        
        # This would integrate with BLAST, HMMER, DIAMOND, etc.
        # For now, we'll create a placeholder structure
        
        protein_annotations = {}
        
        for contig_id, data in combined_results.items():
            # Get genes from both kaiju_genes and additional_genes
            kaiju_genes = data.get('kaiju_genes', [])
            additional_genes = data.get('additional_genes', [])
            all_genes = kaiju_genes + additional_genes
            annotations = []
            
            for gene in all_genes:
                if 'protein_sequence' in gene:
                    # Placeholder for protein annotation
                    annotation = {
                        'gene_id': gene.get('id', gene.get('gene_id', 'unknown')),
                        'length': len(gene['protein_sequence']),
                        'molecular_weight': self._calculate_molecular_weight(gene['protein_sequence']),
                        'isoelectric_point': self._calculate_pi(gene['protein_sequence']),
                        'functional_prediction': 'Unknown',
                        'confidence': gene.get('confidence', 0)
                    }
                    annotations.append(annotation)
            
            protein_annotations[contig_id] = annotations
        
        # Write protein annotations
        annotation_file = output_dir / "protein_annotations.tsv"
        with open(annotation_file, 'w') as f:
            f.write("contig_id\tgene_id\tlength\tmolecular_weight\tisoelectric_point\tfunctional_prediction\tconfidence\n")
            
            for contig_id, annotations in protein_annotations.items():
                for ann in annotations:
                    f.write(f"{contig_id}\t{ann['gene_id']}\t{ann['length']}\t"
                           f"{ann['molecular_weight']:.1f}\t{ann['isoelectric_point']:.2f}\t"
                           f"{ann['functional_prediction']}\t{ann['confidence']:.2f}\n")
        
        return protein_annotations
    
    def _calculate_molecular_weight(self, protein_sequence: str) -> float:
        """Calculate molecular weight of protein."""
        # Simple calculation - in practice, use BioPython
        aa_weights = {
            'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.2,
            'Q': 146.1, 'E': 147.1, 'G': 75.1, 'H': 155.2, 'I': 131.2,
            'L': 131.2, 'K': 146.2, 'M': 149.2, 'F': 165.2, 'P': 115.1,
            'S': 105.1, 'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1
        }
        
        weight = sum(aa_weights.get(aa, 0) for aa in protein_sequence.upper())
        return weight
    
    def _calculate_pi(self, protein_sequence: str) -> float:
        """Calculate isoelectric point of protein."""
        # Simplified calculation - in practice, use proper algorithms
        return 7.0  # Placeholder
    
    def _run_vog_annotation(self, combined_results: Dict[str, List[Dict]], output_dir: Path, contigs_file: str = None) -> Dict[str, Dict]:
        """Run VOG annotation on predicted proteins."""
        if not self.vog_annotator.vog_db_path or not self.vog_annotator.vog_db_path.exists():
            logger.warning("VOG database not available, skipping VOG annotation")
            return {}
        
        logger.info("Running VOG annotation for functional classification")
        
        try:
            # Create protein FASTA file for VOG annotation
            protein_file = output_dir / "all_proteins.faa"
            self._write_protein_fasta(combined_results, protein_file)
            
            # Run VOG annotation
            vog_results = self.vog_annotator.annotate_proteins(str(protein_file), str(output_dir / "vog_annotation"))
            
            # Get viral classifications based on VOG
            viral_classifications = self.vog_annotator.get_viral_classification(vog_results)
            
            # Write species-specific protein files (only if classifications are available)
            if viral_classifications:
                # Use the contigs file path from the method parameter
                logger.info(f"Passing contigs file to species protein files: {contigs_file}")
                self._write_species_protein_files(combined_results, viral_classifications, output_dir, contigs_file)
            else:
                logger.info("No viral classifications available - skipping species-specific protein files")
            
            return {
                "vog_hits": vog_results,
                "viral_classifications": viral_classifications,
                "total_annotated_proteins": len(vog_results)
            }
            
        except Exception as e:
            logger.error(f"Error during VOG annotation: {e}")
            return {}
    
    def _write_protein_fasta(self, combined_results: Dict[str, List[Dict]], output_file: Path):
        """Write all predicted proteins to FASTA file."""
        protein_count = 0
        with open(output_file, 'w') as f:
            for contig_id, data in combined_results.items():
                # Handle Kaiju genes
                kaiju_genes = data.get('kaiju_genes', [])
                for gene in kaiju_genes:
                    if 'protein_sequence' in gene and gene['protein_sequence']:
                        gene_id = gene.get('id', f"{contig_id}_{protein_count}")
                        f.write(f">{gene_id}\n{gene['protein_sequence']}\n")
                        protein_count += 1
                
                # Handle additional genes
                additional_genes = data.get('additional_genes', [])
                for gene in additional_genes:
                    if 'protein_sequence' in gene and gene['protein_sequence']:
                        gene_id = gene.get('id', f"{contig_id}_{protein_count}")
                        f.write(f">{gene_id}\n{gene['protein_sequence']}\n")
                        protein_count += 1
        
        logger.info(f"Wrote {protein_count} protein sequences to {output_file}")

    def _write_species_protein_files(self, combined_results: Dict[str, List[Dict]], viral_classifications: Dict[str, Dict], output_dir: Path, contigs_file: str = None):
        """Write protein files grouped by species classification."""
        species_proteins = {}
        
        # Load VOG annotations for functional information
        vog_annotations = {}
        vog_annotation_file = output_dir / "vog_annotation" / "vog_annotations.tsv"
        if vog_annotation_file.exists():
            try:
                with open(vog_annotation_file, 'r') as f:
                    lines = f.readlines()
                    if len(lines) > 1:  # Skip header
                        for line in lines[1:]:
                            parts = line.strip().split('\t')
                            if len(parts) >= 7:
                                query_id = parts[0]
                                vog_id = parts[1]
                                function = parts[3]
                                category = parts[4]
                                consensus_function = parts[5]
                                
                                # Use consensus_function if available, otherwise function
                                annotation = consensus_function if consensus_function != 'Unknown' else function
                                vog_annotations[query_id] = {
                                    'vog_id': vog_id,
                                    'function': annotation,
                                    'category': category
                                }
                logger.info(f"Loaded {len(vog_annotations)} VOG annotations for protein naming")
            except Exception as e:
                logger.warning(f"Failed to load VOG annotations: {e}")
        
        # Load the protein sequences from all_proteins.faa to get the correct gene IDs
        protein_sequences = {}
        all_proteins_file = output_dir / "all_proteins.faa"
        if all_proteins_file.exists():
            try:
                with open(all_proteins_file, 'r') as f:
                    current_id = None
                    for line in f:
                        line = line.strip()
                        if line.startswith('>'):
                            current_id = line[1:]  # Remove '>' prefix
                        elif current_id and line:
                            protein_sequences[current_id] = line
                logger.info(f"Loaded {len(protein_sequences)} protein sequences from all_proteins.faa")
            except Exception as e:
                logger.warning(f"Failed to load protein sequences: {e}")
        
        # Group proteins by species
        logger.info(f"Viral classifications available: {list(viral_classifications.keys())}")
        for contig_id, data in combined_results.items():
            # Use the viral_type from the combined_results data instead of looking up in viral_classifications
            species_name = data.get('viral_type', 'Unknown')
            logger.info(f"Contig {contig_id} species_name from data: {species_name}")
            
            # Sanitize species name for filename
            safe_species_name = self._sanitize_filename(species_name)
            
            if safe_species_name not in species_proteins:
                species_proteins[safe_species_name] = []
            
            # Collect proteins from this contig
            kaiju_genes = data.get('kaiju_genes', [])
            logger.info(f"Contig {contig_id} has {len(kaiju_genes)} kaiju_genes")
            for gene in kaiju_genes:
                if 'protein_sequence' in gene and gene['protein_sequence']:
                    logger.info(f"Found protein sequence for gene {gene.get('id', 'unknown')}")
                    # Find the matching protein ID from all_proteins.faa by sequence
                    gene_sequence = gene['protein_sequence']
                    matching_protein_id = None
                    
                    for protein_id, sequence in protein_sequences.items():
                        if sequence == gene_sequence:
                            matching_protein_id = protein_id
                            break
                    
                    if matching_protein_id:
                        # Get VOG annotation using the matching protein ID
                        vog_info = vog_annotations.get(matching_protein_id, {})
                        function = vog_info.get('function', 'Unknown_function')
                        vog_id = vog_info.get('vog_id', '')
                        
                        # Create informative protein name
                        if vog_id and function != 'Unknown_function':
                            # Sanitize function name for FASTA header
                            safe_function = function.replace(' ', '_').replace('|', '_').replace(';', '_').replace(':', '_')
                            safe_function = safe_function[:50]  # Limit length
                            protein_name = f"{matching_protein_id}__{vog_id}__{safe_function}"
                        else:
                            protein_name = matching_protein_id
                        
                        species_proteins[safe_species_name].append({
                            'gene_id': matching_protein_id,
                            'protein_name': protein_name,
                            'contig_id': contig_id,
                            'sequence': gene['protein_sequence'],
                            'start': gene.get('start', 0),
                            'end': gene.get('end', 0),
                            'strand': gene.get('strand', '+'),
                            'function': function,
                            'vog_id': vog_id
                        })
            
            # Also collect from additional genes
            additional_genes = data.get('additional_genes', [])
            logger.info(f"Contig {contig_id} has {len(additional_genes)} additional_genes")
            for gene in additional_genes:
                if 'protein_sequence' in gene and gene['protein_sequence']:
                    logger.info(f"Found protein sequence for additional gene {gene.get('id', 'unknown')}")
                    # Find the matching protein ID from all_proteins.faa by sequence
                    gene_sequence = gene['protein_sequence']
                    matching_protein_id = None
                    
                    logger.info(f"Looking for matching protein for additional gene {gene.get('id', 'unknown')} with sequence length {len(gene_sequence)}")
                    for protein_id, sequence in protein_sequences.items():
                        if sequence == gene_sequence:
                            matching_protein_id = protein_id
                            logger.info(f"Found matching protein: {matching_protein_id}")
                            break
                    
                    if not matching_protein_id:
                        logger.warning(f"No matching protein found for additional gene {gene.get('id', 'unknown')} from contig {contig_id}")
                    
                    if matching_protein_id:
                        # Get VOG annotation using the matching protein ID
                        vog_info = vog_annotations.get(matching_protein_id, {})
                        function = vog_info.get('function', 'Unknown_function')
                        vog_id = vog_info.get('vog_id', '')
                        
                        # Create informative protein name
                        if vog_id and function != 'Unknown_function':
                            # Sanitize function name for FASTA header
                            safe_function = function.replace(' ', '_').replace('|', '_').replace(';', '_').replace(':', '_')
                            safe_function = safe_function[:50]  # Limit length
                            protein_name = f"{matching_protein_id}__{vog_id}__{safe_function}"
                        else:
                            protein_name = matching_protein_id
                        
                        species_proteins[safe_species_name].append({
                            'gene_id': matching_protein_id,
                            'protein_name': protein_name,
                            'contig_id': contig_id,
                            'sequence': gene['protein_sequence'],
                            'start': gene.get('start', 0),
                            'end': gene.get('end', 0),
                            'strand': gene.get('strand', '+'),
                            'function': function,
                            'vog_id': vog_id
                        })
        
        # Write protein and nucleotide files for each species
        species_dir = output_dir / "species_proteins"
        species_dir.mkdir(exist_ok=True)
        
        for species_name, proteins in species_proteins.items():
            if proteins:  # Only write if there are proteins
                # Write protein file
                species_protein_file = species_dir / f"{species_name}_proteins.faa"
                with open(species_protein_file, 'w') as f:
                    for protein in proteins:
                        # Write header with functional annotation
                        header = f">{protein['protein_name']} [contig={protein['contig_id']}] [pos={protein['start']}-{protein['end']}] [strand={protein['strand']}]"
                        if protein['function'] != 'Unknown_function':
                            header += f" [function={protein['function']}]"
                        f.write(f"{header}\n")
                        f.write(f"{protein['sequence']}\n")
                
                # Write nucleotide file
                species_nucleotide_file = species_dir / f"{species_name}_genes.fna"
                with open(species_nucleotide_file, 'w') as f:
                    for protein in proteins:
                        # Extract nucleotide sequence from contig
                        nucleotide_seq = self._extract_nucleotide_sequence(
                            protein['contig_id'], 
                            protein['start'], 
                            protein['end'], 
                            protein['strand'],
                            contigs_file
                        )
                        
                        if nucleotide_seq:
                            # Write header with functional annotation
                            header = f">{protein['protein_name']} [contig={protein['contig_id']}] [pos={protein['start']}-{protein['end']}] [strand={protein['strand']}]"
                            if protein['function'] != 'Unknown_function':
                                header += f" [function={protein['function']}]"
                            f.write(f"{header}\n")
                            f.write(f"{nucleotide_seq}\n")
                
                logger.info(f"Wrote {len(proteins)} proteins and genes for {species_name} to {species_dir}")
        
        # Write species summary
        summary_file = species_dir / "species_protein_summary.tsv"
        with open(summary_file, 'w') as f:
            f.write("species_name\tprotein_count\tcontig_count\n")
            for species_name, proteins in species_proteins.items():
                contigs = set(protein['contig_id'] for protein in proteins)
                f.write(f"{species_name}\t{len(proteins)}\t{len(contigs)}\n")
        
        logger.info(f"Species protein and nucleotide files written to {species_dir}")

    def _extract_nucleotide_sequence(self, contig_id: str, start: int, end: int, strand: str, contigs_file: str) -> str:
        """Extract nucleotide sequence from contig based on gene coordinates."""
        try:
            # Load contig sequences from the original contigs file
            if not contigs_file or not Path(contigs_file).exists():
                logger.warning(f"Contigs file not found: {contigs_file}")
                return ""
            
            from Bio import SeqIO
            for record in SeqIO.parse(contigs_file, "fasta"):
                if record.id == contig_id:
                    # Extract the sequence (convert to 0-based indexing)
                    start_idx = start - 1  # Convert to 0-based
                    end_idx = end
                    
                    # Ensure indices are within bounds
                    if start_idx < 0:
                        start_idx = 0
                    if end_idx > len(record.seq):
                        end_idx = len(record.seq)
                    
                    # Extract sequence
                    nucleotide_seq = str(record.seq[start_idx:end_idx])
                    
                    # Reverse complement if on negative strand
                    if strand == '-':
                        from Bio.Seq import Seq
                        nucleotide_seq = str(Seq(nucleotide_seq).reverse_complement())
                    
                    return nucleotide_seq
            
            logger.warning(f"Contig {contig_id} not found in contigs file")
            return ""
            
        except Exception as e:
            logger.warning(f"Error extracting nucleotide sequence for {contig_id}: {e}")
            return ""

    def _generate_gene_summary(
        self,
        combined_results: Dict[str, List[Dict]],
        protein_annotations: Dict[str, Dict],
        vog_annotations: Dict[str, Dict] = None
    ) -> Dict[str, Union[int, float, str]]:
        """Generate summary statistics for gene predictions."""
        
        total_contigs = len(combined_results)
        total_genes = sum(len(data.get('kaiju_genes', [])) for data in combined_results.values())
        
        if total_contigs > 0:
            avg_genes_per_contig = total_genes / total_contigs
        else:
            avg_genes_per_contig = 0
        
        # Calculate gene length statistics
        all_gene_lengths = []
        for data in combined_results.values():
            genes = data.get('kaiju_genes', [])
            for gene in genes:
                length = gene['end'] - gene['start'] + 1
                all_gene_lengths.append(length)
        
        if all_gene_lengths:
            avg_gene_length = sum(all_gene_lengths) / len(all_gene_lengths)
            min_gene_length = min(all_gene_lengths)
            max_gene_length = max(all_gene_lengths)
        else:
            avg_gene_length = min_gene_length = max_gene_length = 0
        
        # Add VOG annotation statistics
        vog_stats = {}
        if vog_annotations:
            vog_stats = {
                "vog_annotated_proteins": vog_annotations.get("total_annotated_proteins", 0),
                "vog_viral_classifications": len(vog_annotations.get("viral_classifications", {})),
                "vog_available": True
            }
        else:
            vog_stats = {
                "vog_annotated_proteins": 0,
                "vog_viral_classifications": 0,
                "vog_available": False
            }
        
        summary = {
            "total_contigs": total_contigs,
            "total_genes": total_genes,
            "avg_genes_per_contig": round(avg_genes_per_contig, 1),
            "avg_gene_length": round(avg_gene_length, 1),
            "min_gene_length": min_gene_length,
            "max_gene_length": max_gene_length,
            "status": "completed",
            **vog_stats
        }
        
        return summary
    
    def cleanup(self) -> None:
        """Clean up temporary files."""
        if self.temp_dir.exists():
            import shutil
            shutil.rmtree(self.temp_dir)
            logger.info("Cleaned up temporary files")
