"""
Command-line interface for Virall - comprehensive viral genome analysis.
"""

import os
import sys
from pathlib import Path
from typing import Optional, Dict

import click
import yaml
from loguru import logger

from .core.assembler import ViralAssembler
from .core.preprocessor import Preprocessor
from .core.viral_identifier import ViralIdentifier
from .core.validator import AssemblyValidator


@click.group()
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose logging')
@click.option('--log-file', help='Log file path')
def main(verbose: bool, log_file: Optional[str]):
    """Virall - Comprehensive viral genome analysis including assembly, classification, gene prediction, and annotation. Available commands: assemble, analyse, classify, preprocess, quantify, validate, annotate. Under development: train-model (machine learning model training)"""
    
    # Configure logging
    if verbose:
        logger.remove()
        logger.add(sys.stderr, level="DEBUG")
    
    if log_file:
        logger.add(log_file, level="DEBUG")
    
    logger.info("Virall started")


@main.command()
@click.option('--short-reads-1', help='Path to first mate of paired-end reads')
@click.option('--short-reads-2', help='Path to second mate of paired-end reads')
@click.option('--nanopore', help='Path to ONT Nanopore long reads')
@click.option('--pacbio', help='Path to PacBio long reads')
@click.option('--single-reads', help='Path to single-end reads')
@click.option('--reference', help='Optional reference genome for guided assembly')
@click.option('--output-dir', '-o', required=True, help='Output directory')
@click.option('--threads', '-t', default=8, help='Number of threads')
@click.option('--memory', default='16G', help='Memory allocation')
@click.option('--config', help='Configuration file path')
@click.option('--min-contig-length', default=1000, help='Minimum contig length')
@click.option('--viral-confidence', default=0.8, help='Viral sequence confidence threshold')
@click.option('--assembly-strategy', 
              type=click.Choice(['hybrid', 'short_only', 'long_only']),
              default='hybrid', help='Assembly strategy')
@click.option('--rna-mode', is_flag=True, help='Enable RNA-specific assembly parameters')
@click.option('--mem-efficient', '-m', is_flag=True, help='Enable memory-efficient mode with read subsampling for large datasets')
@click.option('--single-cell', is_flag=True, help='Enable single-cell sequencing mode')
@click.option('--cell-barcodes', help='Path to cell barcodes file (if not 10X format)')
@click.option('--min-cells', default=100, help='Minimum number of cells to process')
@click.option('--cellranger-path', help='Path to Cell Ranger installation')
@click.option('--barcode-whitelist', help='Path to barcode whitelist file')
@click.option('--index-reads', help='Path to index reads (I1) for single-cell data')
def assemble(
    short_reads_1: Optional[str],
    short_reads_2: Optional[str],
    nanopore: Optional[str],
    pacbio: Optional[str],
    single_reads: Optional[str],
    reference: Optional[str],
    output_dir: str,
    threads: int,
    memory: str,
    config: Optional[str],
    min_contig_length: int,
    viral_confidence: float,
    assembly_strategy: str,
    rna_mode: bool,
    mem_efficient: bool,
    single_cell: bool,
    cell_barcodes: Optional[str],
    min_cells: int,
    cellranger_path: Optional[str],
    barcode_whitelist: Optional[str],
    index_reads: Optional[str]
):
    """Assemble reads and identify viral contigs.
    
    This command performs:
    1. Assembly (SPAdes-only; rnaviral/sc/metaviral modes)
    2. Viral contig identification (Kaiju)
    3. Basic assembly validation
    
    Note: No preprocessing is performed - use preprocessed reads as input.
    Output: Assembled contigs and viral contigs ready for annotation.
    """
    
    # Validate mutually exclusive long-read inputs
    if nanopore and pacbio:
        click.echo("Error: Specify only one of --nanopore or --pacbio", err=True)
        sys.exit(1)

    # Validate input files
    input_files = [short_reads_1, short_reads_2, nanopore, pacbio, single_reads]
    input_files = [f for f in input_files if f is not None]
    
    if not input_files:
        click.echo("Error: No input read files provided", err=True)
        sys.exit(1)
    
    for file_path in input_files:
        if not Path(file_path).exists():
            click.echo(f"Error: Input file not found: {file_path}", err=True)
            sys.exit(1)
    
    # Create configuration
    config_dict = {
        "min_contig_length": min_contig_length,
        "viral_confidence_threshold": viral_confidence,
        "assembly_strategy": assembly_strategy,
        "rna_mode": rna_mode
    }
    
    # Load config file if provided
    if config and Path(config).exists():
        with open(config, 'r') as f:
            file_config = yaml.safe_load(f)
            config_dict.update(file_config)
    
    try:
        # Handle single-cell mode as pooled scRNA-seq: use R2 as single-end, enable RNA mode
        if single_cell:
            click.echo("Single-cell mode enabled - pooling R2 as single-end for scRNA-seq...")
            click.echo("Note: Using SPAdes --rnaviral without --sc for scRNA-seq")
            # Require paired-end inputs to extract cDNA reads (R2)
            if short_reads_1 and short_reads_2:
                # Set RNA mode and mark single-cell intent in config
                rna_mode = True
                config_dict["rna_mode"] = True
                config_dict["single_cell_mode"] = True
                # Treat R2 (cDNA) as single-end input; ignore barcode/UMI R1 for assembly
                single_reads = short_reads_2
                short_reads_1 = None
                short_reads_2 = None
            else:
                click.echo("Error: Single-cell mode requires paired-end reads (R1 and R2)", err=True)
                sys.exit(1)
        
        # Record long-read technology in config
        if nanopore:
            config_dict["long_read_tech"] = "nanopore"
        elif pacbio:
            config_dict["long_read_tech"] = "pacbio"

        # Initialize assembler
        assembler = ViralAssembler(
            output_dir=output_dir,
            threads=threads,
            memory=memory,
            config=config_dict,
            rna_mode=rna_mode or single_cell,  # Enable RNA mode for single-cell
            mem_efficient=mem_efficient
        )
        
        # Run assembly-only workflow (no preprocessing/annotation/quantification)
        # Map provided inputs to assembler's expected keys
        reads_dict: Dict[str, str] = {}
        if short_reads_1 and short_reads_2:
            reads_dict["short_1"] = short_reads_1
            reads_dict["short_2"] = short_reads_2
        if single_reads:
            reads_dict["single"] = single_reads
        # Map long reads
        long_reads_path = nanopore or pacbio
        if long_reads_path:
            reads_dict["long"] = long_reads_path

        # Store read paths in config for downstream quantification steps
        assembler.config["short_reads_1"] = reads_dict.get("short_1")
        assembler.config["short_reads_2"] = reads_dict.get("short_2")
        assembler.config["single_reads"] = reads_dict.get("single")
        assembler.config["long_reads"] = reads_dict.get("long")

        assembly_results = assembler._perform_assembly(reads_dict, reference)

        # Viral contig identification from assembled contigs
        viral_contig_results = assembler._identify_viral_contigs_efficiently(assembly_results)

        # Basic validation
        validation_results = assembler._validate_assemblies(
            viral_contig_results.get("viral_genomes", []),
            assembly_dir=assembler.output_dir / "01_assemblies"
        )

        # Prepare a compact results dict for reporting
        results = {
            "assembly_dir": str(assembler.output_dir),
            "viral_genomes": viral_contig_results.get("viral_genomes", []),
            "viral_contig_info": viral_contig_results,
            "validation_results": validation_results,
            "total_viral_contigs": sum(
                info.get("viral_contig_count", len(info.get("viral_contigs", [])))
                for info in viral_contig_results.get("viral_contig_info", {}).values()
            ),
        }
        
        # Print results summary
        click.echo("\n" + "="*50)
        click.echo("ASSEMBLY AND VIRAL IDENTIFICATION COMPLETED SUCCESSFULLY")
        click.echo("="*50)
        click.echo(f"Output directory: {results['assembly_dir']}")
        
        # Check if this was a reference-guided assembly
        if results.get('assembly_type') == 'reference_guided':
            click.echo(f"\nReference-Guided Assembly Results:")
            click.echo(f"Reference genome: {results.get('reference_genome', 'N/A')}")
            click.echo(f"Status: {results.get('status', 'Unknown')}")
            click.echo(f"Message: {results.get('message', 'N/A')}")
            
            if results.get('status') == 'completed':
                total_viral_contigs = results.get('total_viral_contigs', 0)
                reference_matching_contigs = results.get('reference_matching_contigs', [])
                num_reference_hits = results.get('num_reference_hits', 0)
                
                click.echo(f"Reference-matching contigs: {len(reference_matching_contigs)}")
                click.echo(f"Minimap2 hits to reference: {num_reference_hits}")
                click.echo(f"Total assembled contigs: {total_viral_contigs}")
                
                # Display assembly statistics
                stats = results.get('statistics', {})
                if stats:
                    click.echo(f"\nAssembly Statistics (with reference):")
                    for key, value in stats.items():
                        click.echo(f"  {key}: {value}")
                
                # Display CheckV results
                checkv_results = results.get('checkv_results', {})
                if checkv_results.get('status') == 'completed':
                    click.echo(f"\nCheckV Results:")
                    for key, value in checkv_results.items():
                        if key not in ['status', 'output_directory']:
                            click.echo(f"  {key}: {value}")
                
                # Display complete pipeline results
                if results.get('pipeline_status') == 'completed':
                    click.echo(f"\nComplete Pipeline Results:")
                    
                    # Viral classification results
                    viral_classifications = results.get('viral_classifications', {})
                    if viral_classifications:
                        click.echo(f"  Viral classification: Completed")
                        kaiju_summary = viral_classifications.get('kaiju_summary', {})
                        if kaiju_summary:
                            click.echo(f"    Classified contigs: {kaiju_summary.get('classified_contigs', 0)}")
                            click.echo(f"    Unclassified contigs: {kaiju_summary.get('unclassified_contigs', 0)}")
                    
                    # Gene prediction results
                    gene_predictions = results.get('gene_predictions', {})
                    if gene_predictions:
                        click.echo(f"  Gene prediction: Completed")
                        total_genes = gene_predictions.get('total_genes', 0)
                        annotated_genes = gene_predictions.get('total_annotated_proteins', 0)
                        click.echo(f"    Total genes predicted: {total_genes}")
                        click.echo(f"    Annotated genes: {annotated_genes}")
                    
                    # Quantification results
                    quantification = results.get('quantification', {})
                    if quantification:
                        click.echo(f"  Quantification: Completed")
                        quant_summary = quantification.get('summary', {})
                        if quant_summary:
                            click.echo(f"    Mapped reads: {quant_summary.get('total_mapped_reads', 0)}")
                            click.echo(f"    Average coverage: {quant_summary.get('average_coverage', 0):.2f}x")
                
                elif results.get('pipeline_status') == 'failed':
                    click.echo(f"\nPipeline Error: {results.get('error', 'Unknown error')}")
            else:
                click.echo(f"\nReference genome not detected in sample!")
                click.echo(f"Error: {results.get('error', 'Unknown error')}")
                click.echo(f"\nSuggestions:")
                click.echo(f"  - Check if the reference genome is correct")
                click.echo(f"  - Verify the reference virus is present in your sample")
                click.echo(f"  - Try running without --reference to see all viral contigs")
                return
        else:
            # Display viral contig information for normal assembly
            total_viral_contigs = results.get('total_viral_contigs', 0)
            stats = results.get('statistics') or {}
            if stats:
                click.echo(f"Total contigs: {stats.get('total_contigs', 0)}")
                total_len = stats.get('total_length')
                if total_len is not None:
                    click.echo(f"Total length: {int(total_len):,} bp")
                avg_len = stats.get('average_length', stats.get('average_contig_length'))
                if avg_len is not None:
                    click.echo(f"Average contig length: {float(avg_len):.0f} bp")
        
        if results['viral_genomes']:
            click.echo("\nViral contig files (contigs and scaffolds):")
            for genome_file in results['viral_genomes']:
                click.echo(f"  - {genome_file}")
            
            # Display viral classification information
            if 'viral_contig_info' in results:
                click.echo("\nViral Classification Summary:")
                for assembly_type, info in results['viral_contig_info'].items():
                    if 'classification_summary' in info:
                        summary = info['classification_summary']
                        click.echo(f"\n{assembly_type.title()} Assembly:")
                        click.echo(f"  Viral groups: {summary.get('viral_groups', {})}")
                        click.echo(f"  Average confidence: {summary.get('average_confidence', 0.0)}")
                        click.echo(f"  Average length: {summary.get('average_length', 0):.0f} bp")
                        click.echo(f"  Total hallmark genes: {summary.get('total_hallmark_genes', 0)}")
                        click.echo(f"  High confidence contigs: {summary.get('high_confidence_contigs', 0)}")
                        click.echo(f"  Medium confidence contigs: {summary.get('medium_confidence_contigs', 0)}")
                        click.echo(f"  Low confidence contigs: {summary.get('low_confidence_contigs', 0)}")
                        
                    # Display Kaiju classification results
                    if 'kaiju_classification' in info:
                        kaiju_results = info['kaiju_classification']
                        if kaiju_results.get('status') == 'completed':
                            classifications = kaiju_results.get('classifications', {})
                            click.echo(f"  Kaiju classification: {len(classifications)} contigs classified")
                            click.echo(f"  Kaiju summary: {kaiju_results.get('summary_file', 'N/A')}")
                        else:
                            click.echo(f"  Kaiju classification: Failed - {kaiju_results.get('error', 'Unknown error')}")
                    
                    # Display viral quantification results
                    if 'quantification' in info:
                        quant_results = info['quantification']
                        if quant_results.get('status') == 'completed':
                            quant_data = quant_results.get('quantification_results', {})
                            click.echo(f"  Viral quantification: {len(quant_data)} contigs quantified")
                            click.echo(f"  Abundance file: {quant_results.get('bam_file', 'N/A').replace('.bam', '_abundance.tsv')}")
                        else:
                            click.echo(f"  Viral quantification: Failed - {quant_results.get('error', 'Unknown error')}")
            
            # Display validation results
            if 'validation_results' in results:
                click.echo("\nValidation Results:")
                
                # Assembly quality evaluation removed
                
                # Display individual genome validation results
                if 'validation_results' in results:
                    for genome_file, validation in results['validation_results'].items():
                        # Skip any special validation entries
                            
                        genome_name = Path(genome_file).stem
                        click.echo(f"\n{genome_name}:")
                        
                        # Display CheckV results
                        checkv_results = validation.get('checkv', {})
                        if checkv_results.get('status') == 'completed':
                            click.echo("  CheckV Quality Assessment:")
                            click.echo(f"    Total contigs analyzed: {checkv_results.get('total_contigs', 0)}")
                            click.echo(f"    Complete contigs (≥90%): {checkv_results.get('complete_contigs', 0)}")
                            click.echo(f"    High quality (≥50%): {checkv_results.get('high_quality_contigs', 0)}")
                            click.echo(f"    Medium quality (30-50%): {checkv_results.get('medium_quality_contigs', 0)}")
                            click.echo(f"    Low quality (10-30%): {checkv_results.get('low_quality_contigs', 0)}")
                            click.echo(f"    Undetermined (<10%): {checkv_results.get('undetermined_contigs', 0)}")
                            click.echo(f"    Average completeness: {checkv_results.get('avg_completeness', 0):.1f}%")
                            click.echo(f"    Average contamination: {checkv_results.get('avg_contamination', 0):.1f}%")
                        elif checkv_results.get('status') == 'skipped':
                            click.echo(f"  CheckV: Skipped - {checkv_results.get('reason', 'Unknown reason')}")
                        else:
                            click.echo(f"  CheckV: Failed - {checkv_results.get('error', 'Unknown error')}")
                        
                        # Display assembly statistics
                        stats = validation.get('statistics', {})
                        if stats:
                            click.echo("  Assembly Statistics:")
                            for key, value in stats.items():
                                click.echo(f"    {key}: {value}")
                
                # Display gene prediction results
                if 'gene_predictions' in results and results['gene_predictions'].get('status') == 'completed':
                    click.echo("\nGene Prediction Results:")
                    gene_pred = results['gene_predictions']
                    summary = gene_pred.get('summary', {})
                    click.echo(f"  Total contigs analyzed: {summary.get('total_contigs', 0)}")
                    click.echo(f"  Total genes predicted: {summary.get('total_genes', 0)}")
                    click.echo(f"  Average genes per contig: {summary.get('avg_genes_per_contig', 0)}")
                    click.echo(f"  Average gene length: {summary.get('avg_gene_length', 0)} bp")
                    
                    # VOG annotation results
                    if summary.get('vog_available', False):
                        click.echo(f"\nVOG Functional Annotation:")
                        click.echo(f"  Annotated proteins: {summary.get('vog_annotated_proteins', 0)}")
                        click.echo(f"  Viral classifications: {summary.get('vog_viral_classifications', 0)}")
                    else:
                        click.echo(f"\nVOG Functional Annotation: Not available")
                    
                    click.echo(f"  Gene prediction directory: {gene_pred.get('output_directory', 'N/A')}")
                
                # Display analysis recommendations
                if 'scaffolds' in results['viral_contig_info'] and 'contigs' in results['viral_contig_info']:
                    click.echo("\n" + "="*60)
                    click.echo("ANALYSIS COMPLETE - CHOOSE YOUR RESULTS")
                    click.echo("="*60)
                    click.echo("Both contigs and scaffolds have been analyzed with Kaiju.")
                    click.echo("Choose the best results based on your downstream analysis needs:")
                    click.echo("\nCONTIGS:")
                    click.echo("  - No gaps (N's) - better for gene-based analysis")
                    click.echo("  - More reliable for Kaiju classification")
                    click.echo("  - May be fragmented - shorter sequences")
                    click.echo("\nSCAFFOLDS:")
                    click.echo("  - Longer sequences - more complete genomes")
                    click.echo("  - Better for downstream analysis")
                    click.echo("  - May contain gaps (N's) - can affect gene prediction")
                    click.echo("\nRECOMMENDATION:")
                    click.echo("  • For gene-based analysis: Use CONTIGS")
                    click.echo("  • For genome completeness: Use SCAFFOLDS")
                    click.echo("  • For publication: Use SCAFFOLDS (if gap percentage < 5%)")
                    click.echo("="*60)
        
        click.echo("\nDetailed validation results available in output directory.")
        
    except Exception as e:
        logger.error(f"Assembly failed: {e}")
        click.echo(f"Error: Assembly failed - {e}", err=True)
        sys.exit(1)


@main.command()
@click.option('--short-reads-1', help='Path to first mate of paired-end reads')
@click.option('--short-reads-2', help='Path to second mate of paired-end reads')
@click.option('--nanopore', help='Path to ONT Nanopore long reads')
@click.option('--pacbio', help='Path to PacBio long reads')
@click.option('--single-reads', help='Path to single-end reads')
@click.option('--reference', help='Optional reference genome for guided assembly')
@click.option('--output-dir', '-o', required=True, help='Output directory')
@click.option('--threads', '-t', default=8, help='Number of threads')
@click.option('--memory', default='16G', help='Memory allocation')
@click.option('--config', help='Configuration file path')
@click.option('--min-contig-length', default=1000, help='Minimum contig length')
@click.option('--viral-confidence', default=0.8, help='Viral sequence confidence threshold')
@click.option('--assembly-strategy', 
              type=click.Choice(['hybrid', 'short_only', 'long_only']),
              default='hybrid', help='Assembly strategy')
@click.option('--rna-mode', is_flag=True, help='Enable RNA-specific assembly parameters')
@click.option('--mem-efficient', '-m', is_flag=True, help='Enable memory-efficient mode with read subsampling for large datasets')
@click.option('--single-cell', is_flag=True, help='Enable single-cell sequencing mode')
@click.option('--cell-barcodes', help='Path to cell barcodes file (if not 10X format)')
@click.option('--min-cells', default=100, help='Minimum number of cells to process')
@click.option('--cellranger-path', help='Path to Cell Ranger installation')
@click.option('--barcode-whitelist', help='Path to barcode whitelist file')
@click.option('--index-reads', help='Path to index reads (I1) for single-cell data')
def analyse(
    short_reads_1: Optional[str],
    short_reads_2: Optional[str],
    nanopore: Optional[str],
    pacbio: Optional[str],
    single_reads: Optional[str],
    reference: Optional[str],
    output_dir: str,
    threads: int,
    memory: str,
    config: Optional[str],
    min_contig_length: int,
    viral_confidence: float,
    assembly_strategy: str,
    rna_mode: bool,
    mem_efficient: bool,
    single_cell: bool,
    cell_barcodes: Optional[str],
    min_cells: int,
    cellranger_path: Optional[str],
    barcode_whitelist: Optional[str],
    index_reads: Optional[str]
):
    """Run complete viral genome analysis pipeline from sequencing reads.
    
    This command performs the full analysis workflow:
    1. Preprocessing (quality control, trimming)
    2. Assembly (SPAdes)
    3. Viral contig identification (Kaiju)
    4. Gene prediction and annotation (Prodigal + VOG)
    5. Validation and quantification (CheckV + BWA)
    """
    
    # Validate mutually exclusive long-read inputs
    if nanopore and pacbio:
        click.echo("Error: Specify only one of --nanopore or --pacbio", err=True)
        sys.exit(1)

    # Validate input files
    input_files = [short_reads_1, short_reads_2, nanopore, pacbio, single_reads]
    input_files = [f for f in input_files if f is not None]
    
    if not input_files:
        click.echo("Error: At least one input file must be provided", err=True)
        sys.exit(1)
    
    for file_path in input_files:
        if not os.path.exists(file_path):
            click.echo(f"Error: Input file not found: {file_path}", err=True)
            sys.exit(1)
    
    # Load configuration
    config_dict = {
        "min_contig_length": min_contig_length,
        "viral_confidence_threshold": viral_confidence,
        "assembly_strategy": assembly_strategy,
        "rna_mode": rna_mode
    }
    
    if config and os.path.exists(config):
        with open(config, 'r') as f:
            file_config = yaml.safe_load(f)
            config_dict.update(file_config)
    
    try:
        # Handle single-cell mode as pooled scRNA-seq: use R2 as single-end, enable RNA mode
        if single_cell:
            click.echo("Single-cell mode enabled - pooling R2 as single-end for scRNA-seq...")
            click.echo("Note: Using SPAdes --rnaviral without --sc for scRNA-seq")
            if short_reads_1 and short_reads_2:
                rna_mode = True
                config_dict["rna_mode"] = True
                config_dict["single_cell_mode"] = True
                single_reads = short_reads_2
                short_reads_1 = None
                short_reads_2 = None
            else:
                click.echo("Error: Single-cell mode requires paired-end reads (R1 and R2)", err=True)
                sys.exit(1)
        
        # Record long-read technology in config
        if nanopore:
            config_dict["long_read_tech"] = "nanopore"
        elif pacbio:
            config_dict["long_read_tech"] = "pacbio"

        # Initialize assembler
        assembler = ViralAssembler(
            output_dir=output_dir,
            threads=threads,
            memory=memory,
            config=config_dict,
            rna_mode=rna_mode or single_cell,  # Enable RNA mode for single-cell
            mem_efficient=mem_efficient
        )
        
        # Run complete analysis pipeline
        long_reads_path = nanopore or pacbio
        results = assembler.assemble(
            short_reads_1=short_reads_1,
            short_reads_2=short_reads_2,
            long_reads=long_reads_path,
            single_reads=single_reads,
            reference=reference
        )
        
        # Check if reference-guided assembly failed
        if results.get('reference_guided_failed', False):
            click.echo("\n" + "="*50)
            click.echo("REFERENCE-GUIDED ASSEMBLY FAILED")
            click.echo("="*50)
            click.echo(f"Output directory: {results['assembly_dir']}")
            click.echo(f"Error: {results.get('error', 'Unknown error')}")
            click.echo(f"Message: {results.get('message', 'No additional information')}")
            click.echo("\nThis means no contigs similar to the reference genome were found.")
            click.echo("This could indicate:")
            click.echo("  - The reference genome is not present in the sample")
            click.echo("  - The reference genome is too divergent from the sample")
            click.echo("  - The assembly quality was insufficient")
            click.echo("  - The filtering parameters were too strict")
            return
        
        # Print results summary for successful analysis
        click.echo("\n" + "="*50)
        click.echo("ANALYSIS COMPLETED SUCCESSFULLY")
        click.echo("="*50)
        click.echo(f"Output directory: {results['assembly_dir']}")
        # Display viral contig information
        total_viral_contigs = results.get('total_viral_contigs', 0)
        click.echo(f"Viral contigs found: {total_viral_contigs}")
        click.echo(f"Total contigs: {results['statistics']['total_contigs']}")
        click.echo(f"Total length: {results['statistics']['total_length']:,} bp")
        click.echo(f"Average contig length: {results['statistics']['average_contig_length']:,} bp")
        
        click.echo(f"\nViral contig files (contigs and scaffolds):")
        for genome_file in results['viral_genomes']:
            click.echo(f"  - {genome_file}")
        
        # Display organized output structure
        click.echo(f"\nOrganized output structure:")
        click.echo(f"  01_assemblies/          - Assembly outputs (contigs, scaffolds)")
        click.echo(f"  02_viral_contigs/       - Viral-specific contig files")
        click.echo(f"  03_classifications/     - Classification results (Kaiju)")
        click.echo(f"  04_quality_assessment/  - Quality metrics (CheckV)")
        click.echo(f"  05_gene_predictions/    - Gene prediction and annotation")
        click.echo(f"  06_quantification/      - Read mapping and abundance")
        
        # Display Kaiju classification results
        viral_contig_info = results.get('viral_contig_info', {})
        if viral_contig_info:
            click.echo(f"\nViral Classification Summary:")
            for assembly_type, info in viral_contig_info.items():
                if 'kaiju_classification' in info:
                    kaiju_result = info['kaiju_classification']
                    if kaiju_result.get('status') == 'completed':
                        summary_file = kaiju_result.get('summary_file', 'N/A')
                        click.echo(f"  {assembly_type}: {summary_file}")
                        
                        # Try to read and display classification results
                        if summary_file != 'N/A' and os.path.exists(summary_file):
                            try:
                                with open(summary_file, 'r') as f:
                                    lines = f.readlines()
                                    if len(lines) > 1:  # Skip header
                                        click.echo(f"    Identified viruses:")
                                        for line in lines[1:]:  # Skip header line
                                            parts = line.strip().split('\t')
                                            if len(parts) >= 4:
                                                contig_id, taxon_id, taxon_name, classification = parts[:4]
                                                click.echo(f"      - {contig_id}: {classification}")
                            except Exception as e:
                                click.echo(f"    (Could not read classification details: {e})")
        
        # Display validation results
        if 'validation_results' in results:
            click.echo("\nValidation Results:")
            for genome_file, validation in results['validation_results'].items():
                genome_name = Path(genome_file).stem
                click.echo(f"\n{genome_name}:")
                
                # Display assembly statistics
                stats = validation.get('statistics', {})
                if stats:
                    click.echo(f"  Assembly statistics:")
                    for key, value in stats.items():
                        click.echo(f"    {key}: {value}")
        
        # Display gene prediction results
        if 'gene_prediction_results' in results:
            gene_results = results['gene_prediction_results']
            click.echo(f"\nGene Prediction Results:")
            click.echo(f"  Total contigs analyzed: {gene_results.get('total_contigs', 0)}")
            click.echo(f"  Total genes predicted: {gene_results.get('total_genes', 0)}")
            click.echo(f"  Average genes per contig: {gene_results.get('avg_genes_per_contig', 0):.1f}")
            click.echo(f"  Average gene length: {gene_results.get('avg_gene_length', 0):.0f} bp")
        
        # Display VOG annotation results
        if 'vog_annotations' in results:
            vog_results = results['vog_annotations']
            click.echo(f"\nVOG Functional Annotation:")
            click.echo(f"  Annotated proteins: {vog_results.get('total_annotated_proteins', 0)}")
            
            # Display VOG classifications summary instead of full dictionary
            viral_classifications = vog_results.get('viral_classifications', {})
            if isinstance(viral_classifications, dict):
                unique_vogs = len(set(viral_classifications.keys()))
                click.echo(f"  Unique VOG families: {unique_vogs}")
                if unique_vogs > 0:
                    # Show top 5 most common VOG families
                    from collections import Counter
                    vog_counts = Counter(viral_classifications.values())
                    top_vogs = vog_counts.most_common(5)
                    click.echo("  Top VOG functional categories:")
                    for category, count in top_vogs:
                        click.echo(f"    {category}: {count} genes")
            else:
                click.echo(f"  Viral classifications: {viral_classifications}")
            
            click.echo(f"  Gene prediction directory: {vog_results.get('output_directory', 'N/A')}")
        
        click.echo(f"\nDetailed validation results available in output directory.")
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        click.echo(f"Error: Analysis failed - {e}", err=True)
        sys.exit(1)


@main.command()
@click.option('--viral-contigs', '-v', required=True, help='Path to viral contigs file')
@click.option('--output-dir', '-o', required=True, help='Output directory for annotation results')
@click.option('--threads', '-t', default=8, help='Number of threads')
@click.option('--config', help='Configuration file path')
@click.option('--kaiju-results', help='Path to Kaiju classification results directory (optional)')
def annotate(viral_contigs: str, output_dir: str, threads: int, config: Optional[str], kaiju_results: Optional[str]):
    """Annotate viral contigs with gene prediction and functional annotation.
    
    This command performs:
    1. Gene prediction (Prodigal)
    2. Functional annotation (VOG)
    3. Protein sequence extraction
    4. Species-specific protein grouping (if Kaiju results provided)
    
    Input: Viral contigs file (FASTA)
    Optional: Kaiju classification results directory
    Output: Annotated genes and proteins, optionally grouped by species
    """
    
    if not os.path.exists(viral_contigs):
        click.echo(f"Error: Viral contigs file not found: {viral_contigs}", err=True)
        sys.exit(1)
    
    # Load configuration
    config_dict = {}
    if config and os.path.exists(config):
        with open(config, 'r') as f:
            file_config = yaml.safe_load(f)
            config_dict.update(file_config)
    
    try:
        from .core.gene_predictor import ViralGenePredictor
        
        # Initialize gene predictor
        gene_predictor = ViralGenePredictor()
        
        # Load classification data
        classification_data = {}
        
        if kaiju_results and os.path.exists(kaiju_results):
            # Load existing Kaiju results
            kaiju_summary_file = os.path.join(kaiju_results, 'kaiju_summary.tsv')
            if os.path.exists(kaiju_summary_file):
                logger.info(f"Loading Kaiju classification results from {kaiju_results}")
                classifications = {}
                
                # Parse Kaiju summary file
                with open(kaiju_summary_file, 'r') as f:
                    lines = f.readlines()
                    if len(lines) > 1:  # Skip header
                        for line in lines[1:]:
                            parts = line.strip().split('\t')
                            if len(parts) >= 4:
                                contig_id = parts[0]
                                taxon_id = parts[1]
                                taxon_name = parts[2]
                                classification = parts[3]
                                
                                classifications[contig_id] = {
                                    'taxon_id': taxon_id,
                                    'taxon_name': taxon_name,
                                    'classification': classification
                                }
                
                classification_data = {
                    'kaiju_contigs': {
                        'method': 'kaiju',
                        'classification_method': 'kaiju_contigs_mode',
                        'status': 'completed',
                        'classifications': classifications,
                        'summary_file': kaiju_summary_file
                    }
                }
                logger.info(f"Loaded {len(classifications)} Kaiju classifications")
            else:
                logger.warning(f"Kaiju summary file not found: {kaiju_summary_file}")
                classification_data = {
                    'kaiju_contigs': {
                        'method': 'kaiju',
                        'classification_method': 'kaiju_contigs_mode',
                        'status': 'completed',
                        'classifications': {},
                        'summary_file': ''
                    }
                }
        else:
            # Create mock classification data for Kaiju mode
            classification_data = {
                'kaiju_contigs': {
                    'method': 'kaiju',
                    'classification_method': 'kaiju_contigs_mode',
                    'status': 'completed',
                    'classifications': {},
                    'summary_file': ''
                }
            }
        
        # Extract actual classifications for gene prediction
        viral_classifications = {}
        for classification_group in classification_data.values():
            if 'classifications' in classification_group:
                viral_classifications.update(classification_group['classifications'])
        
        # Run gene prediction and annotation
        results = gene_predictor.predict_genes_comprehensive(
            viral_contigs,
            viral_classifications,
            output_dir
        )
        
        # Print results summary
        click.echo("\n" + "="*50)
        click.echo("ANNOTATION COMPLETED SUCCESSFULLY")
        click.echo("="*50)
        click.echo(f"Output directory: {output_dir}")
        
        if results.get('status') == 'completed':
            summary = results.get('summary', {})
            click.echo(f"Total contigs analyzed: {summary.get('total_contigs', 0)}")
            click.echo(f"Total genes predicted: {summary.get('total_genes', 0)}")
            click.echo(f"Average genes per contig: {summary.get('avg_genes_per_contig', 0):.1f}")
            click.echo(f"Average gene length: {summary.get('avg_gene_length', 0):.0f} bp")
            
            vog_results = results.get('vog_annotations', {})
            annotated_proteins = vog_results.get('total_annotated_proteins', 0)
            viral_classifications = results.get('vog_classifications', vog_results.get('viral_classifications', {}))
            
            click.echo(f"Annotated proteins: {annotated_proteins}")
            
            # Display Kaiju taxonomic classification
            kaiju_classifications = results.get('kaiju_classifications', {})
            if kaiju_classifications:
                from collections import Counter
                kaiju_counts = Counter(kaiju_classifications.values())
                click.echo(f"Kaiju taxonomic classifications:")
                for species, count in kaiju_counts.most_common():
                    click.echo(f"  {species}: {count} contigs")
            else:
                click.echo("No Kaiju taxonomic classifications found")
            
            # Display VOG classifications summary instead of full dictionary
            if isinstance(viral_classifications, dict):
                unique_vogs = len(set(viral_classifications.keys()))
                click.echo(f"Unique VOG families: {unique_vogs}")
                if unique_vogs > 0:
                    # Show top 10 most common VOG families
                    from collections import Counter
                    vog_counts = Counter(viral_classifications.values())
                    top_vogs = vog_counts.most_common(10)
                    click.echo("Top VOG functional categories:")
                    for category, count in top_vogs:
                        click.echo(f"  {category}: {count} genes")
            else:
                click.echo(f"Viral classifications: {viral_classifications}")
        else:
            click.echo(f"Annotation failed: {results.get('error', 'Unknown error')}")
        
        click.echo(f"\nDetailed results available in: {output_dir}")
        
    except Exception as e:
        logger.error(f"Annotation failed: {e}")
        click.echo(f"Error: Annotation failed - {e}", err=True)
        sys.exit(1)

@main.command()
@click.option('--genome', '-g', required=True, help='Genome file to validate')
@click.option('--output-dir', '-o', help='Output directory for validation results')
@click.option('--reference', help='Reference genome for comparison')
@click.option('--assembly-dir', help='Optional assembly directory for comparative statistics')
@click.option('--threads', '-t', default=8, help='Number of threads')
def validate(genome: str, output_dir: Optional[str], reference: Optional[str], assembly_dir: Optional[str], threads: int):
    """Validate and assess quality of assembled viral genomes."""
    
    if not Path(genome).exists():
        click.echo(f"Error: Genome file not found: {genome}", err=True)
        sys.exit(1)
    
    try:
        # Initialize validator
        validator = AssemblyValidator(threads=threads)
        
        # Assembly quality evaluation removed
        
        # Run validation tools
        if output_dir:
            # Create output directory if specified
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            checkv_results = validator.run_checkv(genome, output_dir=f"{output_dir}/checkv")
        else:
            checkv_results = validator.run_checkv(genome)
        stats = validator.calculate_genome_statistics(genome)
        
        # Print results
        click.echo("\n" + "="*50)
        click.echo("GENOME VALIDATION RESULTS")
        click.echo("="*50)
        
        click.echo("\nBasic Statistics:")
        for key, value in stats.items():
            click.echo(f"  {key}: {value}")
        
        # Display CheckV results
        if checkv_results.get("status") == "completed":
            click.echo("\nCheckV Quality Assessment:")
            click.echo(f"  Total contigs analyzed: {checkv_results.get('total_contigs', 0)}")
            click.echo(f"  Complete contigs (≥90%): {checkv_results.get('complete_contigs', 0)}")
            click.echo(f"  High quality (≥50%): {checkv_results.get('high_quality_contigs', 0)}")
            click.echo(f"  Medium quality (30-50%): {checkv_results.get('medium_quality_contigs', 0)}")
            click.echo(f"  Low quality (10-30%): {checkv_results.get('low_quality_contigs', 0)}")
            click.echo(f"  Undetermined (<10%): {checkv_results.get('undetermined_contigs', 0)}")
            click.echo(f"  Average completeness: {checkv_results.get('avg_completeness', 0):.1f}%")
            click.echo(f"  Average contamination: {checkv_results.get('avg_contamination', 0):.1f}%")
        elif checkv_results.get("status") == "skipped":
            click.echo(f"\nCheckV: Skipped - {checkv_results.get('reason', 'Unknown reason')}")
        else:
            click.echo(f"\nCheckV: Failed - {checkv_results.get('error', 'Unknown error')}")
        
        click.echo("\nAssembly Statistics:")
        for key, value in stats.items():
            click.echo(f"  {key}: {value}")
        
        # Generate report if output directory specified
        if output_dir:
            report_file = validator.generate_quality_report([genome], output_file=f"{output_dir}/quality_report.html")
            click.echo(f"\nDetailed report saved to: {report_file}")
        
    except Exception as e:
        logger.error(f"Validation failed: {e}")
        click.echo(f"Error: Validation failed - {e}", err=True)
        sys.exit(1)


@main.command()
@click.option('--viral-contigs', '-v', required=True, help='Path to viral contigs file')
@click.option('--short-reads-1', help='Path to first mate of paired-end reads')
@click.option('--short-reads-2', help='Path to second mate of paired-end reads')
@click.option('--single-reads', help='Path to single-end reads')
@click.option('--nanopore', help='Path to ONT Nanopore long reads')
@click.option('--pacbio', help='Path to PacBio long reads')
@click.option('--output-dir', '-o', required=True, help='Output directory for quantification results')
@click.option('--threads', '-t', default=8, help='Number of threads')
def quantify(
    viral_contigs: str,
    short_reads_1: Optional[str],
    short_reads_2: Optional[str],
    single_reads: Optional[str],
    nanopore: Optional[str],
    pacbio: Optional[str],
    output_dir: str,
    threads: int
):
    """Quantify viral contigs by mapping reads back to them.
    
    This command performs:
    1. Read mapping (BWA for short reads, minimap2 for long reads)
    2. Coverage calculation
    3. Abundance estimation
    
    Input: Viral contigs file + sequencing reads
    Output: Quantification statistics and BAM files
    """
    
    if not Path(viral_contigs).exists():
        click.echo(f"Error: Viral contigs file not found: {viral_contigs}", err=True)
        sys.exit(1)
    
    # Validate that at least one read file is provided
    if nanopore and pacbio:
        click.echo("Error: Specify only one of --nanopore or --pacbio", err=True)
        sys.exit(1)

    read_files = [short_reads_1, short_reads_2, single_reads, nanopore, pacbio]
    read_files = [f for f in read_files if f is not None]
    
    if not read_files:
        click.echo("Error: At least one read file must be provided", err=True)
        sys.exit(1)
    
    # Validate read files exist
    for file_path in read_files:
        if not Path(file_path).exists():
            click.echo(f"Error: Read file not found: {file_path}", err=True)
            sys.exit(1)
    
    try:
        from .core.viral_identifier import ViralIdentifier
        
        # Initialize viral identifier
        viral_identifier = ViralIdentifier()
        
        # Run quantification
        long_reads_path = nanopore or pacbio
        results = viral_identifier.quantify_viral_contigs(
            contigs_file=viral_contigs,
            reads_1=short_reads_1,
            reads_2=short_reads_2,
            single_reads=single_reads,
            long_reads=long_reads_path,
            output_dir=output_dir
        )
        
        # Print results summary
        click.echo("\n" + "="*50)
        click.echo("VIRAL CONTIG QUANTIFICATION COMPLETED")
        click.echo("="*50)
        click.echo(f"Output directory: {output_dir}")
        
        if results.get("status") == "completed":
            quant_data = results.get("quantification_results", {})
            click.echo(f"Contigs quantified: {len(quant_data)}")
            
            if quant_data:
                # Calculate summary statistics
                total_coverage = sum(contig.get("total_coverage", 0) for contig in quant_data.values())
                avg_coverage = total_coverage / len(quant_data) if quant_data else 0
                max_coverage = max(contig.get("total_coverage", 0) for contig in quant_data.values())
                
                click.echo(f"Total coverage: {total_coverage:,.0f} reads")
                click.echo(f"Average coverage per contig: {avg_coverage:,.1f} reads")
                click.echo(f"Maximum coverage: {max_coverage:,.0f} reads")
                
                # Show top 5 most abundant contigs
                sorted_contigs = sorted(quant_data.items(), 
                                      key=lambda x: x[1].get("total_coverage", 0), 
                                      reverse=True)
                
                click.echo("\nTop 5 most abundant contigs:")
                for i, (contig_id, data) in enumerate(sorted_contigs[:5]):
                    coverage = data.get("total_coverage", 0)
                    abundance = data.get("abundance", 0)
                    click.echo(f"  {i+1}. {contig_id}: {coverage:,.0f} reads ({abundance:.2f}%)")
                
                # Show output files
                click.echo(f"\nOutput files:")
                click.echo(f"  BAM file: {results.get('bam_file', 'N/A')}")
                click.echo(f"  Abundance file: {results.get('bam_file', 'N/A').replace('.bam', '_abundance.tsv')}")
            else:
                click.echo("No contigs were quantified")
        else:
            click.echo(f"Quantification failed: {results.get('error', 'Unknown error')}")
        
        click.echo(f"\nDetailed results available in: {output_dir}")
        
    except Exception as e:
        logger.error(f"Quantification failed: {e}")
        click.echo(f"Error: Quantification failed - {e}", err=True)
        sys.exit(1)


@main.command()
@click.option('--viral-contigs', '-v', required=True, help='Path to viral contigs file')
@click.option('--output-dir', '-o', help='Output directory for classification results')
@click.option('--threads', '-t', default=8, help='Number of threads')
def classify(viral_contigs: str, output_dir: Optional[str], threads: int):
    """Classify viral contigs using Kaiju nucleotide-based classification.
    
    This command performs:
    1. Kaiju classification against viral genome database
    2. Taxonomic classification and confidence scoring
    
    Input: Viral contigs file (FASTA)
    Output: Classification results and summary file
    """
    
    if not Path(viral_contigs).exists():
        click.echo(f"Error: Viral contigs file not found: {viral_contigs}", err=True)
        sys.exit(1)
    
    try:
        from .core.viral_identifier import ViralIdentifier
        
        # Initialize viral identifier
        viral_identifier = ViralIdentifier()
        
        # Run Kaiju classification
        results = viral_identifier.classify_viral_contigs(
            contigs_file=viral_contigs,
            output_dir=output_dir
        )
        
        # Print results summary
        click.echo("\n" + "="*50)
        click.echo("VIRAL CONTIG CLASSIFICATION COMPLETED")
        click.echo("="*50)
        
        if results.get("status") == "completed":
            classifications = results.get("classifications", {})
            click.echo(f"Contigs classified: {len(classifications)}")
            click.echo(f"Summary file: {results.get('summary_file', 'N/A')}")
            
            if classifications:
                # Show sample classifications
                click.echo("\nSample classifications:")
                for i, (contig_id, classification) in enumerate(list(classifications.items())[:5]):
                    # For Kaiju, use direct classification (no consensus needed)
                    classification_name = classification.get('classification', 'Unknown')
                    click.echo(f"  {contig_id}: {classification_name}")
                
                if len(classifications) > 5:
                    click.echo(f"  ... and {len(classifications) - 5} more")
                
                # Show classification distribution
                classification_counts = {}
                for classification in classifications.values():
                    classification_name = classification.get('classification', 'Unknown')
                    classification_counts[classification_name] = classification_counts.get(classification_name, 0) + 1
                
                click.echo(f"\nClassification distribution:")
                for classification_name, count in classification_counts.items():
                    click.echo(f"  {classification_name}: {count} contigs")
            else:
                click.echo("No contigs were classified")
        else:
            click.echo(f"Classification failed: {results.get('error', 'Unknown error')}")
                
        if output_dir:
            click.echo(f"\nDetailed results available in: {output_dir}")
        else:
            click.echo(f"\nDetailed results available in: {results.get('summary_file', 'N/A')}")
        
    except Exception as e:
        logger.error(f"Classification failed: {e}")
        click.echo(f"Error: Classification failed - {e}", err=True)
        sys.exit(1)


@main.command()
@click.option('--viral-sequences', required=True, help='File containing viral sequences')
@click.option('--non-viral-sequences', required=True, help='File containing non-viral sequences')
@click.option('--output-model', required=True, help='Output path for trained model')
@click.option('--test-size', default=0.2, help='Fraction of data for testing')
def train_model(viral_sequences: str, non_viral_sequences: str, output_model: str, test_size: float):
    """Train a viral sequence identification model.
    
    Note: This feature is under development and not yet fully implemented.
    This feature will be available in future versions of Virall.
    """
    
    click.echo("Note: The train-model command is under development and not yet fully implemented.")
    click.echo("This feature will be available in future versions of Virall.")
    
    if not Path(viral_sequences).exists():
        click.echo(f"Error: Viral sequences file not found: {viral_sequences}", err=True)
        sys.exit(1)
    
    if not Path(non_viral_sequences).exists():
        click.echo(f"Error: Non-viral sequences file not found: {non_viral_sequences}", err=True)
        sys.exit(1)
    
    try:
        # Load sequences
        from Bio import SeqIO
        
        viral_seqs = [str(record.seq) for record in SeqIO.parse(viral_sequences, "fasta")]
        non_viral_seqs = [str(record.seq) for record in SeqIO.parse(non_viral_sequences, "fasta")]
        
        click.echo(f"Loaded {len(viral_seqs)} viral sequences")
        click.echo(f"Loaded {len(non_viral_seqs)} non-viral sequences")
        
        # Initialize identifier and train model
        identifier = ViralIdentifier()
        results = identifier.train_model(viral_seqs, non_viral_seqs, test_size)
        
        # Save model
        identifier.save_model(output_model)
        
        # Print results
        click.echo("\n" + "="*50)
        click.echo("MODEL TRAINING COMPLETED")
        click.echo("="*50)
        click.echo(f"Accuracy: {results['accuracy']:.3f}")
        click.echo(f"Model saved to: {output_model}")
        
    except Exception as e:
        logger.error(f"Model training failed: {e}")
        click.echo(f"Error: Model training failed - {e}", err=True)
        sys.exit(1)


@main.command()
@click.option('--reads-1', help='First mate of paired-end reads')
@click.option('--reads-2', help='Second mate of paired-end reads')
@click.option('--single-reads', help='Single-end reads file')
@click.option('--nanopore', help='ONT Nanopore long reads file')
@click.option('--pacbio', help='PacBio long reads file')
@click.option('--output-dir', '-o', required=True, help='Output directory')
@click.option('--threads', '-t', default=8, help='Number of threads')
def preprocess(reads_1: Optional[str], reads_2: Optional[str], single_reads: Optional[str], 
               nanopore: Optional[str], pacbio: Optional[str], output_dir: str, threads: int):
    """Preprocess sequencing reads (quality control, trimming, error correction).
    
    This command performs:
    1. Quality control (FastQC)
    2. Adapter trimming (fastp with automatic adapter detection)
    3. Error correction (SPAdes)
    
    Supports: Single-end, paired-end, and long reads
    """
    
    # Validate input files
    if nanopore and pacbio:
        click.echo("Error: Specify only one of --nanopore or --pacbio", err=True)
        sys.exit(1)

    input_files = [reads_1, reads_2, single_reads, nanopore, pacbio]
    input_files = [f for f in input_files if f is not None]
    
    if not input_files:
        click.echo("Error: At least one read file must be provided", err=True)
        sys.exit(1)
    
    for file_path in input_files:
        if not Path(file_path).exists():
            click.echo(f"Error: Read file not found: {file_path}", err=True)
            sys.exit(1)
    
    try:
        # Initialize preprocessor
        preprocessor = Preprocessor(threads=threads)
        
        # Create output directory
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        processed_files = []
        
        # Process paired-end reads
        if reads_1 and reads_2:
            click.echo("Processing paired-end reads...")
            processed_1, processed_2 = preprocessor.process_paired_reads(reads_1, reads_2)
            
            # Copy to output directory
            import shutil
            output_1 = Path(output_dir) / f"trimmed_{Path(reads_1).name}"
            output_2 = Path(output_dir) / f"trimmed_{Path(reads_2).name}"
            shutil.copy2(processed_1, output_1)
            shutil.copy2(processed_2, output_2)
            processed_files.extend([str(output_1), str(output_2)])
            
            click.echo(f"Paired-end reads processed:")
            click.echo(f"  R1: {output_1}")
            click.echo(f"  R2: {output_2}")
        
        # Process single-end reads
        if single_reads:
            click.echo("Processing single-end reads...")
            processed_single = preprocessor.process_single_reads(single_reads)
            
            # Copy to output directory
            import shutil
            output_single = Path(output_dir) / f"trimmed_{Path(single_reads).name}"
            shutil.copy2(processed_single, output_single)
            processed_files.append(str(output_single))
            
            click.echo(f"Single-end reads processed: {output_single}")
        
        # Process long reads
        long_reads = nanopore or pacbio
        if long_reads:
            click.echo("Processing long reads...")
            processed_long = preprocessor.process_long_reads(long_reads)

            # Copy to output directory
            import shutil
            output_long = Path(output_dir) / f"trimmed_{Path(long_reads).name}"
            shutil.copy2(processed_long, output_long)
            processed_files.append(str(output_long))

            click.echo(f"Long reads processed: {output_long}")
        
        click.echo(f"\nPreprocessing completed successfully!")
        click.echo(f"Processed files saved to: {output_dir}")
        
        # Generate quality report and copy to output directory
        quality_report = preprocessor.get_quality_report()
        if quality_report:
            click.echo(f"Quality reports:")
            for sample, report_data in quality_report.items():
                if report_data.get('status') == 'completed':
                    report_file = report_data.get('report_file')
                    if report_file and Path(report_file).exists():
                        # Copy quality report to output directory
                        import shutil
                        output_report = Path(output_dir) / f"{sample}.html"
                        shutil.copy2(report_file, output_report)
                        click.echo(f"  {sample}: {output_report}")
                    else:
                        click.echo(f"  {sample}: {report_file}")
        
    except Exception as e:
        logger.error(f"Preprocessing failed: {e}")
        click.echo(f"Error: Preprocessing failed - {e}", err=True)
        sys.exit(1)


if __name__ == '__main__':
    main()
