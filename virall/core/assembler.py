"""
Main viral genome assembler class implementing hybrid assembly strategies.
"""

import os
import subprocess
import tempfile
import shutil
import glob
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, Any
from loguru import logger

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .preprocessor import Preprocessor
from .viral_identifier import ViralIdentifier
from .validator import AssemblyValidator
from .gene_predictor import ViralGenePredictor
from .rnabloom_assembler import RNABloomAssembler


class ViralAssembler:
    """
    Main class for viral genome assembly from short and long reads.
    
    This class implements a comprehensive pipeline that:
    1. Preprocesses input reads (quality control, trimming, error correction)
    2. Identifies viral sequences using machine learning approaches
    3. Performs hybrid assembly combining short and long reads
    4. Polishes and validates the assembled genomes
    """
    
    def __init__(
        self,
        output_dir: Union[str, Path],
        threads: int = 8,
        memory: str = "16G",
        config: Optional[Dict] = None,
        rna_mode: bool = False,
        mem_efficient: bool = False
    ):
        """
        Initialize the viral assembler.
        
        Args:
            output_dir: Directory for output files
            threads: Number of threads to use
            memory: Memory allocation for assembly
            config: Configuration dictionary
            rna_mode: Whether to use RNA-specific parameters
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.threads = threads
        self.memory = memory
        # Merge config with defaults, ensuring database paths are absolute
        default_config = self._get_default_config()
        if config:
            # Merge user config with defaults
            self.config = {**default_config, **config}
            # Ensure database paths are absolute
            if 'databases' in self.config:
                for db_key, db_path in self.config['databases'].items():
                    if isinstance(db_path, str) and not Path(db_path).is_absolute():
                        # Convert relative paths to absolute
                        software_dir = Path(__file__).parent.parent.parent
                        self.config['databases'][db_key] = str(software_dir / db_path)
        else:
            self.config = default_config
        self.rna_mode = rna_mode or self.config.get("rna_mode", False)
        self.mem_efficient = mem_efficient
        
        # Adjust parameters for RNA mode
        if self.rna_mode:
            if "min_contig_length" not in self.config:
                self.config["min_contig_length"] = 500  # Shorter for RNA
            if "viral_confidence_threshold" not in self.config:
                self.config["viral_confidence_threshold"] = 0.7  # Lower threshold for RNA
        
        # Initialize components
        self.preprocessor = Preprocessor(threads=threads)
        self.viral_identifier = ViralIdentifier(config=self.config, threads=threads, memory=memory)
        self.validator = AssemblyValidator(threads=threads, config=self.config)
        self.gene_predictor = ViralGenePredictor(
            threads=threads, 
            config=self.config,
            vog_db_path=self.config.get('databases', {}).get('vog_db_path')
        )
        
        # Disable RNA-Bloom in favor of SPAdes-only workflow
        self.rnabloom_assembler = None
        
        logger.info(f"ViralAssembler initialized with {threads} threads, {memory} memory")
    
    def _detect_single_cell_data(self, reads: Dict[str, str]) -> bool:
        """
        Detect if the data appears to be single-cell RNA-seq based on read ID patterns.
        
        Args:
            reads: Dictionary containing read file paths
            
        Returns:
            True if single-cell data is detected, False otherwise
        """
        # Check for cell barcode patterns in read IDs
        for read_type, read_path in reads.items():
            if read_type in ["short_1", "short_2"] and read_path:
                try:
                    # Handle gzipped files
                    import gzip
                    if str(read_path).endswith('.gz'):
                        f = gzip.open(read_path, 'rt')
                    else:
                        f = open(read_path, 'r')
                    
                    with f:
                        for i, line in enumerate(f):
                            if i >= 10:  # Check first 10 reads
                                break
                            if line.startswith('@'):
                                # Look for cell barcode patterns in read IDs
                                read_id = line.strip()
                                if '_cell_' in read_id and '_umi_' in read_id:
                                    logger.info("Detected single-cell data based on read ID patterns")
                                    return True
                                # Also check for 10X Genomics patterns
                                if 'CB:Z:' in read_id or 'UB:Z:' in read_id:
                                    logger.info("Detected 10X Genomics single-cell data")
                                    return True
                except Exception as e:
                    logger.debug(f"Could not check read file {read_path}: {e}")
                    continue
        
        # No single-cell patterns detected - return False silently
        return False
    
    def _combine_transcript_files(self, long_file: str, short_file: str, output_file: Path) -> None:
        """
        Combine long and short transcript files for comprehensive viral identification.
        
        Args:
            long_file: Path to long transcripts file
            short_file: Path to short transcripts file  
            output_file: Path to combined output file
        """
        logger.info(f"Combining long and short transcripts for comprehensive viral identification")
        
        from Bio import SeqIO
        
        # Read all sequences from both files
        all_sequences = []
        
        # Add long transcripts
        if Path(long_file).exists():
            for record in SeqIO.parse(long_file, "fasta"):
                all_sequences.append(record)
            logger.info(f"Added {len(list(SeqIO.parse(long_file, 'fasta')))} long transcripts")
        
        # Add short transcripts
        if Path(short_file).exists():
            for record in SeqIO.parse(short_file, "fasta"):
                all_sequences.append(record)
            logger.info(f"Added {len(list(SeqIO.parse(short_file, 'fasta')))} short transcripts")
        
        # Write combined file
        SeqIO.write(all_sequences, output_file, "fasta")
        logger.info(f"Combined {len(all_sequences)} transcripts into {output_file}")
    
    def _find_installation_directory(self) -> Path:
        """Find the actual installation directory where databases are located."""
        # Start from the current file location
        current_file = Path(__file__)
        
        # Try different approaches to find the installation directory
        
        # Method 1: Look for a databases directory in the current package structure
        # Go up from virall/core/ to virall/ to parent (installation root)
        package_dir = current_file.parent.parent  # virall/
        installation_dir = package_dir.parent     # parent of virall/
        
        # Check if databases directory exists in installation
        if (installation_dir / "databases").exists():
            return installation_dir
        
        # Method 2: Look for setup.py or other installation markers
        if (installation_dir / "setup.py").exists():
            return installation_dir
        
        # Method 3: Check if we're in a development installation
        # Look for common development markers
        if (installation_dir / "README.md").exists() or (installation_dir / "requirements.txt").exists():
            return installation_dir
        
        # Method 4: Fall back to current working directory if databases exist there
        cwd = Path.cwd()
        if (cwd / "databases").exists():
            return cwd
        
        # Method 5: Fall back to package directory (original behavior)
        logger.warning("Could not find installation directory, falling back to package directory")
        return installation_dir
    
    def _get_default_config(self) -> Dict:
        """Get default configuration parameters."""
        # Get software directory for database paths - look in installation directory
        # Try to find the actual installation directory, not site-packages
        software_dir = self._find_installation_directory()
        
        return {
            "min_contig_length": 1000,
            "min_coverage": 5,
            "viral_confidence_threshold": 0.8,
            "assembly_strategy": "hybrid",  # hybrid, short_only, long_only
            "polishing_iterations": 3,
            "quality_threshold": 20,
            "trim_adapters": True,
            "error_correction": True,
            "databases": {
                "vog_db_path": str(software_dir / "databases" / "vog_db"),
                "kaiju_db_path": str(software_dir / "databases" / "kaiju_db"),
                "checkv_db_path": str(software_dir / "databases" / "checkv_db")
            },
            "viral_identification": {
                "viral_filtering": True,
                "confidence_threshold": 0.8,
                "use_ml_model": True,
                "use_database_search": True,
                "min_viral_length": 1000
            },
            # Deprecated Flye parameters (kept for backward compat; unused)
            "genome_size": "1m",
            "flye_iterations": 1,
            "flye_min_overlap": 1000,
            "flye_read_mode": "raw",
            "flye_read_error": 0.1,
            # Long-read subsampling parameters
            "max_long_reads": 50000,  # Subsample if more than this many reads
            "long_read_subsample_size": 20000  # Subsample to this many reads
        }
    
    def assemble(
        self,
        short_reads_1: Optional[Union[str, Path]] = None,
        short_reads_2: Optional[Union[str, Path]] = None,
        long_reads: Optional[Union[str, Path]] = None,
        single_reads: Optional[Union[str, Path]] = None,
        reference: Optional[Union[str, Path]] = None
    ) -> Dict[str, Union[str, List[str], Dict]]:
        """
        Main assembly pipeline.
        
        Args:
            short_reads_1: Path to first mate of paired-end reads
            short_reads_2: Path to second mate of paired-end reads  
            long_reads: Path to long reads (PacBio/ONT)
            single_reads: Path to single-end reads
            reference: Optional reference genome for guided assembly
            
        Returns:
            Dictionary containing assembly results and statistics
        """
        logger.info("Starting viral genome analysis pipeline")
        
        # Check if reference-guided assembly should be used
        if reference:
            logger.info(f"Reference-guided assembly mode: {reference}")
            result = self._reference_guided_assembly(
                short_reads_1, short_reads_2, long_reads, single_reads, reference
            )
            
            # If reference-guided assembly failed (no reference detected), stop here
            if result.get('status') == 'failed':
                logger.warning("Reference-guided assembly failed - reference genome not detected in sample")
                return result
            
            return result
        
        # Store read paths for quantification
        self.config["short_reads_1"] = short_reads_1
        self.config["short_reads_2"] = short_reads_2
        self.config["single_reads"] = single_reads
        self.config["long_reads"] = long_reads
        
        # Step 1: Preprocessing
        logger.info("Step 1: Preprocessing reads")
        preprocessed_reads = self._preprocess_reads(
            short_reads_1, short_reads_2, long_reads, single_reads
        )
        # Make preprocessed read paths available for downstream steps (e.g., quantification)
        try:
            if isinstance(preprocessed_reads, dict):
                if preprocessed_reads.get("short_1"):
                    self.config["short_reads_1"] = preprocessed_reads.get("short_1")
                if preprocessed_reads.get("short_2"):
                    self.config["short_reads_2"] = preprocessed_reads.get("short_2")
                if preprocessed_reads.get("single"):
                    self.config["single_reads"] = preprocessed_reads.get("single")
                if preprocessed_reads.get("long"):
                    self.config["long_reads"] = preprocessed_reads.get("long")
        except Exception as e:
            logger.debug(f"Could not propagate preprocessed read paths to config: {e}")
        
        # Step 2: Single assembly of all reads
        logger.info("Step 2: Performing initial assembly of all reads")
        assembly_results = self._perform_assembly(preprocessed_reads, reference)
        
        # Step 3: Efficient viral identification from contigs
        logger.info("Step 3: Identifying viral contigs from assembly")
        viral_contig_results = self._identify_viral_contigs_efficiently(assembly_results)
        
        # Step 4: Gene prediction and annotation
        logger.info("Step 4: Predicting and annotating viral genes")
        
        # Collect classification data from all assembly types
        all_classification_data = {}
        viral_contig_info = viral_contig_results.get("viral_contig_info", {})
        for assembly_type, results in viral_contig_info.items():
            if "classification_data" in results:
                # VirSorter2 removed - no longer needed
                all_classification_data.update(results["classification_data"])
            # Add Kaiju classification results
            if "kaiju_classification" in results:
                kaiju_result = results["kaiju_classification"]
                if kaiju_result.get("status") == "completed":
                    all_classification_data[f"kaiju_{assembly_type}"] = {
                        "method": "kaiju",
                        "classification_method": "kaiju_contigs_mode",
                        "status": "completed",
                        "classifications": kaiju_result.get("classifications", {}),
                        "summary_file": kaiju_result.get("summary_file", "")
                    }
        
        gene_prediction_results = self._predict_viral_genes(
            viral_contig_results.get("viral_genomes", []),
            all_classification_data
        )
        
        # Step 5: Validation and quality assessment
        logger.info("Step 5: Validating viral assemblies")
        validation_results = self._validate_assemblies(
            viral_contig_results.get("viral_genomes", []),
            assembly_dir=self.output_dir / "01_assemblies"
        )
        
        # Calculate total viral contigs across all assembly types
        total_viral_contigs = 0
        for assembly_type, info in viral_contig_results.get("viral_contig_info", {}).items():
            if "viral_contig_count" in info:
                total_viral_contigs += info.get("viral_contig_count", 0)
            elif "viral_contigs" in info:
                total_viral_contigs += len(info.get("viral_contigs", []))
        
        # Compile results
        results = {
            "assembly_dir": str(self.output_dir),
            "viral_genomes": viral_contig_results.get("viral_genomes", []),
            "viral_contig_info": viral_contig_results,
            "gene_predictions": gene_prediction_results,
            "validation_results": validation_results,
            "statistics": self._generate_statistics(validation_results),
            "total_viral_contigs": total_viral_contigs
        }
        
        # Clean up temporary files
        self.cleanup_temp_files()
        
        logger.info("Efficient assembly pipeline completed successfully")
        return results
    
    def assemble_and_identify(
        self,
        short_reads_1: Optional[str] = None,
        short_reads_2: Optional[str] = None,
        long_reads: Optional[str] = None,
        single_reads: Optional[str] = None,
        reference: Optional[str] = None
    ) -> Dict[str, Union[str, List[str], Dict]]:
        """
        Assemble reads and identify viral contigs (stops before gene prediction).
        
        Args:
            short_reads_1: Path to first mate of paired-end reads
            short_reads_2: Path to second mate of paired-end reads
            long_reads: Path to long reads (PacBio/ONT)
            single_reads: Path to single-end reads
            reference: Optional reference genome for guided assembly
            
        Returns:
            Dictionary containing assembly results and viral contig information
        """
        logger.info("Starting viral genome assembly and identification pipeline")
        
        # Store read paths for later use
        self.config["short_reads_1"] = short_reads_1
        self.config["short_reads_2"] = short_reads_2
        self.config["single_reads"] = single_reads
        self.config["long_reads"] = long_reads
        
        # Step 1: Preprocessing
        logger.info("Step 1: Preprocessing reads")
        preprocessed_reads = self._preprocess_reads(
            short_reads_1, short_reads_2, long_reads, single_reads
        )
        
        # Step 2: Assembly
        logger.info("Step 2: Performing initial assembly of all reads")
        assembly_results = self._perform_assembly(
            preprocessed_reads, reference
        )
        
        # Step 3: Viral contig identification
        logger.info("Step 3: Identifying viral contigs from assembly")
        viral_contig_results = self._identify_viral_contigs_efficiently(assembly_results)
        
        # Calculate statistics
        viral_genomes = viral_contig_results.get("viral_genomes", [])
        total_contigs = sum(
            viral_contig_results.get("viral_contig_info", {}).get(assembly_type, {}).get("viral_contig_count", 0)
            for assembly_type in ["contigs", "scaffolds"]
        )
        total_length = 0
        for genome_file in viral_genomes:
            if os.path.exists(genome_file):
                stats = self.validator.calculate_genome_statistics(genome_file)
                total_length += stats.get("total_length", 0)
        
        results = {
            "status": "completed",
            "assembly_dir": str(self.output_dir),
            "viral_genomes": viral_genomes,
            "viral_contig_info": viral_contig_results.get("viral_contig_info", {}),
            "statistics": {
                "total_contigs": total_contigs,
                "total_length": total_length,
                "average_contig_length": total_length // max(total_contigs, 1)
            }
        }
        
        # Clean up temporary files
        self.cleanup_temp_files()
        
        logger.info("Assembly and viral identification completed successfully")
        return results
    
    def _preprocess_reads(
        self,
        short_reads_1: Optional[Union[str, Path]],
        short_reads_2: Optional[Union[str, Path]], 
        long_reads: Optional[Union[str, Path]],
        single_reads: Optional[Union[str, Path]]
    ) -> Dict[str, str]:
        """Preprocess all input reads."""
        preprocessed = {}
        
        if short_reads_1 and short_reads_2:
            logger.info("Preprocessing paired-end reads")
            preprocessed["short_1"], preprocessed["short_2"] = self.preprocessor.process_paired_reads(
                short_reads_1, short_reads_2
            )
        
        if long_reads:
            logger.info("Preprocessing long reads")
            preprocessed["long"] = self.preprocessor.process_long_reads(long_reads)
        
        if single_reads:
            logger.info("Preprocessing single-end reads")
            preprocessed["single"] = self.preprocessor.process_single_reads(single_reads)
        
        # Persist FastQC reports (if any) from temp dir to output QC directory
        try:
            qc_dir = self.output_dir / "00_qc" / "fastqc"
            qc_dir.mkdir(parents=True, exist_ok=True)
            # Copy any FastQC outputs (html/zip) from the preprocessor temp dir
            for ext in ("*.html", "*.zip"):
                for report in (self.preprocessor.temp_dir).glob(ext):
                    destination = qc_dir / report.name
                    try:
                        shutil.copy2(report, destination)
                        logger.debug(f"Copied FastQC report: {report} -> {destination}")
                    except Exception as e:
                        logger.warning(f"Failed to copy FastQC report {report}: {e}")
        except Exception as e:
            logger.warning(f"Failed to persist FastQC reports: {e}")
        
        return preprocessed
    
    def _identify_viral_sequences(self, reads: Dict[str, str]) -> Dict[str, str]:
        """Identify and filter viral sequences from reads."""
        viral_reads = {}
        
        for read_type, read_path in reads.items():
            logger.info(f"Identifying viral sequences in {read_type} reads")
            viral_path = self.viral_identifier.identify_viral_reads(
                read_path, 
                confidence_threshold=self.config["viral_confidence_threshold"]
            )
            viral_reads[read_type] = viral_path
        
        return viral_reads
    
    def _identify_viral_contigs_efficiently(
        self, 
        assembly_results: Dict[str, str]
    ) -> Dict[str, Union[str, List[str], Dict]]:
        """
        Efficiently identify viral contigs from assembly results.
        This replaces the old read-based approach with direct contig analysis.
        
        Args:
            assembly_results: Dictionary containing assembly output files
            
        Returns:
            Dictionary with viral contig information and genomes
        """
        viral_genomes = []
        viral_contig_info = {}
        
        for assembly_type, contig_file in assembly_results.items():
            if not os.path.exists(contig_file):
                logger.warning(f"Assembly file not found: {contig_file}")
                continue
            
            # Skip assembly graph files (not suitable for viral identification)
            if assembly_type == "assembly_graph":
                logger.info(f"Skipping {assembly_type} - not suitable for viral identification")
                continue
                
            logger.info(f"Identifying viral contigs in {assembly_type} assembly")
            
            # Use the efficient viral identification method
            viral_contigs_dir = self.output_dir / "02_viral_contigs"
            viral_contigs_dir.mkdir(parents=True, exist_ok=True)
            viral_output_file = viral_contigs_dir / f"viral_{assembly_type}.fasta"
            viral_results = self.viral_identifier.identify_viral_contigs_from_assembly(
                contig_file,
                confidence_threshold=self.config["viral_confidence_threshold"],
                output_file=str(viral_output_file)
            )
            
            # Check for viral contigs - handle both normal and RNA mode results
            viral_contig_count = 0
            if viral_results:
                if "viral_contig_count" in viral_results:
                    viral_contig_count = viral_results.get("viral_contig_count", 0)
                elif "viral_contigs" in viral_results:
                    viral_contig_count = len(viral_results.get("viral_contigs", []))
            
            if viral_results and viral_contig_count > 0:
                viral_genomes.append(str(viral_output_file))
                viral_contig_info[assembly_type] = viral_results
                logger.info(f"Found {viral_contig_count} viral contigs in {assembly_type}")
                
                # Run Kaiju-based contig classification and store summary
                try:
                    classifications_dir = self.output_dir / "03_classifications"
                    classifications_dir.mkdir(parents=True, exist_ok=True)
                    kaiju_dir = classifications_dir / f"kaiju_{assembly_type}"
                    kaiju_results = self.viral_identifier.classify_viral_contigs(
                        contigs_file=str(viral_output_file),
                        output_dir=str(kaiju_dir)
                    )
                    viral_contig_info[assembly_type]["kaiju_classification"] = kaiju_results
                    if kaiju_results.get("status") == "completed":
                        logger.info(
                            f"Kaiju classification completed for {assembly_type}. "
                            f"Summary: {kaiju_results.get('summary_file', 'N/A')}"
                        )
                    else:
                        logger.warning(
                            f"Kaiju classification failed for {assembly_type}: "
                            f"{kaiju_results.get('error', 'Unknown error')}"
                        )
                except Exception as e:
                    logger.warning(f"Failed to run Kaiju classification for {assembly_type}: {e}")

                # Run viral contig quantification
                try:
                    quantification_dir = self.output_dir / "06_quantification"
                    quantification_dir.mkdir(parents=True, exist_ok=True)
                    quant_dir = quantification_dir / assembly_type
                    quant_results = self.viral_identifier.quantify_viral_contigs(
                        contigs_file=str(viral_output_file),
                        reads_1=self.config.get("short_reads_1"),
                        reads_2=self.config.get("short_reads_2"),
                        single_reads=self.config.get("single_reads"),
                        long_reads=self.config.get("long_reads"),
                        output_dir=str(quant_dir),
                        classification_data=viral_results.get("kaiju_classifications", {})
                    )
                    viral_contig_info[assembly_type]["quantification"] = quant_results
                    if quant_results.get("status") == "completed":
                        num_quantified = len(quant_results.get("quantification_results", {}))
                        logger.info(
                            f"Viral quantification completed for {assembly_type}. "
                            f"Quantified {num_quantified} contigs. "
                            f"Abundance file: {quant_dir}/contig_abundance.tsv"
                        )
                    else:
                        logger.warning(
                            f"Viral quantification failed for {assembly_type}: "
                            f"{quant_results.get('error', 'Unknown error')}"
                        )
                except Exception as e:
                    logger.warning(f"Failed to run viral quantification for {assembly_type}: {e}")
                
                # Log classification summary if available
                if "classification_summary" in viral_results:
                    summary = viral_results["classification_summary"]
                    logger.info(f"Viral classification summary for {assembly_type}:")
                    logger.info(f"  Viral groups: {summary.get('viral_groups', {})}")
                    logger.info(f"  Average confidence: {summary.get('average_confidence', 0.0)}")
                    logger.info(f"  Average length: {summary.get('average_length', 0)} bp")
                    logger.info(f"  Total hallmark genes: {summary.get('total_hallmark_genes', 0)}")
                    logger.info(f"  High confidence contigs: {summary.get('high_confidence_contigs', 0)}")
            else:
                logger.warning(f"No viral contigs found in {assembly_type}")
        
        return {
            "viral_genomes": viral_genomes,
            "viral_contig_info": viral_contig_info
        }
    
    def _perform_assembly(
        self, 
        reads: Dict[str, str], 
        reference: Optional[Union[str, Path]]
    ) -> Dict[str, str]:
        """Perform assembly based on available read types."""
        assembly_dir = self.output_dir / "01_assemblies"
        assembly_dir.mkdir(parents=True, exist_ok=True)
        
        strategy = self.config["assembly_strategy"]
        
        # Auto-detect strategy if hybrid is requested but not possible
        if strategy == "hybrid":
            if "short_1" in reads and "long" in reads:
                return self._hybrid_assembly(reads, assembly_dir, reference)
            elif "short_1" in reads or "single" in reads:
                logger.info("Auto-detecting strategy: hybrid not possible, using short_only")
                return self._short_read_assembly(reads, assembly_dir, reference)
            elif "long" in reads:
                logger.info("Auto-detecting strategy: hybrid not possible, using long_only")
                return self._long_read_assembly(reads, assembly_dir, reference)
            else:
                raise ValueError("No suitable reads found for assembly")
        
        # Use explicitly requested strategy
        if strategy == "short_only" and ("short_1" in reads or "single" in reads):
            return self._short_read_assembly(reads, assembly_dir, reference)
        elif strategy == "long_only" and "long" in reads:
            return self._long_read_assembly(reads, assembly_dir, reference)
        else:
            raise ValueError(f"Invalid assembly strategy '{strategy}' for available reads")
    
    def _hybrid_assembly(self, reads: Dict[str, str], output_dir: Path, reference: Optional[Union[str, Path]] = None) -> Dict[str, str]:
        """Perform hybrid assembly using SPAdes only."""
        
        logger.info("Running hybrid assembly with SPAdes")
        
        # Convert memory from "16G" format to just number for SPAdes
        memory_value = self.memory.replace('G', '').replace('g', '')
        
        cmd = [
            "spades.py",
            "--threads", str(self.threads),
            "--memory", memory_value,
            "-o", str(output_dir)
        ]
        
        # Add reference genome for guided assembly if provided
        if reference:
            cmd.extend(["--trusted-contigs", str(reference)])
            logger.info(f"Using reference genome for guided assembly: {reference}")
        
        # Decide SPAdes mode flags
        is_single_cell = False
        try:
            is_single_cell = bool(self.config.get("single_cell_mode")) or self._detect_single_cell_data(reads)
        except Exception:
            is_single_cell = bool(self.config.get("single_cell_mode"))
        if self.rna_mode:
            cmd.append("--rnaviral")
            # Do not use SPAdes --sc (genomic single-cell mode) for scRNA-seq
            logger.info("Using SPAdes rna-viral mode")
        else:
            # Use regular SPAdes for DNA assemblies (not --metaviral)
            # Metaviral mode can be too strict and fail to produce final contigs
            logger.info("Using SPAdes standard mode for DNA assembly")
            # Note: --metaviral removed for better compatibility with single-end reads
        
        # Add read files
        if "short_1" in reads and "short_2" in reads:
            cmd.extend(["-1", reads["short_1"], "-2", reads["short_2"]])
        if "long" in reads:
            long_flag = "--pacbio" if str(self.config.get("long_read_tech", "")).lower() == "pacbio" else "--nanopore"
            cmd.extend([long_flag, reads["long"]])
        if "single" in reads:
            cmd.extend(["-s", reads["single"]])
        
        # Run SPAdes
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"SPAdes failed: {result.stderr}")
        
        # In RNA mode, prefer transcripts; if missing, prefer top-level contigs/scaffolds, then fallback
        if self.rna_mode:
            contigs_path = output_dir / "transcripts.fasta"
            if not contigs_path.exists():
                # Prefer standard SPAdes top-level contigs.fasta if present
                top_contigs = output_dir / "contigs.fasta"
                if top_contigs.exists():
                    logger.warning("transcripts.fasta not found, using top-level contigs.fasta")
                    contigs_path = top_contigs
                else:
                    # Fallback to highest-K final_contigs.fasta
                    k_dirs = [d for d in output_dir.iterdir() if d.is_dir() and d.name.startswith('K')]
                    if k_dirs:
                        k_dirs.sort(key=lambda x: int(x.name[1:]) if x.name[1:].isdigit() else 0)
                        highest_k_dir = k_dirs[-1]
                        fallback_contigs = highest_k_dir / "final_contigs.fasta"
                        if fallback_contigs.exists():
                            logger.warning(f"transcripts.fasta not found, using {fallback_contigs} as fallback")
                            contigs_path = fallback_contigs
            # Choose scaffolds file preference in RNA mode
            scaffolds_path = output_dir / "hard_filtered_transcripts.fasta"
            if not scaffolds_path.exists():
                std_scaffolds = output_dir / "scaffolds.fasta"
                if std_scaffolds.exists():
                    scaffolds_path = std_scaffolds
            return {
                "contigs": str(contigs_path),
                "scaffolds": str(scaffolds_path),
                "assembly_graph": str(output_dir / "assembly_graph.fastg")
            }
        else:
            return {
                "contigs": str(output_dir / "contigs.fasta"),
                "scaffolds": str(output_dir / "scaffolds.fasta"),
                "assembly_graph": str(output_dir / "assembly_graph.fastg")
            }
    
    def _short_read_assembly(self, reads: Dict[str, str], output_dir: Path, reference: Optional[Union[str, Path]] = None) -> Dict[str, str]:
        """Perform short-read only assembly with SPAdes."""
        
        logger.info("Running short-read assembly with SPAdes")
        
        # Convert memory from "16G" format to just number for SPAdes
        memory_value = self.memory.replace('G', '').replace('g', '')
        
        cmd = [
            "spades.py",
            "--threads", str(self.threads),
            "--memory", memory_value,
            "-o", str(output_dir)
        ]
        
        # Add reference genome for guided assembly if provided
        if reference:
            cmd.extend(["--trusted-contigs", str(reference)])
            logger.info(f"Using reference genome for guided assembly: {reference}")
        
        # Decide SPAdes mode flags
        is_single_cell = False
        try:
            is_single_cell = bool(self.config.get("single_cell_mode")) or self._detect_single_cell_data(reads)
        except Exception:
            is_single_cell = bool(self.config.get("single_cell_mode"))
        if self.rna_mode:
            cmd.append("--rnaviral")
            # Do not use SPAdes --sc (genomic single-cell mode) for scRNA-seq
            logger.info("Using SPAdes rna-viral mode")
        else:
            # Use regular SPAdes for DNA assemblies (not --metaviral)
            # Metaviral mode can be too strict and fail to produce final contigs
            logger.info("Using SPAdes standard mode for DNA assembly")
            # Note: --metaviral removed for better compatibility with single-end reads
        
        if "short_1" in reads and "short_2" in reads:
            cmd.extend(["-1", reads["short_1"], "-2", reads["short_2"]])
        if "single" in reads:
            cmd.extend(["-s", reads["single"]])
        
        # Debug: Log the command being run
        logger.info(f"Running SPAdes command: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"SPAdes stderr: {result.stderr}")
            logger.error(f"SPAdes stdout: {result.stdout}")
            raise RuntimeError(f"SPAdes failed: {result.stderr}")
        
        # Log what files SPAdes actually produced (for debugging)
        logger.debug(f"SPAdes completed. Checking output directory: {output_dir}")
        output_files = list(output_dir.glob("*.fasta")) + list(output_dir.glob("*.fa"))
        if output_files:
            logger.debug(f"Found {len(output_files)} FASTA files in output: {[f.name for f in output_files]}")
        else:
            logger.warning(f"No FASTA files found in {output_dir}. Checking subdirectories...")
            for subdir in output_dir.iterdir():
                if subdir.is_dir():
                    sub_files = list(subdir.glob("*.fasta")) + list(subdir.glob("*.fa"))
                    if sub_files:
                        logger.debug(f"Found files in {subdir.name}: {[f.name for f in sub_files]}")
        
        # In RNA mode, prefer transcripts; if missing, prefer top-level contigs/scaffolds, then fallback
        if self.rna_mode:
            contigs_path = output_dir / "transcripts.fasta"
            if not contigs_path.exists():
                top_contigs = output_dir / "contigs.fasta"
                if top_contigs.exists():
                    logger.warning("transcripts.fasta not found, using top-level contigs.fasta")
                    contigs_path = top_contigs
                else:
                    k_dirs = [d for d in output_dir.iterdir() if d.is_dir() and d.name.startswith('K')]
                    if k_dirs:
                        k_dirs.sort(key=lambda x: int(x.name[1:]) if x.name[1:].isdigit() else 0)
                        highest_k_dir = k_dirs[-1]
                        fallback_contigs = highest_k_dir / "final_contigs.fasta"
                        if fallback_contigs.exists():
                            logger.warning(f"transcripts.fasta not found, using {fallback_contigs} as fallback")
                            contigs_path = fallback_contigs
            scaffolds_path = output_dir / "hard_filtered_transcripts.fasta"
            if not scaffolds_path.exists():
                std_scaffolds = output_dir / "scaffolds.fasta"
                if std_scaffolds.exists():
                    scaffolds_path = std_scaffolds
            return {
                "contigs": str(contigs_path),
                "scaffolds": str(scaffolds_path)
            }
        else:
            # DNA mode: check for standard SPAdes output files
            contigs_path = output_dir / "contigs.fasta"
            scaffolds_path = output_dir / "scaffolds.fasta"
            
            # If standard files don't exist, check for fallback locations
            if not contigs_path.exists():
                logger.warning("contigs.fasta not found in output directory, checking for alternative locations...")
                # Check K-mer directories for final_contigs.fasta
                k_dirs = [d for d in output_dir.iterdir() if d.is_dir() and d.name.startswith('K')]
                if k_dirs:
                    k_dirs.sort(key=lambda x: int(x.name[1:]) if x.name[1:].isdigit() else 0)
                    highest_k_dir = k_dirs[-1]
                    fallback_contigs = highest_k_dir / "final_contigs.fasta"
                    if fallback_contigs.exists():
                        logger.warning(f"Using {fallback_contigs} as fallback for contigs")
                        contigs_path = fallback_contigs
                    else:
                        # Check for before_rr.fasta as last resort
                        before_rr = output_dir / "before_rr.fasta"
                        if before_rr.exists():
                            logger.warning(f"Using {before_rr} as fallback for contigs")
                            contigs_path = before_rr
            
            if not scaffolds_path.exists():
                logger.warning("scaffolds.fasta not found in output directory")
                # If scaffolds don't exist, use contigs as scaffolds
                if contigs_path.exists():
                    logger.warning("Using contigs.fasta as scaffolds (scaffolds not available)")
                    scaffolds_path = contigs_path
            
            return {
                "contigs": str(contigs_path) if contigs_path.exists() else str(output_dir / "contigs.fasta"),
                "scaffolds": str(scaffolds_path) if scaffolds_path.exists() else str(output_dir / "scaffolds.fasta")
            }
    
    def _long_read_assembly(self, reads: Dict[str, str], output_dir: Path, reference: Optional[Union[str, Path]] = None) -> Dict[str, str]:
        """Perform long-read only assembly using Flye (ONT or PacBio)."""
        logger.info("Running long-read only assembly with Flye")

        # Determine long-read technology and input file
        long_reads_file = reads.get("long")
        if not long_reads_file:
            raise ValueError("No long reads provided for long-read assembly")

        subsampled = self._subsample_long_reads_if_needed(long_reads_file)
        input_file = subsampled if subsampled else long_reads_file

        tech = str(self.config.get("long_read_tech", "")).lower()
        if tech not in ("nanopore", "pacbio"):
            logger.warning("long_read_tech not set; defaulting to Nanopore flags for Flye")
            tech = "nanopore"

        flye_dir = output_dir / "flye_assembly"
        flye_dir.mkdir(parents=True, exist_ok=True)

        # Choose Flye input flag by technology
        if tech == "pacbio":
            input_flag = "--pacbio-raw"
        else:
            input_flag = "--nano-raw"

        flye_cmd = [
            "flye",
            input_flag, input_file,
            "--out-dir", str(flye_dir),
            "--threads", str(self.threads),
            "--meta"
        ]

        logger.info(f"Flye command: {' '.join(flye_cmd)}")
        result = subprocess.run(flye_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"Flye failed: {result.stderr}")
            # Fallback to simple longest-reads selection if Flye fails
            return self._long_read_assembly_simple_fallback(reads, output_dir)

        contigs_file = flye_dir / "assembly.fasta"
        if not contigs_file.exists():
            logger.warning("Flye produced no contigs; using fallback strategy")
            return self._long_read_assembly_simple_fallback(reads, output_dir)

        return {
            "contigs": str(contigs_file),
            "scaffolds": str(contigs_file)
        }
    
    
    def _long_read_assembly_simple_fallback(self, reads: Dict[str, str], output_dir: Path, reference: Optional[Union[str, Path]] = None) -> Dict[str, str]:
        """Simple fallback: select longest reads as contigs when all assemblers fail."""
        logger.info("Running simple read selection strategy (last resort)")
        
        from Bio import SeqIO
        import gzip
        
        # Read all sequences and sort by length
        sequences = []
        if self._is_gzipped(str(reads["long"])):
            fh = gzip.open(reads["long"], "rt")
        else:
            fh = open(reads["long"], "rt")
        with fh as handle:
            for record in SeqIO.parse(handle, "fastq"):
                sequences.append(record)
        
        # Sort by length (descending)
        sequences.sort(key=lambda x: len(x.seq), reverse=True)
        
        # Select top 100 longest sequences as "contigs"
        selected_sequences = sequences[:100]
        
        # Write contigs
        contigs_file = output_dir / "contigs.fasta"
        scaffolds_file = output_dir / "scaffolds.fasta"
        
        with open(contigs_file, "w") as handle:
            for i, record in enumerate(selected_sequences):
                record.id = f"contig_{i+1}_length_{len(record.seq)}"
                record.description = f"Selected read {i+1} (length: {len(record.seq)})"
                SeqIO.write(record, handle, "fasta")
        
        # Copy contigs to scaffolds (same file for simplicity)
        import shutil
        shutil.copy2(contigs_file, scaffolds_file)
        
        logger.info(f"Simple fallback completed: selected {len(selected_sequences)} longest reads as contigs")
        return {
            "contigs": str(contigs_file),
            "scaffolds": str(scaffolds_file)
        }
    
    def _subsample_long_reads_if_needed(self, reads_file: str) -> Optional[str]:
        """Subsample long reads if dataset is too large for memory constraints."""
        from Bio import SeqIO
        import gzip
        import tempfile
        
        # Count total reads first
        total_reads = 0
        if self._is_gzipped(reads_file):
            fh = gzip.open(reads_file, "rt")
        else:
            fh = open(reads_file, "rt")
        with fh as handle:
            for _ in SeqIO.parse(handle, "fastq"):
                total_reads += 1
        
        # Check if subsampling is needed based on memory-efficient mode
        if self.mem_efficient:
            max_reads = self.config.get("max_long_reads", 50000)
            subsample_size = self.config.get("long_read_subsample_size", 20000)
            
            if total_reads > max_reads:
                logger.info(f"Large dataset detected ({total_reads} reads). Subsampling to {subsample_size} longest reads for memory efficiency...")
                
                # Read all sequences and sort by length
                sequences = []
                if self._is_gzipped(reads_file):
                    fh2 = gzip.open(reads_file, "rt")
                else:
                    fh2 = open(reads_file, "rt")
                with fh2 as handle:
                    for record in SeqIO.parse(handle, "fastq"):
                        sequences.append(record)
                
                # Sort by length (descending) and take top N
                sequences.sort(key=lambda x: len(x.seq), reverse=True)
                selected_sequences = sequences[:subsample_size]
                
                # Write subsampled reads to temporary file
                temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fastq.gz', delete=False)
                with gzip.open(temp_file.name, "wt") as handle:
                    for record in selected_sequences:
                        SeqIO.write(record, handle, "fastq")
                
                logger.info(f"Subsampled to {len(selected_sequences)} reads (saved to {temp_file.name})")
                return temp_file.name
        
        return None

    def _is_gzipped(self, file_path: str) -> bool:
        """Return True if the file is gzipped based on magic bytes."""
        try:
            with open(file_path, "rb") as fh:
                return fh.read(2) == b"\x1f\x8b"
        except Exception:
            return file_path.endswith(".gz")
    
    def _extract_viral_genomes(self, assembly_results: Dict[str, str]) -> List[str]:
        """Extract and polish viral genomes from assembly results."""
        viral_genomes = []
        
        for assembly_type, contig_file in assembly_results.items():
            if not os.path.exists(contig_file):
                continue
                
            logger.info(f"Extracting viral genomes from {assembly_type}")
            
            # Load contigs
            contigs = list(SeqIO.parse(contig_file, "fasta"))
            
            # Filter by length and identify viral sequences
            viral_contigs = []
            for contig in contigs:
                if len(contig.seq) >= self.config["min_contig_length"]:
                    # Use viral identifier to check if contig is viral
                    is_viral = self.viral_identifier.is_viral_sequence(
                        str(contig.seq),
                        confidence_threshold=self.config["viral_confidence_threshold"]
                    )
                    
                    if is_viral:
                        viral_contigs.append(contig)
            
            # Save viral genomes
            if viral_contigs:
                viral_file = self.output_dir / f"viral_genomes_{assembly_type}.fasta"
                SeqIO.write(viral_contigs, viral_file, "fasta")
                viral_genomes.append(str(viral_file))
                
                logger.info(f"Found {len(viral_contigs)} viral contigs in {assembly_type}")
        
        return viral_genomes
    
    def _predict_viral_genes(
        self,
        viral_genomes: List[str],
        classification_data: Dict[str, Dict]
    ) -> Dict[str, Dict]:
        """Predict and annotate genes in viral genomes."""
        if not viral_genomes:
            logger.info("No viral genomes found for gene prediction")
            return {"status": "skipped", "reason": "no_viral_genomes"}
        
        logger.info(f"Predicting genes in {len(viral_genomes)} viral genomes")
        
        # Flatten classification data to get contig-level classifications
        flattened_classifications = {}
        for key, data in classification_data.items():
            if isinstance(data, dict) and "classifications" in data:
                # Extract contig-level classifications from nested structure
                classifications = data.get("classifications", {})
                flattened_classifications.update(classifications)
        
        # If no flattened classifications found, try to extract from the data structure
        if not flattened_classifications:
            for key, data in classification_data.items():
                if isinstance(data, dict) and "classifications" in data:
                    classifications = data.get("classifications", {})
                    if isinstance(classifications, dict):
                        flattened_classifications.update(classifications)
        
        logger.info(f"Using {len(flattened_classifications)} contig-level classifications for gene prediction")
        
        # Combine all viral contigs for gene prediction
        viral_contigs_dir = self.output_dir / "02_viral_contigs"
        viral_contigs_dir.mkdir(parents=True, exist_ok=True)
        combined_contigs_file = viral_contigs_dir / "all_viral_contigs.fasta"
        self._combine_viral_contigs(viral_genomes, combined_contigs_file)
        
        # Run comprehensive gene prediction
        gene_predictions_dir = self.output_dir / "05_gene_predictions"
        gene_predictions_dir.mkdir(parents=True, exist_ok=True)
        gene_prediction_results = self.gene_predictor.predict_genes_comprehensive(
            str(combined_contigs_file),
            flattened_classifications,
            str(gene_predictions_dir)
        )
        
        # Clean up temporary files
        self.gene_predictor.cleanup()
        
        return gene_prediction_results
    
    def _combine_viral_contigs(self, viral_genomes: List[str], output_file: Path) -> None:
        """Combine all viral contigs into a single file."""
        all_contigs = []
        
        for genome_file in viral_genomes:
            for record in SeqIO.parse(genome_file, "fasta"):
                all_contigs.append(record)
        
        SeqIO.write(all_contigs, output_file, "fasta")
        logger.info(f"Combined {len(all_contigs)} viral contigs into {output_file}")
    
    def _validate_assemblies(self, viral_genomes: List[str], assembly_dir: Optional[Path] = None) -> Dict[str, Dict]:
        """Validate and assess quality of assembled viral genomes."""
        validation_results = {}
        
        # Assembly quality evaluation removed
        
        # Individual viral genome validation
        for genome_file in viral_genomes:
            logger.info(f"Validating {genome_file}")
            
            # Calculate basic statistics (useful for viral contigs)
            stats = self.validator.calculate_genome_statistics(genome_file)
            
            # Run CheckV for viral contig quality assessment
            # Create a specific output directory for CheckV to avoid temp directory issues
            quality_assessment_dir = self.output_dir / "04_quality_assessment"
            quality_assessment_dir.mkdir(parents=True, exist_ok=True)
            checkv_output_dir = quality_assessment_dir / "checkv_results" / Path(genome_file).stem
            checkv_results = self.validator.run_checkv(genome_file, output_dir=str(checkv_output_dir))
            
            validation_results[genome_file] = {
                "checkv": checkv_results,
                "statistics": stats
            }
        
        return validation_results
    
    def _generate_statistics(self, validation_results: Dict) -> Dict:
        """Generate overall assembly statistics."""
        if not validation_results:
            return {
                "total_contigs": 0,
                "total_length": 0,
                "average_length": 0,
                "average_contig_length": 0
            }
        
        total_genomes = len(validation_results)
        total_contigs = sum(
            result.get("statistics", {}).get("num_contigs", 0)
            for result in validation_results.values()
        )
        total_length = sum(
            result.get("statistics", {}).get("total_length", 0)
            for result in validation_results.values()
        )
        
        return {
            "total_viral_genomes": total_genomes,
            "total_contigs": total_contigs,
            "total_length": total_length,
            "average_length": total_length / total_contigs if total_contigs > 0 else 0,
            "average_contig_length": total_length / total_contigs if total_contigs > 0 else 0
        }
    
    def _reference_guided_assembly(
        self,
        short_reads_1: Optional[Union[str, Path]] = None,
        short_reads_2: Optional[Union[str, Path]] = None,
        long_reads: Optional[Union[str, Path]] = None,
        single_reads: Optional[Union[str, Path]] = None,
        reference: Optional[Union[str, Path]] = None
    ) -> Dict[str, Union[str, List[str], Dict]]:
        """
        Perform reference-guided assembly using SPAdes.
        
        Args:
            short_reads_1: Path to first mate of paired-end reads
            short_reads_2: Path to second mate of paired-end reads
            long_reads: Path to long reads (PacBio/ONT)
            single_reads: Path to single-end reads
            reference: Reference genome for guided assembly
            
        Returns:
            Dictionary containing reference-guided assembly results
        """
        logger.info("Starting reference-guided assembly")
        
        # Validate reference file
        if not os.path.exists(str(reference)):
            raise FileNotFoundError(f"Reference genome not found: {reference}")
        
        # Step 1: Preprocess reads
        logger.info("Step 1: Preprocessing reads for reference-guided assembly")
        preprocessed_reads = self._preprocess_reads(
            short_reads_1, short_reads_2, long_reads, single_reads
        )
        
        # Step 2: Choose appropriate assembly strategy based on input data
        has_short_reads = "short_1" in preprocessed_reads or "single" in preprocessed_reads
        has_long_reads = "long" in preprocessed_reads
        
        if has_short_reads and has_long_reads:
            # Both short and long reads - build long-read consensus against reference and polish with short reads
            logger.info("Step 2: Running reference-guided long-read assembly + short-read polishing (minimap2 + bcftools + pilon)")
            guided_assembly_results = self._run_reference_guided_long_then_polish(
                preprocessed_reads, reference
            )
        elif has_short_reads and not has_long_reads:
            # Only short reads - use SPAdes reference-guided assembly
            logger.info("Step 2: Running reference-guided assembly with SPAdes (short reads only)")
            guided_assembly_results = self._run_reference_guided_spades(
                preprocessed_reads, reference
            )
        elif has_long_reads and not has_short_reads:
            # Only long reads - build consensus against reference (no short-read polishing)
            logger.info("Step 2: Running reference-guided long-read consensus (minimap2 + bcftools)")
            guided_assembly_results = self._run_reference_guided_long_only_consensus(
                preprocessed_reads, reference
            )
        else:
            logger.warning("No suitable reads found for reference-guided assembly")
            guided_assembly_results = None
        
        if not guided_assembly_results:
            logger.warning("Reference-guided assembly failed - no reference genome detected")
            return {
                "status": "failed",
                "error": "Reference-guided assembly failed - no reference genome detected",
                "assembly_dir": str(self.output_dir),
                "viral_genomes": [],
                "total_viral_contigs": 0,
                "statistics": {
                    "total_contigs": 0,
                    "total_length": 0,
                    "average_contig_length": 0
                },
                "reference_guided_failed": True
            }
        
        # Step 3: Evaluate reference-guided assembly
        logger.info("Step 3: Evaluating reference-guided assembly")
        evaluation_results = self._evaluate_reference_assembly(
            guided_assembly_results, reference
        )
        
        # Step 4: Run complete pipeline on filtered contigs
        logger.info("Step 4: Running complete pipeline on reference-matching contigs")
        pipeline_results = self._run_complete_pipeline_on_filtered_contigs(
            guided_assembly_results, preprocessed_reads
        )
        
        # Merge evaluation and pipeline results
        if pipeline_results:
            evaluation_results.update(pipeline_results)
        
        # Clean up temporary files
        self.cleanup_temp_files()
        
        logger.info("Reference-guided assembly with complete pipeline completed successfully")
        return evaluation_results
    
    def _run_complete_pipeline_on_filtered_contigs(
        self, 
        guided_assembly_results: Dict[str, str], 
        preprocessed_reads: Dict[str, str]
    ) -> Dict[str, Any]:
        """
        Run the complete viral analysis pipeline on reference-matching contigs.
        
        Args:
            guided_assembly_results: Results from reference-guided assembly
            preprocessed_reads: Preprocessed read files
            
        Returns:
            Dictionary with complete pipeline results
        """
        try:
            # Get the filtered contigs file
            filtered_contigs_file = guided_assembly_results.get("contigs")
            if not filtered_contigs_file or not os.path.exists(filtered_contigs_file):
                logger.warning("No valid filtered contigs file found for pipeline")
                return {}
            
            logger.info(f"Running complete pipeline on {filtered_contigs_file}")
            
            # Step 1: Viral classification with Kaiju
            logger.info("Running viral classification on filtered contigs")
            classifications_dir = self.output_dir / "03_classifications"
            classifications_dir.mkdir(parents=True, exist_ok=True)
            kaiju_dir = classifications_dir / "kaiju_reference_guided"
            
            viral_results = self.viral_identifier.identify_viral_contigs_from_assembly(
                filtered_contigs_file
            )
            
            # Step 2: Gene prediction and annotation
            logger.info("Running gene prediction and annotation on filtered contigs")
            gene_predictions_dir = self.output_dir / "05_gene_predictions"
            gene_predictions_dir.mkdir(parents=True, exist_ok=True)
            
            # Get classification data for gene prediction
            classification_data = viral_results.get("kaiju_classifications", {})
            
            gene_prediction_results = self.gene_predictor.predict_genes_comprehensive(
                filtered_contigs_file,
                classification_data,
                str(gene_predictions_dir)
            )
            
            # Step 3: Quantification
            logger.info("Running quantification on filtered contigs")
            quantification_dir = self.output_dir / "06_quantification"
            quantification_dir.mkdir(parents=True, exist_ok=True)
            quant_dir = quantification_dir / "reference_guided"
            
            quant_results = self.viral_identifier.quantify_viral_contigs(
                contigs_file=filtered_contigs_file,
                reads_1=preprocessed_reads.get("short_1"),
                reads_2=preprocessed_reads.get("short_2"),
                single_reads=preprocessed_reads.get("single"),
                long_reads=preprocessed_reads.get("long"),
                output_dir=str(quant_dir),
                classification_data=classification_data
            )
            
            # Step 4: Organize viral contigs
            logger.info("Organizing viral contigs")
            viral_contigs_dir = self.output_dir / "02_viral_contigs"
            viral_contigs_dir.mkdir(parents=True, exist_ok=True)
            
            # Copy filtered contigs to viral contigs directory
            import shutil
            viral_contigs_file = viral_contigs_dir / "reference_matching_viral_contigs.fasta"
            shutil.copy2(filtered_contigs_file, viral_contigs_file)
            
            # Return complete pipeline results
            return {
                "viral_classifications": viral_results,
                "gene_predictions": gene_prediction_results,
                "quantification": quant_results,
                "viral_contigs_file": str(viral_contigs_file),
                "pipeline_status": "completed"
            }
            
        except Exception as e:
            logger.error(f"Error running complete pipeline on filtered contigs: {e}")
            return {
                "pipeline_status": "failed",
                "error": str(e)
            }
    
    def _run_reference_guided_flye(
        self, 
        reads: Dict[str, str], 
        reference: Union[str, Path]
    ) -> Optional[Dict[str, str]]:
        """
        Run Flye for de novo assembly with long reads, then filter against reference.
        
        Args:
            reads: Dictionary containing preprocessed read files
            reference: Reference genome file
            
        Returns:
            Dictionary containing assembly results or None if failed
        """
        output_dir = self.output_dir / "reference_guided_assembly"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Running Flye de novo assembly with long reads, then filtering against reference: {reference}")
        
        # Check if we need to subsample for memory constraints
        subsampled_reads = self._subsample_long_reads_if_needed(reads["long"])
        input_file = subsampled_reads if subsampled_reads else reads["long"]
        
        # Step 1: Run Flye de novo assembly
        flye_output_dir = output_dir / "flye_assembly"
        flye_cmd = [
            "flye",
            "--nano-raw", input_file,
            "--out-dir", str(flye_output_dir),
            "--threads", str(self.threads)
        ]
        
        # Always use metagenome mode for viral assembly
        flye_cmd.extend(["--meta"])
        if self.rna_mode:
            logger.info("Using Flye metagenome mode for RNA viral assembly")
        else:
            logger.info("Using Flye metagenome mode for viral assembly")
        
        logger.info(f"Flye command: {' '.join(flye_cmd)}")
        
        result = subprocess.run(flye_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.warning(f"Flye assembly failed: {result.stderr}")
            return None
        
        # Check if assembly produced output
        contigs_file = flye_output_dir / "assembly.fasta"
        if not contigs_file.exists():
            logger.warning("Flye assembly produced no contigs")
            return None
        
        # Step 2: Filter contigs against reference using minimap2
        logger.info("Filtering assembly contigs against reference genome using minimap2")
        filtered_contigs_file = output_dir / "filtered_contigs.fasta"
        
        # Run minimap2 to find contigs similar to reference
        minimap2_results = self._filter_contigs_against_reference_minimap2(
            contigs_file, reference, filtered_contigs_file
        )
        
        if not minimap2_results or not minimap2_results.get('has_hits', False):
            logger.warning("No contigs found similar to reference genome")
            return None
        
        # Copy the filtered contigs to the main output
        final_contigs_file = output_dir / "contigs.fasta"
        shutil.copy2(filtered_contigs_file, final_contigs_file)
        
        logger.info(f"Reference-guided assembly completed: {minimap2_results.get('num_hits', 0)} contigs similar to reference")
        
        return {
            "contigs": str(final_contigs_file),
            "scaffolds": str(final_contigs_file),  # Flye doesn't separate contigs/scaffolds
            "reference_hits": minimap2_results
        }
    
    def _filter_contigs_against_reference_minimap2(
        self, 
        contigs_file: Path, 
        reference: Union[str, Path], 
        output_file: Path
    ) -> Optional[Dict]:
        """
        Filter contigs against reference using minimap2.
        
        Args:
            contigs_file: Path to contigs file
            reference: Path to reference genome
            output_file: Path to save filtered contigs
            
        Returns:
            Dictionary with minimap2 results or None if failed
        """
        try:
            # Run minimap2 alignment
            minimap2_output = output_file.parent / "minimap2_results.paf"
            minimap2_cmd = [
                "minimap2",
                "-c",  # Output CIGAR in PAF
                "-x", "asm5",  # Use asm5 preset for assembly-to-assembly alignment
                str(reference),
                str(contigs_file)
            ]
            
            logger.info(f"Running minimap2: {' '.join(minimap2_cmd)}")
            result = subprocess.run(minimap2_cmd, capture_output=True, text=True)
            if result.returncode != 0:
                logger.warning(f"minimap2 alignment failed: {result.stderr}")
                return None
            
            # Write minimap2 output to file
            with open(minimap2_output, 'w') as f:
                f.write(result.stdout)
            
            # Parse minimap2 results and filter contigs
            hits = set()
            all_alignments = []
            if minimap2_output.exists():
                with open(minimap2_output, 'r') as f:
                    for line in f:
                        if line.strip():
                            parts = line.strip().split('\t')
                            if len(parts) >= 12:
                                contig_id = parts[0]
                                contig_length = int(parts[1])
                                contig_start = int(parts[2])
                                contig_end = int(parts[3])
                                identity = float(parts[9]) / float(parts[10]) if int(parts[10]) > 0 else 0
                                alignment_length = int(parts[10])
                                
                                # Store alignment info for debugging
                                all_alignments.append({
                                    'contig_id': contig_id,
                                    'contig_length': contig_length,
                                    'identity': identity,
                                    'alignment_length': alignment_length,
                                    'coverage': alignment_length / contig_length if contig_length > 0 else 0
                                })
                                
                                # Filter criteria: significant alignment with good identity
                                alignment_coverage = alignment_length / contig_length if contig_length > 0 else 0
                                if identity >= 0.6 and alignment_coverage >= 0.3:  # Relaxed: 60% identity, 30% coverage
                                    hits.add(contig_id)
            
            # Log alignment statistics for debugging
            if all_alignments:
                logger.info(f"Found {len(all_alignments)} minimap2 alignments")
                if all_alignments:
                    best_identity = max(a['identity'] for a in all_alignments)
                    best_coverage = max(a['coverage'] for a in all_alignments)
                    logger.info(f"Best alignment: {best_identity:.2%} identity, {best_coverage:.2%} coverage")
            
            if not hits:
                logger.warning("No significant minimap2 hits found with relaxed criteria (≥60% identity, ≥30% coverage)")
                logger.info("Consider checking if the reference genome is appropriate for this sample")
                return None
            
            # Filter contigs to keep only those with hits
            filtered_contigs = []
            with open(contigs_file, 'r') as f:
                current_contig = None
                current_seq = []
                
                for line in f:
                    if line.startswith('>'):
                        # Save previous contig if it had hits
                        if current_contig and current_contig in hits:
                            filtered_contigs.append(f">{current_contig}\n")
                            filtered_contigs.append(''.join(current_seq) + '\n')
                        
                        # Start new contig
                        current_contig = line.strip()[1:].split()[0]  # Get contig ID
                        current_seq = []
                    else:
                        current_seq.append(line.strip())
                
                # Don't forget the last contig
                if current_contig and current_contig in hits:
                    filtered_contigs.append(f">{current_contig}\n")
                    filtered_contigs.append(''.join(current_seq) + '\n')
            
            # Write filtered contigs
            with open(output_file, 'w') as f:
                f.writelines(filtered_contigs)
            
            logger.info(f"Filtered {len(hits)} contigs with significant hits to reference using minimap2")
            
            return {
                'has_hits': True,
                'num_hits': len(hits),
                'hits': list(hits),
                'minimap2_file': str(minimap2_output)
            }
            
        except Exception as e:
            logger.warning(f"Error filtering contigs against reference with minimap2: {e}")
            return None
    
    def _run_reference_guided_spades(
        self, 
        reads: Dict[str, str], 
        reference: Union[str, Path]
    ) -> Optional[Dict[str, str]]:
        """
        Run SPAdes with reference-guided assembly.
        
        Args:
            reads: Dictionary containing preprocessed read files
            reference: Reference genome file
            
        Returns:
            Dictionary containing assembly results or None if failed
        """
        output_dir = self.output_dir / "reference_guided_assembly"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Running SPAdes reference-guided assembly with {reference}")
        
        # Build SPAdes command for reference-guided assembly
        cmd = [
            "spades.py",
            "--threads", str(self.threads),
            # SPAdes expects integer GB for memory; convert if value like "16G"
            "--memory", str(int(str(self.memory).rstrip('Gg')) if str(self.memory).rstrip('Gg').isdigit() else 16),
            "-o", str(output_dir)
        ]
        
        # Add read files - SPAdes requires at least one short read library
        has_short_reads = ("short_1" in reads and "short_2" in reads) or "single" in reads
        
        if not has_short_reads:
            logger.error("SPAdes requires at least one short read library (paired-end or single-end)")
            return None
            
        if "short_1" in reads and "short_2" in reads:
            cmd.extend(["-1", reads["short_1"], "-2", reads["short_2"]])
        if "single" in reads:
            cmd.extend(["-s", reads["single"]])
        if "long" in reads:
            long_flag = "--pacbio" if str(self.config.get("long_read_tech", "")).lower() == "pacbio" else "--nanopore"
            cmd.extend([long_flag, reads["long"]])
        
        # Add reference genome for guided assembly
        cmd.extend(["--trusted-contigs", str(reference)])
        
        # Add mode flags (rnaviral/sc/standard)
        is_single_cell = bool(self.config.get("single_cell_mode", False))
        if self.rna_mode:
            cmd.append("--rnaviral")
            # Do not use SPAdes --sc (genomic single-cell mode) for scRNA-seq
            logger.info("Using SPAdes rna-viral mode")
        else:
            # Use regular SPAdes for DNA assemblies (not --metaviral)
            # Metaviral mode can be too strict and fail to produce final contigs
            logger.info("Using SPAdes standard mode for DNA assembly")
            # Note: --metaviral removed for better compatibility
        
        logger.info(f"SPAdes command: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.warning(f"SPAdes reference-guided assembly failed: {result.stderr}")
            return None
        
        # Check if assembly produced output (handle RNA mode)
        if self.rna_mode:
            contigs_file = output_dir / "transcripts.fasta"
            scaffolds_file = output_dir / "hard_filtered_transcripts.fasta"
        else:
            contigs_file = output_dir / "contigs.fasta"
            scaffolds_file = output_dir / "scaffolds.fasta"
        
        if not contigs_file.exists():
            logger.warning("SPAdes reference-guided assembly produced no contigs/transcripts")
            return None
        
        # Check if the assembly contains significant hits to the reference
        reference_hits = self._check_reference_hits(contigs_file, reference)
        if not reference_hits or not reference_hits.get('has_hits', False):
            logger.warning("Minimap2 alignment failed or no significant hits to reference genome found in assembly")
            logger.info("Proceeding with all assembled contigs/transcripts (reference filtering disabled)")
            # Create a dummy reference_hits object for compatibility
            reference_hits = {
                'has_hits': True,
                'num_hits': 0,
                'matching_contigs': [],
                'alignment_results': []
            }
        
        # Filter contigs to only include those matching the reference
        filtered_contigs_file = output_dir / "reference_matching_contigs.fasta"
        
        # If we have specific matching contigs from Minimap2, filter them
        if reference_hits['matching_contigs']:
            if not self._filter_contigs_by_reference_hits(
                contigs_file, 
                reference_hits['matching_contigs'], 
                filtered_contigs_file
            ):
                logger.warning("Failed to filter contigs by reference hits")
                return None
        else:
            # If Minimap2 failed or no specific matches, use all contigs
            logger.info("Using all assembled contigs/transcripts (no reference filtering)")
            import shutil
            shutil.copy2(contigs_file, filtered_contigs_file)
        
        logger.info("Reference-guided assembly completed successfully")
        return {
            "contigs": str(filtered_contigs_file),
            "scaffolds": str(filtered_contigs_file),  # Use filtered contigs for both
            "original_contigs": str(contigs_file),
            "original_scaffolds": str(scaffolds_file) if scaffolds_file.exists() else str(contigs_file),
            "reference_hits": reference_hits
        }

    def _run_reference_guided_long_only_consensus(
        self,
        reads: Dict[str, str],
        reference: Union[str, Path]
    ) -> Optional[Dict[str, str]]:
        """
        Reference-guided consensus using only long reads (no short-read polishing):
        minimap2 + samtools + bcftools consensus.
        Uses map-ont for Nanopore and map-hifi for PacBio.
        """
        try:
            output_dir = self.output_dir / "reference_guided_assembly"
            output_dir.mkdir(parents=True, exist_ok=True)

            long_reads = reads.get("long")
            if not long_reads:
                logger.error("Long reads are required for reference-guided consensus")
                return None

            aligned_bam = output_dir / "aligned.bam"
            tech = str(self.config.get("long_read_tech", "")).lower()
            preset = "map-hifi" if tech == "pacbio" else "map-ont"

            # Align and sort
            minimap2_cmd = [
                "minimap2", "-t", str(self.threads), "-ax", preset, str(reference), str(long_reads)
            ]
            sort_cmd = ["samtools", "sort", "-o", str(aligned_bam)]

            logger.info(f"Running minimap2: {' '.join(minimap2_cmd)} | {' '.join(sort_cmd)}")
            p1 = subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE, text=False)
            p2 = subprocess.run(sort_cmd, stdin=p1.stdout, capture_output=True)
            p1.stdout.close()  # type: ignore
            if p1.wait() != 0 or p2.returncode != 0:
                logger.error(f"Alignment failed: {p2.stderr.decode('utf-8', 'ignore')}")
                return None

            # Index BAM
            result = subprocess.run(["samtools", "index", str(aligned_bam)], capture_output=True, text=True)
            if result.returncode != 0:
                logger.error(f"samtools index failed: {result.stderr}")
                return None

            # Variant calling and consensus
            vcf_gz = output_dir / "calls.vcf.gz"
            mpileup_cmd = ["bcftools", "mpileup", "-Ou", "-f", str(reference), str(aligned_bam)]
            call_cmd = ["bcftools", "call", "-mv", "-Oz", "-o", str(vcf_gz)]
            logger.info(f"Calling variants: {' '.join(mpileup_cmd)} | {' '.join(call_cmd)}")
            p3 = subprocess.Popen(mpileup_cmd, stdout=subprocess.PIPE, text=False)
            p4 = subprocess.run(call_cmd, stdin=p3.stdout, capture_output=True)
            p3.stdout.close()  # type: ignore
            if p3.wait() != 0 or p4.returncode != 0:
                logger.error(f"bcftools call failed: {p4.stderr.decode('utf-8', 'ignore')}")
                return None

            # Index VCF and build consensus
            subprocess.run(["bcftools", "index", str(vcf_gz)], check=False)
            consensus_fa = output_dir / "consensus.fa"
            consensus_cmd = ["bcftools", "consensus", "-f", str(reference), str(vcf_gz)]
            logger.info(f"Building consensus: {' '.join(consensus_cmd)} > {consensus_fa}")
            with open(consensus_fa, "w") as fh:
                result = subprocess.run(consensus_cmd, stdout=fh, capture_output=True, text=True)
            if result.returncode != 0 or not consensus_fa.exists():
                logger.error(f"bcftools consensus failed: {result.stderr}")
                return None

            return {
                "contigs": str(consensus_fa),
                "scaffolds": str(consensus_fa),
                "long_alignment_bam": str(aligned_bam),
                "variants_vcf": str(vcf_gz)
            }
        except Exception as e:
            logger.error(f"Reference-guided long-read consensus failed: {e}")
            return None

    def _run_reference_guided_long_then_polish(
        self,
        reads: Dict[str, str],
        reference: Union[str, Path]
    ) -> Optional[Dict[str, str]]:
        """
        Build a draft assembly by aligning long reads to the reference (minimap2 + bcftools)
        and polish with short reads using Pilon.
        """
        try:
            output_dir = self.output_dir / "reference_guided_assembly"
            output_dir.mkdir(parents=True, exist_ok=True)

            long_reads = reads.get("long")
            if not long_reads:
                logger.error("Long reads are required for reference-guided long-read assembly")
                return None

            # 1) Long-read consensus against reference
            long_bam = output_dir / "long_aligned.bam"
            long_bam_index = output_dir / "long_aligned.bam.bai"
            tech = str(self.config.get("long_read_tech", "")).lower()
            preset = "map-pb" if tech == "pacbio" else "map-ont"

            minimap2_cmd = [
                "minimap2", "-t", str(self.threads), "-ax", preset, str(reference), str(long_reads)
            ]
            samtools_sort_cmd = ["samtools", "sort", "-o", str(long_bam)]

            logger.info(f"Running minimap2 for long reads: {' '.join(minimap2_cmd)} | {' '.join(samtools_sort_cmd)}")
            # Pipe minimap2 -> samtools sort
            minimap_proc = subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE, text=False)
            sort_proc = subprocess.run(samtools_sort_cmd, stdin=minimap_proc.stdout, capture_output=True)
            minimap_proc.stdout.close()  # type: ignore
            ret1 = minimap_proc.wait()
            if ret1 != 0 or sort_proc.returncode != 0:
                logger.error(f"Long-read alignment failed: {sort_proc.stderr.decode('utf-8', 'ignore')}")
                return None

            # Index BAM
            result = subprocess.run(["samtools", "index", str(long_bam)], capture_output=True, text=True)
            if result.returncode != 0:
                logger.error(f"samtools index failed: {result.stderr}")
                return None

            # Call variants and build consensus
            variants_bcf = output_dir / "variants.bcf"
            variants_norm_bcf = output_dir / "variants.norm.bcf"
            draft_fa = output_dir / "draft_assembly.fa"

            mpileup_cmd = ["bcftools", "mpileup", "-Ou", "-f", str(reference), str(long_bam)]
            call_cmd = ["bcftools", "call", "-mv", "-Ob", "-o", str(variants_bcf)]
            logger.info(f"Calling variants: {' '.join(mpileup_cmd)} | {' '.join(call_cmd)}")
            mpileup_proc = subprocess.Popen(mpileup_cmd, stdout=subprocess.PIPE, text=False)
            call_proc = subprocess.run(call_cmd, stdin=mpileup_proc.stdout, capture_output=True)
            mpileup_proc.stdout.close()  # type: ignore
            ret2 = mpileup_proc.wait()
            if ret2 != 0 or call_proc.returncode != 0:
                logger.error(f"bcftools call failed: {call_proc.stderr.decode('utf-8', 'ignore')}")
                return None

            norm_cmd = ["bcftools", "norm", "-f", str(reference), "-Ob", "-o", str(variants_norm_bcf), str(variants_bcf)]
            logger.info(f"Normalizing variants: {' '.join(norm_cmd)}")
            result = subprocess.run(norm_cmd, capture_output=True, text=True)
            if result.returncode != 0:
                logger.warning(f"bcftools norm failed, proceeding with raw variants: {result.stderr}")
                variants_for_consensus = variants_bcf
            else:
                variants_for_consensus = variants_norm_bcf

            # Index BCF
            subprocess.run(["bcftools", "index", str(variants_for_consensus)], check=False)

            consensus_cmd = ["bcftools", "consensus", "-f", str(reference), str(variants_for_consensus)]
            logger.info(f"Building draft consensus: {' '.join(consensus_cmd)} > {draft_fa}")
            with open(draft_fa, "w") as dfh:
                result = subprocess.run(consensus_cmd, stdout=dfh, capture_output=True, text=True)
            if result.returncode != 0 or not draft_fa.exists():
                logger.error(f"bcftools consensus failed: {result.stderr}")
                return None

            # 2) Polish with short reads using Pilon
            short_1 = reads.get("short_1")
            short_2 = reads.get("short_2")
            single = reads.get("single")

            # Build short-read alignment
            subprocess.run(["bwa", "index", str(draft_fa)], check=False)
            short_bam = output_dir / "short_aligned.bam"
            if short_1 and short_2:
                bwa_cmd = ["bwa", "mem", str(draft_fa), str(short_1), str(short_2)]
            elif single:
                bwa_cmd = ["bwa", "mem", str(draft_fa), str(single)]
            else:
                logger.error("Short reads are required for polishing with Pilon")
                return None

            logger.info(f"Aligning short reads to draft: {' '.join(bwa_cmd)} | samtools sort -o {short_bam}")
            bwa_proc = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE, text=False)
            sort2_proc = subprocess.run(["samtools", "sort", "-o", str(short_bam)], stdin=bwa_proc.stdout, capture_output=True)
            bwa_proc.stdout.close()  # type: ignore
            ret3 = bwa_proc.wait()
            if ret3 != 0 or sort2_proc.returncode != 0:
                logger.error(f"Short-read alignment failed: {sort2_proc.stderr.decode('utf-8', 'ignore')}")
                return None
            subprocess.run(["samtools", "index", str(short_bam)], check=False)

            polished_prefix = output_dir / "polished_assembly"
            pilon_cmd = [
                "pilon",
                "--genome", str(draft_fa),
                "--frags", str(short_bam),
                "--output", str(polished_prefix.name),
                "--outdir", str(output_dir),
            ]
            logger.info(f"Running Pilon: {' '.join(pilon_cmd)}")
            result = subprocess.run(pilon_cmd, capture_output=True, text=True)
            if result.returncode != 0:
                logger.error(f"Pilon failed: {result.stderr}")
                return None

            polished_fa = output_dir / f"{polished_prefix.name}.fasta"
            if not polished_fa.exists():
                # Some Pilon builds output .fa
                alt = output_dir / f"{polished_prefix.name}.fa"
                polished_fa = alt if alt.exists() else draft_fa

            return {
                "contigs": str(polished_fa),
                "scaffolds": str(polished_fa),
                "original_contigs": str(draft_fa),
                "long_alignment_bam": str(long_bam),
                "short_alignment_bam": str(short_bam)
            }
        except Exception as e:
            logger.error(f"Reference-guided long+short polishing pipeline failed: {e}")
            return None
    
    def _check_reference_hits(
        self, 
        contigs_file: Path, 
        reference: Union[str, Path]
    ) -> Dict[str, Union[bool, List[str], Dict]]:
        """
        Check if the assembly contains significant hits to the reference genome.
        
        Args:
            contigs_file: Path to assembled contigs
            reference: Reference genome file
            
        Returns:
            Dictionary containing hit information and matching contig IDs
        """
        try:
            # Use Minimap2 for better sequence alignment (faster and more accurate than BLAST)
            import tempfile
            import os
            
            # Create temporary directory for Minimap2
            temp_dir = tempfile.mkdtemp()
            ref_index = os.path.join(temp_dir, "reference.mmi")
            alignment_output = os.path.join(temp_dir, "alignment.paf")
            
            try:
                # Build Minimap2 index
                index_cmd = [
                    "minimap2",
                    "-d", ref_index,
                    str(reference)
                ]
                
                result = subprocess.run(index_cmd, capture_output=True, text=True)
                if result.returncode != 0:
                    logger.warning(f"Failed to build Minimap2 index: {result.stderr}")
                    return {
                        'has_hits': False,
                        'num_hits': 0,
                        'matching_contigs': [],
                        'alignment_results': []
                    }
                
                # Run Minimap2 alignment
                align_cmd = [
                    "minimap2",
                    "-c",  # Output CIGAR instead of MD
                    "-N", "1",  # Report at most 1 secondary alignment
                    "--secondary=no",  # No secondary alignments
                    "-x", "map-ont" if self.rna_mode else "asm",  # Preset for RNA or DNA
                    ref_index,
                    str(contigs_file)
                ]
                
                result = subprocess.run(align_cmd, capture_output=True, text=True)
                if result.returncode != 0:
                    logger.warning(f"Minimap2 alignment failed: {result.stderr}")
                    return {
                        'has_hits': False,
                        'num_hits': 0,
                        'matching_contigs': [],
                        'alignment_results': []
                    }
                
                # Parse Minimap2 results
                matching_contigs = []
                alignment_results = []
                
                for line in result.stdout.strip().split('\n'):
                    if line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 12:
                            contig_id = parts[0]
                            contig_length = int(parts[1])
                            alignment_start = int(parts[2])
                            alignment_end = int(parts[3])
                            strand = parts[4]
                            ref_name = parts[5]
                            ref_length = int(parts[6])
                            ref_start = int(parts[7])
                            ref_end = int(parts[8])
                            matches = int(parts[9])
                            alignment_length = int(parts[10])
                            quality = int(parts[11])
                            
                            # Calculate identity and coverage
                            identity = matches / alignment_length if alignment_length > 0 else 0
                            coverage = alignment_length / contig_length if contig_length > 0 else 0
                            
                            # Filter for significant hits: >80% identity and >50% coverage
                            if identity > 0.8 and coverage > 0.5:
                                matching_contigs.append(contig_id)
                                alignment_results.append({
                                    'contig_id': contig_id,
                                    'identity': identity,
                                    'coverage': coverage,
                                    'alignment_length': alignment_length,
                                    'matches': matches,
                                    'contig_length': contig_length
                                })
                
                logger.info(f"Found {len(matching_contigs)} contigs with significant alignment to reference")
                
                return {
                    'has_hits': len(matching_contigs) > 0,
                    'num_hits': len(matching_contigs),
                    'matching_contigs': matching_contigs,
                    'alignment_results': alignment_results
                }
                
            finally:
                # Clean up temporary files
                import shutil
                shutil.rmtree(temp_dir, ignore_errors=True)
                
        except Exception as e:
            logger.warning(f"Error checking reference hits: {e}")
            return {
                'has_hits': False,
                'num_hits': 0,
                'matching_contigs': [],
                'alignment_results': []
            }
    
    def _filter_contigs_by_reference_hits(
        self,
        contigs_file: Path,
        matching_contigs: List[str],
        output_file: Path
    ) -> bool:
        """
        Filter contigs to only include those that match the reference genome.
        
        Args:
            contigs_file: Path to original contigs file
            matching_contigs: List of contig IDs that match the reference
            output_file: Path to save filtered contigs
            
        Returns:
            True if filtering successful, False otherwise
        """
        try:
            logger.info(f"Filtering contigs to include only {len(matching_contigs)} reference-matching contigs")
            
            # Load all contigs
            all_contigs = list(SeqIO.parse(contigs_file, "fasta"))
            
            # Filter to only matching contigs
            filtered_contigs = []
            for contig in all_contigs:
                if contig.id in matching_contigs:
                    filtered_contigs.append(contig)
            
            if not filtered_contigs:
                logger.warning("No matching contigs found after filtering")
                return False
            
            # Write filtered contigs
            SeqIO.write(filtered_contigs, output_file, "fasta")
            logger.info(f"Saved {len(filtered_contigs)} reference-matching contigs to {output_file}")
            
            return True
            
        except Exception as e:
            logger.error(f"Error filtering contigs by reference hits: {e}")
            return False
    
    def _evaluate_reference_assembly(
        self, 
        assembly_results: Dict[str, str], 
        reference: Union[str, Path]
    ) -> Dict[str, Union[str, List[str], Dict]]:
        """
        Evaluate the reference-guided assembly results.
        
        Args:
            assembly_results: Dictionary containing assembly files
            reference: Reference genome file
            
        Returns:
            Dictionary containing evaluation results
        """
        logger.info("Evaluating reference-guided assembly")
        
        # Run CheckV with reference
        quality_assessment_dir = self.output_dir / "04_quality_assessment"
        quality_assessment_dir.mkdir(parents=True, exist_ok=True)
        
        checkv_results = self.validator.run_checkv(
            assembly_results["contigs"],
            output_dir=str(quality_assessment_dir / "reference_checkv")
        )
        
        # Calculate basic statistics
        stats = self.validator.calculate_genome_statistics(assembly_results["contigs"])
        
        # Determine if assembly was successful
        success = (
            stats.get("num_contigs", 0) > 0
        )
        
        # Get reference hit information
        reference_hits = assembly_results.get("reference_hits", {})
        
        return {
            "status": "completed" if success else "failed",
            "assembly_type": "reference_guided",
            "reference_genome": str(reference),
            "assembly_dir": str(self.output_dir),
            "assembly_files": assembly_results,
            "checkv_results": checkv_results,
            "statistics": stats,
            "viral_genomes": [assembly_results["contigs"]] if success else [],
            "total_viral_contigs": stats.get("num_contigs", 0) if success else 0,
            "reference_matching_contigs": reference_hits.get("matching_contigs", []),
            "num_reference_hits": reference_hits.get("num_hits", 0),
            "alignment_results": reference_hits.get("alignment_results", []),
            "message": f"Reference genome successfully assembled and filtered - {len(reference_hits.get('matching_contigs', []))} contigs match reference" if success else "Reference genome not detected in sample"
        }
    
    def cleanup_temp_files(self):
        """Clean up temporary files and directories created during analysis."""
        try:
            # Find all temporary directories with viral_ prefix
            temp_patterns = [
                "/tmp/viral_*",
                "/tmp/viral_preprocess_*",
                "/tmp/viral_identify_*", 
                "/tmp/viral_validate_*",
                "/tmp/viral_genes_*"
            ]
            
            total_cleaned = 0
            for pattern in temp_patterns:
                temp_dirs = glob.glob(pattern)
                for temp_dir in temp_dirs:
                    try:
                        if os.path.exists(temp_dir) and os.path.isdir(temp_dir):
                            shutil.rmtree(temp_dir)
                            total_cleaned += 1
                            logger.debug(f"Cleaned up temporary directory: {temp_dir}")
                    except Exception as e:
                        logger.warning(f"Failed to clean up {temp_dir}: {e}")
            
            if total_cleaned > 0:
                logger.info(f"Cleaned up {total_cleaned} temporary directories")
            else:
                logger.debug("No temporary directories found to clean up")
                
        except Exception as e:
            logger.warning(f"Error during temporary file cleanup: {e}")
