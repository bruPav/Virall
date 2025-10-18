"""
Single-cell sequencing data preprocessing module.
Handles 10X Genomics and other single-cell data formats.
"""

import os
import subprocess
import tempfile
import shutil
import gzip
from pathlib import Path
from typing import Tuple, Optional, Dict, List, Set, Union
from loguru import logger

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .preprocessor import Preprocessor


class SingleCellPreprocessor(Preprocessor):
    """
    Specialized preprocessor for single-cell sequencing data.
    Extends the base Preprocessor with single-cell specific functionality.
    """
    
    def __init__(self, threads: int = 8, cellranger_path: str = None):
        """
        Initialize the single-cell preprocessor.
        
        Args:
            threads: Number of threads to use for processing
            cellranger_path: Path to Cell Ranger installation
        """
        super().__init__(threads)
        self.cellranger_path = cellranger_path or "cellranger"
        self.cell_barcodes = {}
        self.umi_length = 12  # Standard UMI length for 10X
        self.barcode_length = 16  # Standard barcode length for 10X
        
        logger.info(f"SingleCellPreprocessor initialized with {threads} threads")
    
    def process_10x_data(
        self, 
        r1_file: str, 
        r2_file: str, 
        i1_file: str = None,
        sample_id: str = "sample",
        output_dir: str = None,
        min_cells: int = 100
    ) -> Dict[str, str]:
        """
        Process 10X Genomics single-cell data.
        
        Args:
            r1_file: Path to R1 reads (barcodes + UMIs)
            r2_file: Path to R2 reads (cDNA sequences)
            i1_file: Path to I1 reads (sample indices, optional)
            sample_id: Sample identifier
            output_dir: Output directory for processed files
            min_cells: Minimum number of cells to process
            
        Returns:
            Dictionary containing processed file paths and cell information
        """
        logger.info(f"Processing 10X single-cell data for sample {sample_id}")
        
        if output_dir is None:
            output_dir = self.temp_dir / "single_cell_processed"
        else:
            output_dir = Path(output_dir)
        
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Step 1: Extract and validate cell barcodes
        logger.info("Step 1: Extracting cell barcodes from R1 reads")
        cell_barcodes = self.extract_cell_barcodes(r1_file)
        
        if len(cell_barcodes) < min_cells:
            logger.warning(f"Only {len(cell_barcodes)} cells found, minimum required: {min_cells}")
            logger.warning("Consider lowering min_cells threshold or checking data quality")
        
        # Step 2: Process reads with cell information
        logger.info("Step 2: Processing reads with cell barcode information")
        processed_reads = self._process_reads_with_cells(
            r1_file, r2_file, i1_file, cell_barcodes, output_dir
        )
        
        # Step 3: Generate cell metadata
        logger.info("Step 3: Generating cell metadata")
        cell_metadata = self._generate_cell_metadata(cell_barcodes, output_dir)
        
        # Step 4: Quality control
        logger.info("Step 4: Running single-cell quality control")
        qc_results = self._run_single_cell_qc(processed_reads, cell_metadata, output_dir)
        
        return {
            "processed_reads": processed_reads,
            "cell_metadata": cell_metadata,
            "qc_results": qc_results,
            "num_cells": len(cell_barcodes),
            "output_dir": str(output_dir)
        }
    
    def extract_cell_barcodes(self, r1_file: str) -> Dict[str, str]:
        """
        Extract and validate cell barcodes from R1 reads.
        
        Args:
            r1_file: Path to R1 FASTQ file
            
        Returns:
            Dictionary mapping cell barcodes to read counts
        """
        logger.info(f"Extracting cell barcodes from {r1_file}")
        
        barcode_counts = {}
        total_reads = 0
        
        # Open file (handle gzipped files)
        if str(r1_file).endswith('.gz'):
            file_handle = gzip.open(r1_file, "rt")
        else:
            file_handle = open(r1_file, "rt")
        
        try:
            for record in SeqIO.parse(file_handle, "fastq"):
                total_reads += 1
                
                # Extract barcode (first 16 bases for 10X)
                barcode = str(record.seq[:self.barcode_length])
                
                # Basic validation (only A, T, C, G)
                if all(base in 'ATCG' for base in barcode):
                    barcode_counts[barcode] = barcode_counts.get(barcode, 0) + 1
                
                # Limit processing for large files
                if total_reads > 1000000:  # Process first 1M reads for barcode extraction
                    break
        finally:
            file_handle.close()
        
        # Filter barcodes with sufficient reads (at least 10 reads)
        filtered_barcodes = {
            barcode: count for barcode, count in barcode_counts.items() 
            if count >= 10
        }
        
        logger.info(f"Extracted {len(filtered_barcodes)} cell barcodes from {total_reads} reads")
        self.cell_barcodes = filtered_barcodes
        
        return filtered_barcodes
    
    def _process_reads_with_cells(
        self, 
        r1_file: str, 
        r2_file: str, 
        i1_file: str,
        cell_barcodes: Dict[str, str],
        output_dir: Path
    ) -> Dict[str, str]:
        """
        Process reads while maintaining cell barcode information.
        
        Args:
            r1_file: Path to R1 reads
            r2_file: Path to R2 reads
            i1_file: Path to I1 reads (optional)
            cell_barcodes: Dictionary of valid cell barcodes
            output_dir: Output directory
            
        Returns:
            Dictionary containing processed file paths
        """
        logger.info("Processing reads with cell barcode information")
        
        # Create output files
        processed_r1 = output_dir / "processed_R1.fastq"
        processed_r2 = output_dir / "processed_R2.fastq"
        cell_assignments = output_dir / "cell_assignments.tsv"
        
        # Open input files
        r1_handle = gzip.open(r1_file, "rt") if str(r1_file).endswith('.gz') else open(r1_file, "rt")
        r2_handle = gzip.open(r2_file, "rt") if str(r2_file).endswith('.gz') else open(r2_file, "rt")
        
        # Open output files
        with open(processed_r1, "w") as out_r1, open(processed_r2, "w") as out_r2, \
             open(cell_assignments, "w") as out_assignments:
            
            # Write header for assignments file
            out_assignments.write("read_id\tcell_barcode\tumi\tis_valid_cell\n")
            
            processed_reads = 0
            valid_reads = 0
            
            try:
                for r1_record, r2_record in zip(
                    SeqIO.parse(r1_handle, "fastq"), 
                    SeqIO.parse(r2_handle, "fastq")
                ):
                    processed_reads += 1
                    
                    # Extract barcode and UMI from R1
                    barcode = str(r1_record.seq[:self.barcode_length])
                    umi = str(r1_record.seq[self.barcode_length:self.barcode_length + self.umi_length])
                    
                    # Check if this is a valid cell
                    is_valid_cell = barcode in cell_barcodes
                    
                    if is_valid_cell:
                        valid_reads += 1
                        
                        # Write processed reads
                        SeqIO.write(r1_record, out_r1, "fastq")
                        SeqIO.write(r2_record, out_r2, "fastq")
                        
                        # Write assignment
                        out_assignments.write(f"{r1_record.id}\t{barcode}\t{umi}\tTrue\n")
                    
                    # Limit processing for testing
                    if processed_reads > 100000:  # Process first 100K reads
                        break
                        
            finally:
                r1_handle.close()
                r2_handle.close()
        
        logger.info(f"Processed {processed_reads} reads, {valid_reads} from valid cells")
        
        return {
            "processed_r1": str(processed_r1),
            "processed_r2": str(processed_r2),
            "cell_assignments": str(cell_assignments)
        }
    
    def _generate_cell_metadata(
        self, 
        cell_barcodes: Dict[str, str], 
        output_dir: Path
    ) -> str:
        """
        Generate cell metadata file.
        
        Args:
            cell_barcodes: Dictionary of cell barcodes and read counts
            output_dir: Output directory
            
        Returns:
            Path to cell metadata file
        """
        logger.info("Generating cell metadata")
        
        metadata_file = output_dir / "cell_metadata.tsv"
        
        with open(metadata_file, "w") as f:
            f.write("cell_barcode\tread_count\tcell_id\n")
            
            for i, (barcode, count) in enumerate(cell_barcodes.items()):
                cell_id = f"cell_{i+1:06d}"
                f.write(f"{barcode}\t{count}\t{cell_id}\n")
        
        logger.info(f"Generated metadata for {len(cell_barcodes)} cells")
        return str(metadata_file)
    
    def _run_single_cell_qc(
        self, 
        processed_reads: Dict[str, str], 
        cell_metadata: str,
        output_dir: Path
    ) -> Dict[str, any]:
        """
        Run single-cell specific quality control.
        
        Args:
            processed_reads: Dictionary of processed read files
            cell_metadata: Path to cell metadata file
            output_dir: Output directory
            
        Returns:
            Dictionary containing QC results
        """
        logger.info("Running single-cell quality control")
        
        # Load cell metadata
        metadata_df = pd.read_csv(cell_metadata, sep='\t')
        
        # Calculate basic statistics
        total_cells = len(metadata_df)
        total_reads = metadata_df['read_count'].sum()
        avg_reads_per_cell = metadata_df['read_count'].mean()
        median_reads_per_cell = metadata_df['read_count'].median()
        
        # Generate QC report
        qc_file = output_dir / "single_cell_qc_report.txt"
        
        with open(qc_file, "w") as f:
            f.write("Single-Cell Quality Control Report\n")
            f.write("=" * 40 + "\n")
            f.write(f"Total cells: {total_cells}\n")
            f.write(f"Total reads: {total_reads:,}\n")
            f.write(f"Average reads per cell: {avg_reads_per_cell:.1f}\n")
            f.write(f"Median reads per cell: {median_reads_per_cell:.1f}\n")
            f.write(f"Min reads per cell: {metadata_df['read_count'].min()}\n")
            f.write(f"Max reads per cell: {metadata_df['read_count'].max()}\n")
        
        logger.info(f"QC report saved to {qc_file}")
        
        return {
            "qc_file": str(qc_file),
            "total_cells": total_cells,
            "total_reads": total_reads,
            "avg_reads_per_cell": avg_reads_per_cell,
            "median_reads_per_cell": median_reads_per_cell
        }
    
    def demultiplex_by_cell(
        self, 
        r1_file: str, 
        r2_file: str,
        cell_barcodes: Dict[str, str],
        output_dir: str
    ) -> Dict[str, str]:
        """
        Demultiplex reads by cell barcode (creates separate files per cell).
        
        Args:
            r1_file: Path to R1 reads
            r2_file: Path to R2 reads
            cell_barcodes: Dictionary of cell barcodes
            output_dir: Output directory for demultiplexed files
            
        Returns:
            Dictionary mapping cell barcodes to file paths
        """
        logger.info(f"Demultiplexing reads by cell barcode")
        
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create output files for each cell
        cell_files = {}
        for barcode in cell_barcodes.keys():
            cell_id = f"cell_{barcode}"
            cell_files[barcode] = {
                "r1": output_dir / f"{cell_id}_R1.fastq",
                "r2": output_dir / f"{cell_id}_R2.fastq"
            }
        
        # Open input files
        r1_handle = gzip.open(r1_file, "rt") if str(r1_file).endswith('.gz') else open(r1_file, "rt")
        r2_handle = gzip.open(r2_file, "rt") if str(r2_file).endswith('.gz') else open(r2_file, "rt")
        
        # Open output files
        file_handles = {}
        for barcode, files in cell_files.items():
            file_handles[barcode] = {
                "r1": open(files["r1"], "w"),
                "r2": open(files["r2"], "w")
            }
        
        try:
            for r1_record, r2_record in zip(
                SeqIO.parse(r1_handle, "fastq"), 
                SeqIO.parse(r2_handle, "fastq")
            ):
                # Extract barcode
                barcode = str(r1_record.seq[:self.barcode_length])
                
                if barcode in cell_barcodes:
                    # Write to appropriate cell files
                    SeqIO.write(r1_record, file_handles[barcode]["r1"], "fastq")
                    SeqIO.write(r2_record, file_handles[barcode]["r2"], "fastq")
        
        finally:
            # Close all files
            r1_handle.close()
            r2_handle.close()
            for handles in file_handles.values():
                handles["r1"].close()
                handles["r2"].close()
        
        logger.info(f"Demultiplexed reads for {len(cell_files)} cells")
        return cell_files
    
    def cleanup(self) -> None:
        """Clean up temporary files."""
        super().cleanup()
        logger.info("SingleCellPreprocessor cleanup completed")
