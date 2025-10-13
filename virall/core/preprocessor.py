"""
Data preprocessing module for quality control, trimming, and error correction.
"""

import os
import subprocess
import tempfile
import shutil
import gzip
from pathlib import Path
from typing import Tuple, Optional, Dict, List
from loguru import logger

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class Preprocessor:
    """
    Handles preprocessing of sequencing reads including quality control,
    adapter trimming, and error correction for both short and long reads.
    """
    
    def __init__(self, threads: int = 8):
        """
        Initialize the preprocessor.
        
        Args:
            threads: Number of threads to use for processing
        """
        self.threads = threads
        self.temp_dir = Path(tempfile.mkdtemp(prefix="viral_preprocess_"))
        
        logger.info(f"Preprocessor initialized with {threads} threads")
    
    def process_paired_reads(
        self, 
        reads_1: str, 
        reads_2: str,
        quality_threshold: int = 20,
        min_length: int = 50
    ) -> Tuple[str, str]:
        """
        Process paired-end reads with quality control and trimming.
        
        Args:
            reads_1: Path to first mate reads
            reads_2: Path to second mate reads
            quality_threshold: Minimum quality score threshold
            min_length: Minimum read length after trimming
            
        Returns:
            Tuple of (processed_reads_1, processed_reads_2)
        """
        logger.info(f"Processing paired-end reads: {reads_1}, {reads_2}")
        
        # Quality control with FastQC
        self._run_fastqc(reads_1, reads_2)
        
        # Trimming with Trimmomatic
        trimmed_1, trimmed_2 = self._trim_paired_reads(
            reads_1, reads_2, quality_threshold, min_length
        )
        
        # Error correction (optional)
        if self._should_correct_errors(reads_1):
            corrected_1, corrected_2 = self._correct_paired_reads(trimmed_1, trimmed_2)
            return corrected_1, corrected_2
        
        return trimmed_1, trimmed_2
    
    def process_long_reads(
        self, 
        reads: str,
        quality_threshold: int = 7,
        min_length: int = 1000
    ) -> str:
        """
        Process long reads (PacBio/ONT) with quality control and trimming.
        
        Args:
            reads: Path to long reads
            quality_threshold: Minimum quality score threshold
            min_length: Minimum read length after filtering
            
        Returns:
            Path to processed reads
        """
        logger.info(f"Processing long reads: {reads}")
        
        # Quality control
        self._run_fastqc(reads)
        
        # Adapter trimming (for ONT reads)
        if self._is_ont_reads(reads):
            trimmed_reads = self._trim_ont_adapters(reads)
        else:
            trimmed_reads = reads
        
        # Quality filtering
        filtered_reads = self._filter_long_reads(
            trimmed_reads, quality_threshold, min_length
        )
        
        # Error correction (optional)
        if self._should_correct_long_reads():
            corrected_reads = self._correct_long_reads(filtered_reads)
            return corrected_reads
        
        return filtered_reads
    
    def process_single_reads(
        self, 
        reads: str,
        quality_threshold: int = 20,
        min_length: int = 50
    ) -> str:
        """
        Process single-end reads with quality control and trimming.
        
        Args:
            reads: Path to single-end reads
            quality_threshold: Minimum quality score threshold
            min_length: Minimum read length after trimming
            
        Returns:
            Path to processed reads
        """
        logger.info(f"Processing single-end reads: {reads}")
        
        # Quality control
        self._run_fastqc(reads)
        
        # Trimming
        trimmed_reads = self._trim_single_reads(reads, quality_threshold, min_length)
        
        # Error correction (optional)
        if self._should_correct_errors(reads):
            corrected_reads = self._correct_single_reads(trimmed_reads)
            return corrected_reads
        
        return trimmed_reads
    
    def _run_fastqc(self, *read_files: str) -> None:
        """Run FastQC quality control on read files."""
        logger.info("Running FastQC quality control")
        
        for read_file in read_files:
            if not os.path.exists(read_file):
                logger.warning(f"Read file not found: {read_file}")
                continue
                
            cmd = [
                "fastqc",
                "--threads", str(self.threads),
                "--outdir", str(self.temp_dir),
                read_file
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                logger.warning(f"FastQC failed for {read_file}: {result.stderr}")
            else:
                logger.info(f"FastQC completed for {read_file}")
    
    def _trim_paired_reads(
        self, 
        reads_1: str, 
        reads_2: str,
        quality_threshold: int,
        min_length: int
    ) -> Tuple[str, str]:
        """Trim paired-end reads with Trimmomatic."""
        logger.info("Trimming paired-end reads with Trimmomatic")
        
        output_1 = self.temp_dir / "trimmed_1.fastq.gz"
        output_2 = self.temp_dir / "trimmed_2.fastq.gz"
        unpaired_1 = self.temp_dir / "unpaired_1.fastq.gz"
        unpaired_2 = self.temp_dir / "unpaired_2.fastq.gz"
        
        cmd = [
            "trimmomatic",
            "PE",
            "-threads", str(self.threads),
            reads_1, reads_2,
            str(output_1), str(unpaired_1),
            str(output_2), str(unpaired_2),
            f"LEADING:{quality_threshold}",
            f"TRAILING:{quality_threshold}",
            f"SLIDINGWINDOW:4:{quality_threshold}",
            f"MINLEN:{min_length}",
            "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10"
        ]
        
        # Set Java heap space for Trimmomatic
        env = os.environ.copy()
        env['_JAVA_OPTIONS'] = '-Xmx8g'
        
        result = subprocess.run(cmd, capture_output=True, text=True, env=env)
        if result.returncode != 0:
            raise RuntimeError(f"Trimmomatic failed: {result.stderr}")
        
        logger.info("Paired-end trimming completed")
        return str(output_1), str(output_2)
    
    def _trim_single_reads(
        self, 
        reads: str,
        quality_threshold: int,
        min_length: int
    ) -> str:
        """Trim single-end reads with Trimmomatic."""
        logger.info("Trimming single-end reads with Trimmomatic")
        
        output = self.temp_dir / "trimmed_single.fastq.gz"
        
        cmd = [
            "trimmomatic",
            "SE",
            "-threads", str(self.threads),
            reads,
            str(output),
            f"LEADING:{quality_threshold}",
            f"TRAILING:{quality_threshold}",
            f"SLIDINGWINDOW:4:{quality_threshold}",
            f"MINLEN:{min_length}",
            "ILLUMINACLIP:TruSeq3-SE.fa:2:30:10"
        ]
        
        # Set Java heap space for Trimmomatic
        env = os.environ.copy()
        env['_JAVA_OPTIONS'] = '-Xmx8g'
        
        result = subprocess.run(cmd, capture_output=True, text=True, env=env)
        if result.returncode != 0:
            raise RuntimeError(f"Trimmomatic failed: {result.stderr}")
        
        logger.info("Single-end trimming completed")
        return str(output)
    
    def _trim_ont_adapters(self, reads: str) -> str:
        """Trim ONT adapters using Porechop."""
        logger.info("Trimming ONT adapters with Porechop")
        
        output = self.temp_dir / "trimmed_ont.fastq.gz"
        
        cmd = [
            "porechop",
            "-i", reads,
            "-o", str(output),
            "--threads", str(self.threads)
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.warning(f"Porechop failed: {result.stderr}")
            return reads  # Return original if trimming fails
        
        logger.info("ONT adapter trimming completed")
        return str(output)
    
    def _filter_long_reads(
        self, 
        reads: str,
        quality_threshold: int,
        min_length: int
    ) -> str:
        """Filter long reads by quality and length."""
        logger.info("Filtering long reads by quality and length")
        
        output = self.temp_dir / "filtered_long.fastq.gz"
        filtered_records = []
        
        for record in SeqIO.parse(reads, "fastq"):
            # Simple quality filtering based on average quality
            if len(record.seq) >= min_length:
                # Calculate average quality
                if hasattr(record, "letter_annotations") and "phred_quality" in record.letter_annotations:
                    avg_quality = sum(record.letter_annotations["phred_quality"]) / len(record.letter_annotations["phred_quality"])
                    if avg_quality >= quality_threshold:
                        filtered_records.append(record)
                else:
                    # If no quality scores, keep the read
                    filtered_records.append(record)
        
        # Write filtered reads
        with gzip.open(output, "wt") as f:
            SeqIO.write(filtered_records, f, "fastq")
        
        logger.info(f"Filtered {len(filtered_records)} long reads")
        return str(output)
    
    def _correct_paired_reads(self, reads_1: str, reads_2: str) -> Tuple[str, str]:
        """Error correction for paired-end reads using Rcorrector."""
        logger.info("Correcting errors in paired-end reads")
        
        # For now, return trimmed reads (error correction can be added later)
        return reads_1, reads_2
    
    def _correct_single_reads(self, reads: str) -> str:
        """Error correction for single-end reads."""
        logger.info("Correcting errors in single-end reads")
        
        # For now, return trimmed reads (error correction can be added later)
        return reads
    
    def _correct_long_reads(self, reads: str) -> str:
        """Error correction for long reads using Nanopolish."""
        logger.info("Correcting errors in long reads with Nanopolish")
        
        # For now, return filtered reads (error correction can be added later)
        return reads
    
    def _should_correct_errors(self, reads: str) -> bool:
        """Determine if error correction should be performed."""
        # Simple heuristic: correct if file size is reasonable
        file_size = os.path.getsize(reads) / (1024 * 1024)  # MB
        return file_size < 1000  # Only correct if < 1GB
    
    def _should_correct_long_reads(self) -> bool:
        """Determine if long read error correction should be performed."""
        return False  # Disabled for now due to computational cost
    
    def _is_ont_reads(self, reads: str) -> bool:
        """Determine if reads are from Oxford Nanopore Technology."""
        # Simple heuristic: check file name or first few reads
        filename = Path(reads).name.lower()
        return "ont" in filename or "nanopore" in filename or "minion" in filename
    
    def get_quality_report(self) -> Dict:
        """Get quality control report from FastQC results."""
        report = {}
        
        # Parse FastQC results (simplified)
        fastqc_dir = self.temp_dir
        for file in fastqc_dir.glob("*.html"):
            report[file.stem] = {
                "status": "completed",
                "report_file": str(file)
            }
        
        return report
    
    def cleanup(self) -> None:
        """Clean up temporary files."""
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
            logger.info("Cleaned up temporary files")
