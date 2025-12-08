"""
Data preprocessing module for quality control, trimming, and error correction.
"""

import os
import subprocess
import tempfile
import shutil
import gzip
from pathlib import Path
from typing import Tuple, Optional, Dict, Optional, Dict, List, Set
from loguru import logger

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class Preprocessor:
    """
    Handles preprocessing of sequencing reads including quality control,
    adapter trimming, and error correction for both short and long reads.
    """
    
    def __init__(self, threads: int = 8, config: Optional[Dict] = None):
        """
        Initialize the preprocessor.
        
        Args:
            threads: Number of threads to use for processing
            config: Optional configuration dictionary
        """
        self.threads = threads
        self.temp_dir = Path(tempfile.mkdtemp(prefix="viral_preprocess_"))
        self.config = config or {}
        
        logger.info(f"Preprocessor initialized with {threads} threads")
    
    def process_paired_reads(
        self, 
        reads_1: str, 
        reads_2: str,
        quality_threshold: Optional[int] = None,
        min_length: Optional[int] = None
    ) -> Tuple[str, str]:
        """
        Process paired-end reads with quality control and trimming.
        
        Args:
            reads_1: Path to first mate reads
            reads_2: Path to second mate reads
            quality_threshold: Minimum quality score threshold (default: 20, or from config)
            min_length: Minimum read length after trimming (default: 50, or from config)
            
        Returns:
            Tuple of (processed_reads_1, processed_reads_2)
        """
        logger.info(f"Processing paired-end reads: {reads_1}, {reads_2}")
        
        # Get thresholds from config if not provided
        if quality_threshold is None:
            quality_threshold = (
                self.config.get("quality_threshold") or
                self.config.get("quality_control", {}).get("quality_threshold", 20)
            )
        if min_length is None:
            min_length = (
                self.config.get("short_read_min_length") or
                self.config.get("quality_control", {}).get("min_read_length", 50)
            )
        
        logger.info(f"Using quality_threshold={quality_threshold}, min_length={min_length} for short reads")
        
        # Quality control with FastQC
        self._run_fastqc(reads_1, reads_2)
        
        # Trimming with fastp (automatic adapter detection)
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
        quality_threshold: Optional[int] = None,
        min_length: Optional[int] = None
    ) -> str:
        """
        Process long reads (PacBio/ONT) with quality control and trimming using fastplong.
        
        Args:
            reads: Path to long reads
            quality_threshold: Minimum quality score threshold (default: 7, or from config)
            min_length: Minimum read length after filtering (default: 1000, or from config)
            
        Returns:
            Path to processed reads
        """
        logger.info(f"Processing long reads: {reads}")
        
        # Get thresholds from config if not provided
        # Support both flat and nested config structures
        if quality_threshold is None:
            quality_threshold = (
                self.config.get("long_read_quality_threshold") or
                self.config.get("quality_control", {}).get("long_read_quality_threshold", 7)
            )
        if min_length is None:
            min_length = (
                self.config.get("long_read_min_length") or
                self.config.get("quality_control", {}).get("long_read_min_length", 1000)
            )
        
        logger.info(f"Using quality_threshold={quality_threshold}, min_length={min_length} for long reads")
        
        # Use fastplong for all-in-one QC, adapter detection, and trimming
        try:
            trimmed_reads = self._trim_long_reads_fastplong(reads, quality_threshold, min_length)
        except Exception as e:
            logger.error(f"fastplong failed ({e})")
            # If fastplong fails, return original reads (no trimming)
            logger.warning("Skipping trimming - using original reads")
            trimmed_reads = reads
        
        # Quality filtering
        trimmed_reads = self._filter_long_reads(
            trimmed_reads, quality_threshold, min_length
        )
        
        # Error correction (optional)
        if self._should_correct_long_reads():
            corrected_reads = self._correct_long_reads(trimmed_reads)
            return corrected_reads
        
        return trimmed_reads
    
    def process_single_reads(
        self, 
        reads: str,
        quality_threshold: Optional[int] = None,
        min_length: Optional[int] = None
    ) -> str:
        """
        Process single-end reads with quality control and trimming.
        
        Args:
            reads: Path to single-end reads
            quality_threshold: Minimum quality score threshold (default: 20, or from config)
            min_length: Minimum read length after trimming (default: 50, or from config)
            
        Returns:
            Path to processed reads
        """
        logger.info(f"Processing single-end reads: {reads}")
        
        # Get thresholds from config if not provided
        if quality_threshold is None:
            quality_threshold = (
                self.config.get("quality_threshold") or
                self.config.get("quality_control", {}).get("quality_threshold", 20)
            )
        if min_length is None:
            min_length = (
                self.config.get("short_read_min_length") or
                self.config.get("quality_control", {}).get("min_read_length", 50)
            )
        
        logger.info(f"Using quality_threshold={quality_threshold}, min_length={min_length} for single reads")
        
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
        """Trim paired-end reads with fastp (automatic adapter detection)."""
        logger.info("Trimming paired-end reads with fastp (automatic adapter detection)")
        
        output_1 = self.temp_dir / "trimmed_1.fastq.gz"
        output_2 = self.temp_dir / "trimmed_2.fastq.gz"
        output_1_plain = self.temp_dir / "trimmed_1.fastq"
        output_2_plain = self.temp_dir / "trimmed_2.fastq"
        
        # fastp HTML and JSON reports
        fastp_html = self.temp_dir / "fastp_report.html"
        fastp_json = self.temp_dir / "fastp_report.json"
        
        # Build fastp command with automatic adapter detection
        cmd = [
            "fastp",
            "-i", reads_1,
            "-I", reads_2,
            "-o", str(output_1),
            "-O", str(output_2),
            "--thread", str(self.threads),
            "--qualified_quality_phred", str(quality_threshold),
            "--length_required", str(min_length),
            "--html", str(fastp_html),
            "--json", str(fastp_json),
            # Automatic adapter detection and trimming
            "--detect_adapter_for_pe",
            # Quality filtering (qualified_quality_phred enables quality filtering)
            # Per-read quality cutting
            "--cut_front",  # Remove low quality bases from 5' end
            "--cut_tail",   # Remove low quality bases from 3' end
            "--cut_mean_quality", str(quality_threshold),
            "--cut_window_size", "4",
            # Overlap analysis and correction for PE data
            "--correction",
            "--overlap_len_require", "30",
            "--overlap_diff_limit", "5",
            # PolyG tail trimming (common in Illumina)
            "--trim_poly_g",
            "--poly_g_min_len", "10"
        ]
        
        # Check if fastp is available
        if not self._check_tool_available("fastp"):
            raise RuntimeError("fastp is not available. Please install fastp for read trimming.")
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"fastp failed: {result.stderr}")
        
        logger.info("fastp paired-end trimming completed (with automatic adapter detection)")
        logger.info(f"Quality reports saved: {fastp_html}, {fastp_json}")

        # Also write plain FASTQ copies alongside gzipped outputs
        try:
            with gzip.open(output_1, "rt") as fin, open(output_1_plain, "w") as fout:
                shutil.copyfileobj(fin, fout)
            with gzip.open(output_2, "rt") as fin, open(output_2_plain, "w") as fout:
                shutil.copyfileobj(fin, fout)
        except Exception as e:
            logger.warning(f"Failed to create plain FASTQ copies: {e}")
        return str(output_1), str(output_2)
    
    def _trim_single_reads(
        self, 
        reads: str,
        quality_threshold: int,
        min_length: int
    ) -> str:
        """Trim single-end reads with fastp (automatic adapter detection)."""
        logger.info("Trimming single-end reads with fastp (automatic adapter detection)")
        
        output = self.temp_dir / "trimmed_single.fastq.gz"
        output_plain = self.temp_dir / "trimmed_single.fastq"
        
        # fastp HTML and JSON reports
        fastp_html = self.temp_dir / "fastp_report_single.html"
        fastp_json = self.temp_dir / "fastp_report_single.json"
        
        # Build fastp command with automatic adapter detection
        cmd = [
            "fastp",
            "-i", reads,
            "-o", str(output),
            "--thread", str(self.threads),
            "--qualified_quality_phred", str(quality_threshold),
            "--length_required", str(min_length),
            "--html", str(fastp_html),
            "--json", str(fastp_json),
            # Automatic adapter detection (enabled by default for single-end reads)
            # Quality filtering (qualified_quality_phred enables quality filtering)
            # Per-read quality cutting
            "--cut_front",  # Remove low quality bases from 5' end
            "--cut_tail",   # Remove low quality bases from 3' end
            "--cut_mean_quality", str(quality_threshold),
            "--cut_window_size", "4",
            # PolyG tail trimming (common in Illumina)
            "--trim_poly_g",
            "--poly_g_min_len", "10"
        ]
        
        # Check if fastp is available
        if not self._check_tool_available("fastp"):
            raise RuntimeError("fastp is not available. Please install fastp for read trimming.")
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"fastp failed: {result.stderr}")
        
        logger.info("fastp single-end trimming completed (with automatic adapter detection)")
        logger.info(f"Quality reports saved: {fastp_html}, {fastp_json}")

        # Also write plain FASTQ copy alongside gzipped output
        try:
            with gzip.open(output, "rt") as fin, open(output_plain, "w") as fout:
                shutil.copyfileobj(fin, fout)
        except Exception as e:
            logger.warning(f"Failed to create plain FASTQ copy: {e}")
        return str(output)
    
    def _trim_long_reads_fastplong(
        self,
        reads: str,
        quality_threshold: int,
        min_length: int
    ) -> str:
        """Trim long reads using fastplong (automatic adapter detection for ONT/PacBio)."""
        logger.info("Trimming long reads with fastplong (automatic adapter detection)")
        
        output = self.temp_dir / "trimmed_long.fastq"
        output_gz = self.temp_dir / "trimmed_long.fastq.gz"
        
        # fastplong HTML and JSON reports
        fastplong_html = self.temp_dir / "fastplong_report.html"
        fastplong_json = self.temp_dir / "fastplong_report.json"
        
        # Build fastplong command with automatic adapter detection
        cmd = [
            "fastplong",
            "-i", reads,
            "-o", str(output),
            "--thread", str(self.threads),
            "--qualified_quality_phred", str(quality_threshold),
            "--length_required", str(min_length),
            "--html", str(fastplong_html),
            "--json", str(fastplong_json),
            # Automatic adapter detection is enabled by default in fastplong
            # Quality filtering
            # Per-read quality cutting
            "--cut_front",  # Remove low quality bases from 5' end
            "--cut_tail",   # Remove low quality bases from 3' end
            "--cut_mean_quality", str(quality_threshold),
            "--cut_window_size", "10"  # Larger window for long reads
        ]
        
        # Check if fastplong is available
        if not self._check_tool_available("fastplong"):
            logger.error("fastplong not available")
            raise RuntimeError("fastplong not available")
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"fastplong failed: {result.stderr}")
            raise RuntimeError(f"fastplong failed: {result.stderr}")
        
        logger.info("fastplong trimming completed (with automatic adapter detection)")
        logger.info(f"Quality reports saved: {fastplong_html}, {fastplong_json}")

        # Fix caption mismatches: fastplong adds prefixes to headers when splitting reads
        # but doesn't update the + line caption, causing Flye errors
        logger.info("Fixing FASTQ caption mismatches from fastplong...")
        self._fix_fastq_captions(output)
        
        # Also create gzipped version
        try:
            with open(output, "rb") as fin, gzip.open(output_gz, "wb") as fout:
                fout.writelines(fin)
        except Exception as e:
            logger.warning(f"Failed to create gzipped version: {e}")
        
        return str(output)
    
    def _should_trim_ont_by_detection(self, reads: str) -> bool:
        """Detect if ONT reads should be trimmed based on filename heuristic.

        Returns True when filename suggests ONT reads.
        """
        # Use filename heuristic to detect ONT reads
        return self._is_ont_reads(reads)
    
    def _filter_long_reads(
        self, 
        reads: str,
        quality_threshold: int,
        min_length: int
    ) -> str:
        """Filter long reads by quality and length."""
        logger.info("Filtering long reads by quality and length")
        
        output = self.temp_dir / "filtered_long.fastq"
        output_gz = self.temp_dir / "filtered_long.fastq.gz"
        filtered_records = []
        
        # Open input with gzip if needed
        if str(reads).endswith('.gz'):
            in_handle = gzip.open(reads, "rt")
        else:
            in_handle = open(reads, "rt")
        for record in SeqIO.parse(in_handle, "fastq"):
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
        
        # Write filtered reads as plain FASTQ
        with open(output, "w") as f:
            SeqIO.write(filtered_records, f, "fastq")

        # Also write gzipped FASTQ
        try:
            with gzip.open(output_gz, "wt") as f:
                SeqIO.write(filtered_records, f, "fastq")
        except Exception as e:
            logger.warning(f"Failed to write gzipped filtered long reads: {e}")
        
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
    
    def _fix_fastq_captions(self, fastq_file: Path) -> None:
        """
        Fix FASTQ caption mismatches caused by fastplong.
        
        When fastplong splits reads by adapters, it adds prefixes like
        'split-by-adapter-right-' to the header but doesn't update the
        + line caption, causing Flye to fail with "Sequence and quality
        captions differ" error.
        
        This method normalizes the FASTQ by removing captions from +
        lines (making them just '+') when they don't match the header.
        """
        try:
            fixed_count = 0
            total_reads = 0
            
            # Create temporary file for output
            temp_output = fastq_file.with_suffix('.fixed.tmp')
            
            with open(fastq_file, 'r') as fin, open(temp_output, 'w') as fout:
                lines = []
                for line in fin:
                    lines.append(line)
                    
                    # Process in groups of 4 lines (one read)
                    if len(lines) == 4:
                        header = lines[0].rstrip('\n\r')
                        sequence = lines[1].rstrip('\n\r')
                        plus_line = lines[2].rstrip('\n\r')
                        quality = lines[3].rstrip('\n\r')
                        
                        # Extract captions
                        header_caption = header[1:] if header.startswith('@') else None
                        plus_caption = plus_line[1:] if plus_line.startswith('+') and len(plus_line) > 1 else None
                        
                        # Check if captions mismatch
                        if plus_caption and header_caption and plus_caption != header_caption:
                            # Remove caption from + line (make it just '+')
                            fout.write(header + '\n')
                            fout.write(sequence + '\n')
                            fout.write('+\n')  # No caption
                            fout.write(quality + '\n')
                            fixed_count += 1
                        else:
                            # Write as-is
                            for l in lines:
                                fout.write(l)
                        
                        lines = []
                        total_reads += 1
                
                # Handle any remaining lines (shouldn't happen in valid FASTQ)
                if lines:
                    logger.warning(f"Found {len(lines)} incomplete lines at end of FASTQ file")
                    for l in lines:
                        fout.write(l)
            
            # Replace original with fixed version
            shutil.move(str(temp_output), str(fastq_file))
            
            if fixed_count > 0:
                logger.info(f"Fixed {fixed_count} reads with caption mismatches out of {total_reads} total reads")
            else:
                logger.debug("No caption mismatches found in FASTQ file")
                
        except Exception as e:
            logger.warning(f"Failed to fix FASTQ captions: {e}. Original file will be used.")
    
    def _check_tool_available(self, tool: str) -> bool:
        """Check if a tool is available in PATH."""
        result = subprocess.run(
            ["which", tool],
            capture_output=True,
            text=True
        )
        return result.returncode == 0
    
    def get_quality_report(self) -> Dict:
        """Get quality control report from FastQC/fastp/fastplong results."""
        report = {}
        
        # Check for fastp reports (short reads)
        for file in self.temp_dir.glob("fastp_report*.html"):
            report[file.stem] = {
                "status": "completed",
                "report_file": str(file),
                "type": "fastp"
            }
        
        # Check for fastplong reports (long reads)
        for file in self.temp_dir.glob("fastplong_report*.html"):
            report[file.stem] = {
                "status": "completed",
                "report_file": str(file),
                "type": "fastplong"
            }
        
        # Also include FastQC reports if present
        for file in self.temp_dir.glob("*.html"):
            if "fastp" not in file.stem and "fastplong" not in file.stem:
                report[file.stem] = {
                    "status": "completed",
                    "report_file": str(file),
                    "type": "fastqc"
                }
        
        return report
    
    def cleanup(self) -> None:
        """Clean up temporary files."""
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
            logger.info("Cleaned up temporary files")
