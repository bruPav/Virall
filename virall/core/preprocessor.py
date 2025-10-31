"""
Data preprocessing module for quality control, trimming, and error correction.
"""

import os
import subprocess
import tempfile
import shutil
import gzip
from pathlib import Path
from typing import Tuple, Optional, Dict, List, Set
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
        quality_threshold: int = 7,
        min_length: int = 1000
    ) -> str:
        """
        Process long reads (PacBio/ONT) with quality control and trimming using fastplong.
        
        Args:
            reads: Path to long reads
            quality_threshold: Minimum quality score threshold
            min_length: Minimum read length after filtering
            
        Returns:
            Path to processed reads
        """
        logger.info(f"Processing long reads: {reads}")
        
        # Use fastplong for all-in-one QC, adapter detection, and trimming
        try:
            trimmed_reads = self._trim_long_reads_fastplong(reads, quality_threshold, min_length)
        except Exception as e:
            logger.warning(f"fastplong failed ({e}); falling back to Porechop/manual filtering")
            # Fallback to old method
            trimmed_reads = reads
            try:
                if self._should_trim_ont_by_detection(reads):
                    trimmed_reads = self._trim_ont_adapters(reads)
            except Exception as e2:
                logger.warning(f"ONT adapter detection also failed ({e2}); continuing without adapter trimming")
            
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
        
        # Check if fastplong is available, fallback to Porechop
        if not self._check_tool_available("fastplong"):
            logger.warning("fastplong not available, will fall back to Porechop/manual filtering")
            raise RuntimeError("fastplong not available")
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"fastplong failed: {result.stderr}")
            raise RuntimeError(f"fastplong failed: {result.stderr}")
        
        logger.info("fastplong trimming completed (with automatic adapter detection)")
        logger.info(f"Quality reports saved: {fastplong_html}, {fastplong_json}")

        # Also create gzipped version
        try:
            with open(output, "rb") as fin, gzip.open(output_gz, "wb") as fout:
                fout.writelines(fin)
        except Exception as e:
            logger.warning(f"Failed to create gzipped version: {e}")
        
        return str(output)
    
    def _trim_ont_adapters(self, reads: str) -> str:
        """Fallback: Trim ONT adapters using Porechop (write plain + gz FASTQ)."""
        logger.info("Using Porechop fallback for ONT adapter trimming")
        
        output = self.temp_dir / "trimmed_ont.fastq"
        output_gz = self.temp_dir / "trimmed_ont.fastq.gz"
        
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
        
        logger.info("ONT adapter trimming completed with Porechop")

        # Also gzip a copy
        try:
            with open(output, "rt") as fin, gzip.open(output_gz, "wt") as fout:
                shutil.copyfileobj(fin, fout)
        except Exception as e:
            logger.warning(f"Failed to create gzipped ONT-trimmed copy: {e}")
        return str(output)

    def _should_trim_ont_by_detection(self, reads: str) -> bool:
        """Detect ONT adapters using Porechop's adapter catalog.

        Returns True when known ONT adapter motifs are observed in a sample of reads.
        Falls back to filename heuristic if the catalog is unavailable.
        """
        # Quick hint from filename
        if self._is_ont_reads(reads):
            return True
        adapters = self._load_known_ont_adapters()
        if not adapters:
            return False
        return self._reads_contain_any_sequence(reads, adapters, max_reads=2000)

    def _load_known_ont_adapters(self) -> Set[str]:
        """Load ONT adapter sequences from Porechop's adapters.py if available."""
        try:
            import importlib
            adapters_mod = importlib.import_module("porechop.adapters")
        except Exception:
            return set()
        candidates: Set[str] = set()

        def add_seq(s: str) -> None:
            s_clean = s.strip().upper()
            if len(s_clean) >= 12 and all(c in "ACGTN" for c in s_clean):
                candidates.add(s_clean)

        for name in dir(adapters_mod):
            if name.startswith("__"):
                continue
            try:
                obj = getattr(adapters_mod, name)
            except Exception:
                continue
            if isinstance(obj, str):
                add_seq(obj)
            elif isinstance(obj, (list, tuple, set)):
                for item in obj:
                    if isinstance(item, str):
                        add_seq(item)
                    elif isinstance(item, dict):
                        for v in item.values():
                            if isinstance(v, str):
                                add_seq(v)
            elif isinstance(obj, dict):
                for v in obj.values():
                    if isinstance(v, str):
                        add_seq(v)
                    elif isinstance(v, (list, tuple, set)):
                        for x in v:
                            if isinstance(x, str):
                                add_seq(x)
        return candidates

    def _reads_contain_any_sequence(self, reads: str, motifs: Set[str], max_reads: int = 2000) -> bool:
        """Scan the first max_reads for presence of any motif (simple substring search)."""
        try:
            motifs_upper = {m.upper() for m in motifs}
            count = 0
            for record in SeqIO.parse(reads, "fastq"):
                seq_upper = str(record.seq).upper()
                for m in motifs_upper:
                    if m in seq_upper:
                        logger.info("Detected ONT adapter motifs in reads (via Porechop catalog)")
                        return True
                count += 1
                if count >= max_reads:
                    break
        except Exception:
            return False
        return False
    
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
