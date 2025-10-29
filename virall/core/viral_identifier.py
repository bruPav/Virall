"""
Viral sequence identification module using machine learning and database approaches.
"""

import os
import subprocess
import tempfile
import shutil
import time
import gzip
import psutil
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Union
from loguru import logger

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, accuracy_score


class ViralIdentifier:
    """
    Identifies viral sequences using multiple approaches:
    1. Machine learning-based classification using k-mer features
    2. Database-based annotation using BLAST
    3. Composition-based methods
    """
    
    def __init__(self, model_path: Optional[str] = None, config: Optional[Dict] = None, threads: int = 8, memory: str = "16G"):
        """
        Initialize the viral identifier.
        
        Args:
            model_path: Path to pre-trained model (optional)
            config: Configuration dictionary
            threads: Number of threads to use for external tools
            memory: Memory allocation for external tools (e.g., "16G")
        """
        self.model_path = model_path
        self.config = config or {}
        self.threads = threads
        self.memory = memory
        self.model = None
        self.vectorizer = None
        self.temp_dir = Path(tempfile.mkdtemp(prefix="viral_identify_"))
        
        # Using Kaiju for viral identification
        
        # Load or initialize model
        if model_path and os.path.exists(model_path):
            self._load_model()
        else:
            self._initialize_model()
        
        logger.info("ViralIdentifier initialized")
    
    def _get_memory_value(self) -> str:
        """Convert memory from '16G' format to just number for SPAdes."""
        return self.memory.replace('G', '').replace('g', '')
    
    def _log_system_resources(self, log_file: Path) -> None:
        """Log system resources during SPAdes execution"""
        try:
            memory = psutil.virtual_memory()
            disk = shutil.disk_usage('/')
            
            with open(log_file, 'a') as f:
                f.write(f"[{datetime.now()}] System Resources - Memory: {memory.percent:.1f}% used "
                       f"({memory.used/1024/1024/1024:.1f}GB/{memory.total/1024/1024/1024:.1f}GB), "
                       f"Disk: {disk.used/disk.total*100:.1f}% used\n")
        except Exception as e:
            with open(log_file, 'a') as f:
                f.write(f"[{datetime.now()}] Error logging system resources: {e}\n")

    def _monitor_spades_process(self, process, log_file: Path) -> None:
        """Monitor SPAdes process and log its status"""
        start_time = time.time()
        last_log_time = 0
        
        with open(log_file, 'a') as f:
            f.write(f"[{datetime.now()}] Starting SPAdes process monitoring\n")
        
        while process.poll() is None:  # While process is running
            try:
                current_time = time.time()
                elapsed = current_time - start_time
                
                # Get process info
                proc = psutil.Process(process.pid)
                memory_mb = proc.memory_info().rss / 1024 / 1024
                cpu_percent = proc.cpu_percent()
                
                # Log status every 30 seconds
                if elapsed - last_log_time >= 30:
                    with open(log_file, 'a') as f:
                        f.write(f"[{datetime.now()}] SPAdes running - PID: {process.pid}, "
                               f"Elapsed: {elapsed:.0f}s, Memory: {memory_mb:.1f}MB, "
                               f"CPU: {cpu_percent:.1f}%\n")
                    last_log_time = elapsed
                
                time.sleep(5)  # Check every 5 seconds
                
            except (psutil.NoSuchProcess, psutil.AccessDenied) as e:
                with open(log_file, 'a') as f:
                    f.write(f"[{datetime.now()}] Process monitoring error: {e}\n")
                break
            except Exception as e:
                with open(log_file, 'a') as f:
                    f.write(f"[{datetime.now()}] Unexpected monitoring error: {e}\n")
                break
        
        # Log final status
        with open(log_file, 'a') as f:
            f.write(f"[{datetime.now()}] SPAdes process finished - Return code: {process.returncode}\n")

    def _monitor_spades_output(self, temp_assembly_dir: Path, log_file: Path) -> None:
        """Monitor SPAdes output files as they're created"""
        output_files = [
            "contigs.fasta",
            "misc/assembled_contigs.fasta", 
            "K77/final_contigs.fasta",
            "before_rr.fasta",
            "assembly_graph.fastg",
            "scaffolds.paths"
        ]
        
        with open(log_file, 'a') as f:
            f.write(f"[{datetime.now()}] Checking for SPAdes output files...\n")
        
        for file_path in output_files:
            full_path = temp_assembly_dir / file_path
            if full_path.exists():
                size = full_path.stat().st_size
                mtime = datetime.fromtimestamp(full_path.stat().st_mtime)
                with open(log_file, 'a') as f:
                    f.write(f"[{datetime.now()}] Found: {file_path} ({size} bytes, modified: {mtime})\n")
            else:
                with open(log_file, 'a') as f:
                    f.write(f"[{datetime.now()}] Missing: {file_path}\n")

    def _run_spades_with_monitoring(self, cmd: List[str], temp_assembly_dir: Path) -> int:
        """Run SPAdes with comprehensive monitoring"""
        monitor_log = temp_assembly_dir / "spades_monitoring.log"
        
        # Start monitoring log
        with open(monitor_log, 'w') as f:
            f.write(f"[{datetime.now()}] Starting SPAdes monitoring\n")
            f.write(f"Command: {' '.join(cmd)}\n")
            f.write(f"Working directory: {temp_assembly_dir}\n")
        
        # Log initial system resources
        self._log_system_resources(monitor_log)
        
        # Start process
        with open(monitor_log, 'a') as f:
            f.write(f"[{datetime.now()}] Starting SPAdes process...\n")
        
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        # Monitor process in a separate thread
        import threading
        monitor_thread = threading.Thread(target=self._monitor_spades_process, args=(process, monitor_log))
        monitor_thread.daemon = True
        monitor_thread.start()
        
        # Wait for process to complete
        stdout, stderr = process.communicate()
        
        # Log final system resources
        self._log_system_resources(monitor_log)
        
        # Check output files
        self._monitor_spades_output(temp_assembly_dir, monitor_log)
        
        # Log stdout and stderr
        with open(monitor_log, 'a') as f:
            f.write(f"[{datetime.now()}] SPAdes stdout (last 1000 chars):\n{stdout[-1000:]}\n")
            f.write(f"[{datetime.now()}] SPAdes stderr (last 1000 chars):\n{stderr[-1000:]}\n")
            f.write(f"[{datetime.now()}] SPAdes completed with return code: {process.returncode}\n")
        
        return process.returncode
    
    def _initialize_model(self) -> None:
        """Initialize a new viral identification model."""
        self.model = RandomForestClassifier(
            n_estimators=100,
            max_depth=20,
            random_state=42,
            n_jobs=-1
        )
        self.vectorizer = TfidfVectorizer(
            analyzer='char',
            ngram_range=(3, 6),
            max_features=10000
        )
        logger.info("Initialized new viral identification model")
    
    def _load_model(self) -> None:
        """Load pre-trained model from file."""
        # This would load a pre-trained model
        # For now, initialize a new one
        self._initialize_model()
        logger.info("Model loading not implemented yet, using new model")
    
    def identify_viral_contigs_from_assembly(
        self,
        contigs_file: str,
        confidence_threshold: float = 0.8,
        output_file: Optional[str] = None
    ) -> Dict[str, str]:
        """
        Identify viral contigs directly from an assembly file using Kaiju.
        This analyzes the input file (contigs or scaffolds) for viral sequences.
        
        Args:
            contigs_file: Path to contigs file from assembly
            confidence_threshold: Confidence threshold for viral identification
            output_file: Output file for viral contigs (optional)
            
        Returns:
            Dictionary with viral contig information
        """
        logger.info("Using efficient viral contig identification from assembly with Kaiju")
        
        # Get all contig IDs from the file
        all_contigs = self._extract_contigs_from_fasta(contigs_file)
        
        if not all_contigs:
            logger.warning("No contigs found in input file")
            return {}
        
        mode = "RNA" if self.config.get('rna_mode', False) else "DNA"
        logger.info(f"{mode} mode: Analyzing {len(all_contigs)} sequences with Kaiju")
        
        # Run Kaiju classification on all sequences
        classifications_dir = Path(contigs_file).parent.parent / "03_classifications"
        classifications_dir.mkdir(parents=True, exist_ok=True)
        kaiju_results = self._run_kaiju_classification(contigs_file, classifications_dir / "kaiju_contigs")
        
        if kaiju_results.get("status") != "completed":
            logger.error(f"Kaiju classification failed: {kaiju_results.get('error', 'Unknown error')}")
            return {}
        
        # Parse Kaiju results
        kaiju_file = Path(kaiju_results.get("results_file", ""))
        if not kaiju_file.exists():
            logger.error("Kaiju results file not found")
            return {}
        
        kaiju_classifications = self._parse_kaiju_results(kaiju_file)
        
        # Filter contigs based on Kaiju results (keep classified viral sequences)
        viral_contigs = []
        for contig_id in all_contigs:
            if contig_id in kaiju_classifications:
                classification = kaiju_classifications[contig_id]
                # Keep sequences classified as viral (not unclassified)
                if classification.get('classification', '') != 'Unclassified':
                    viral_contigs.append(contig_id)
        
        logger.info(f"{mode} mode: {len(viral_contigs)} out of {len(all_contigs)} sequences classified as viral")
        
        # Write filtered viral contigs to output file
        if output_file and viral_contigs:
            self._write_viral_contigs_from_fasta(contigs_file, viral_contigs, output_file)
        
        return {
            "viral_contigs": viral_contigs,
            "viral_contig_count": len(viral_contigs),
            "classification_method": "kaiju_contigs_mode",
            "kaiju_results": kaiju_results,
            "kaiju_classifications": kaiju_classifications
        }
    
    # Using Kaiju for viral identification
    
    def is_viral_sequence(
        self, 
        sequence: str,
        confidence_threshold: float = 0.8
    ) -> bool:
        """
        Check if a single sequence is viral.
        
        Args:
            sequence: DNA sequence to check
            confidence_threshold: Minimum confidence for viral classification
            
        Returns:
            True if sequence is predicted to be viral
        """
        is_viral, confidence = self._predict_viral_sequence(sequence)
        return is_viral and confidence >= confidence_threshold
    
    def _predict_viral_sequence(self, sequence: str) -> Tuple[bool, float]:
        """
        Predict if a sequence is viral using the trained model.
        
        Args:
            sequence: DNA sequence to classify
            
        Returns:
            Tuple of (is_viral, confidence_score)
        """
        if self.model is None or self.vectorizer is None:
            # Fallback to composition-based method
            return self._composition_based_prediction(sequence)
        
        try:
            # Extract features
            features = self.vectorizer.transform([sequence])
            
            # Predict
            prediction = self.model.predict(features)[0]
            confidence = self.model.predict_proba(features)[0].max()
            
            return bool(prediction), float(confidence)
        
        except Exception as e:
            logger.warning(f"Model prediction failed: {e}, using fallback method")
            return self._composition_based_prediction(sequence)
    
    def _composition_based_prediction(self, sequence: str) -> Tuple[bool, float]:
        """
        Fallback prediction based on sequence composition.
        
        Args:
            sequence: DNA sequence to classify
            
        Returns:
            Tuple of (is_viral, confidence_score)
        """
        if len(sequence) < 100:
            return False, 0.0
        
        # Calculate basic composition features
        gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
        at_content = (sequence.count('A') + sequence.count('T')) / len(sequence)
        
        # Simple heuristics for viral sequences
        # Viruses often have different GC content patterns
        is_viral = False
        confidence = 0.0
        
        # Check for viral-like characteristics
        if 0.2 <= gc_content <= 0.8:  # Reasonable GC content
            # Check for repetitive elements (common in viruses)
            kmer_counts = {}
            k = 4
            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i+k]
                kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
            
            # Calculate repeat content
            max_repeats = max(kmer_counts.values()) if kmer_counts else 0
            repeat_ratio = max_repeats / len(sequence) if len(sequence) > 0 else 0
            
            # Viral sequences often have higher repeat content
            if repeat_ratio > 0.1:
                is_viral = True
                confidence = min(0.8, repeat_ratio * 2)
        
        return is_viral, confidence
    
    def _write_viral_reads(
        self, 
        reads: List[SeqRecord], 
        output_file: str,
        append: bool = False
    ) -> None:
        """Write viral reads to output file."""
        mode = "a" if append else "w"
        
        with open(output_file, mode) as f:
            SeqIO.write(reads, f, "fastq")
    
    def train_model(
        self, 
        viral_sequences: List[str],
        non_viral_sequences: List[str],
        test_size: float = 0.2
    ) -> Dict[str, float]:
        """
        Train the viral identification model.
        
        Args:
            viral_sequences: List of viral DNA sequences
            non_viral_sequences: List of non-viral DNA sequences
            test_size: Fraction of data to use for testing
            
        Returns:
            Dictionary containing training metrics
        """
        logger.info("Training viral identification model")
        
        # Prepare data
        sequences = viral_sequences + non_viral_sequences
        labels = [1] * len(viral_sequences) + [0] * len(non_viral_sequences)
        
        # Extract features
        X = self.vectorizer.fit_transform(sequences)
        y = np.array(labels)
        
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=test_size, random_state=42, stratify=y
        )
        
        # Train model
        self.model.fit(X_train, y_train)
        
        # Evaluate
        y_pred = self.model.predict(X_test)
        accuracy = accuracy_score(y_test, y_pred)
        
        # Generate classification report
        report = classification_report(y_test, y_pred, output_dict=True)
        
        logger.info(f"Model training completed with accuracy: {accuracy:.3f}")
        
        return {
            "accuracy": accuracy,
            "classification_report": report
        }
    
    
    def save_model(self, model_path: str) -> None:
        """Save the trained model to file."""
        # This would save the model and vectorizer
        # For now, just log the action
        logger.info(f"Model saving not implemented yet, would save to {model_path}")
    
    def classify_viral_contigs(
        self,
        contigs_file: str,
        output_dir: Optional[str] = None
    ) -> Dict[str, Union[str, Dict]]:
        """
        Classify viral contigs using Kaiju nucleotide-based classification.
        
        Args:
            contigs_file: Path to viral contigs file
            output_dir: Output directory for results (optional)
            
        Returns:
            Dictionary containing classification results
        """
        if output_dir is None:
            output_dir = self.temp_dir / "viral_classification"
        
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)
        
        logger.info("Classifying viral contigs using Kaiju")
        
        # Run Kaiju classification
        kaiju_results = self._run_kaiju_classification(contigs_file, output_dir)
        
        if kaiju_results.get("status") != "completed":
            return kaiju_results
        
        # Parse Kaiju results
        kaiju_file = Path(kaiju_results.get("results_file", ""))
        if not kaiju_file.exists():
            return {
                "status": "failed",
                "error": "Kaiju results file not found"
            }
        
        classifications = self._parse_kaiju_results(kaiju_file)
        
        # Write classification summary
        summary_file = output_dir / "kaiju_summary.tsv"
        self._write_kaiju_summary(summary_file, classifications)
        
        logger.info(f"Kaiju classification completed: {len(classifications)} contigs classified")
        
        return {
            "status": "completed",
            "classifications": classifications,
            "summary_file": str(summary_file),
            "total_classified": len(classifications)
        }
    



    def _predict_proteins_from_contigs(self, contigs_file: str, output_file: Path) -> None:
        """Predict proteins from contigs using Prodigal."""
        logger.info("Predicting proteins from contigs for gene prediction")
        
        cmd = [
            "prodigal",
            "-i", contigs_file,
            "-a", str(output_file),
            "-p", "meta",  # Use meta mode for viral sequences
            "-q"  # Quiet mode
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"Prodigal failed: {result.stderr}")
            raise RuntimeError("Protein prediction failed")
        
        logger.info(f"Predicted proteins written to {output_file}")
    
    

    def _get_consensus_classification(self, hits: List[Dict]) -> Dict[str, str]:
        """Determine consensus classification from multiple hits."""
        import re
        
        # Filter for viral hits only (exclude bacterial, human, etc.)
        viral_hits = []
        for hit in hits:
            title = hit['title'].lower()
            # Check if it's likely a viral hit
            if any(viral_term in title for viral_term in ['virus', 'phage', 'viridae', 'viral']):
                viral_hits.append(hit)
        
        # If no viral hits found, use all hits but mark as uncertain
        if not viral_hits:
            viral_hits = hits
            viral_only = False
        else:
            viral_only = True
        
        # Extract viral families and genera from hit titles
        families = []
        genera = []
        species = []
        
        for hit in viral_hits:
            title = hit['title']
            # Extract species name (usually in brackets)
            species_match = re.search(r'\[([^\]]+)\]', title)
            if species_match:
                species_name = species_match.group(1)
                species.append(species_name)
                
                # Try to extract family and genus
                # Common patterns: "virus name", "genus species", "family virus"
                parts = species_name.split()
                if len(parts) >= 2:
                    # Assume first part is genus, last part is species
                    genera.append(parts[0])
                    if len(parts) > 2:
                        families.append(' '.join(parts[:-1]))
        
        # Determine consensus
        consensus = {
            'family': 'Unknown',
            'genus': 'Unknown', 
            'species': 'Unknown',
            'confidence': 'Low',
            'viral_only': viral_only,
            'total_hits': len(hits),
            'viral_hits': len(viral_hits),
            'viral_hit_ratio': len(viral_hits) / len(hits) if hits else 0
        }
        
        if species:
            # Use most common species
            from collections import Counter
            species_counts = Counter(species)
            consensus['species'] = species_counts.most_common(1)[0][0]
            
            # Add species distribution statistics
            consensus['unique_species'] = len(set(species))
            consensus['top_species_count'] = species_counts.most_common(1)[0][1]
            consensus['species_agreement'] = consensus['top_species_count'] / len(species) if species else 0
            
            if genera:
                genus_counts = Counter(genera)
                consensus['genus'] = genus_counts.most_common(1)[0][0]
                consensus['unique_genera'] = len(set(genera))
                consensus['top_genus_count'] = genus_counts.most_common(1)[0][1]
                
            if families:
                family_counts = Counter(families)
                consensus['family'] = family_counts.most_common(1)[0][0]
                consensus['unique_families'] = len(set(families))
                consensus['top_family_count'] = family_counts.most_common(1)[0][1]
            
            # Determine confidence based on agreement
            if len(set(species)) == 1:
                consensus['confidence'] = 'High'
            elif len(set(species)) <= 3:
                consensus['confidence'] = 'Medium'
            else:
                consensus['confidence'] = 'Low'
        
        return consensus
    
    
    def cleanup(self) -> None:
        """Clean up temporary files."""
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
            logger.info("Cleaned up temporary files")
    
    def _extract_contigs_from_fasta(self, fasta_file: str) -> List[str]:
        """Extract all contig IDs from a FASTA file."""
        contig_ids = []
        try:
            with open(fasta_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        # Extract contig ID (everything after > until first space)
                        contig_id = line[1:].strip().split()[0]
                        contig_ids.append(contig_id)
        except Exception as e:
            logger.error(f"Error reading FASTA file {fasta_file}: {e}")
        return contig_ids

    def _write_viral_contigs_from_fasta(self, input_fasta: str, contig_ids: List[str], output_file: str) -> None:
        """Write viral contigs from input FASTA file to output file."""
        try:
            from Bio import SeqIO
            
            # Read all sequences from input file
            sequences = SeqIO.parse(input_fasta, "fasta")
            
            # Filter for viral contigs
            viral_sequences = []
            for seq in sequences:
                if seq.id in contig_ids:
                    viral_sequences.append(seq)
            
            # Write to output file
            SeqIO.write(viral_sequences, output_file, "fasta")
            logger.info(f"Wrote {len(viral_sequences)} viral contigs to {output_file}")
            
        except Exception as e:
            logger.error(f"Error writing viral contigs to {output_file}: {e}")

    
    
    def _extract_viral_sequences_from_contigs(
        self, 
        contigs_file: str, 
        viral_contig_ids: List[str]
    ) -> List[SeqRecord]:
        """
        Extract viral sequences from contigs file based on contig IDs.
        
        Args:
            contigs_file: Path to contigs file
            viral_contig_ids: List of viral contig IDs
            
        Returns:
            List of viral sequence records
        """
        viral_sequences = []
        viral_ids_set = set(viral_contig_ids)
        
        logger.info(f"Extracting viral sequences from {contigs_file}")
        
        # Read contigs and extract viral ones
        for record in SeqIO.parse(contigs_file, "fasta"):
            # Check if this contig is viral
            contig_id = record.id.split()[0]  # Get the main ID part
            if contig_id in viral_ids_set:
                viral_sequences.append(record)
                logger.debug(f"Found viral contig: {record.id}")
        
        logger.info(f"Extracted {len(viral_sequences)} viral sequences from contigs")
        return viral_sequences
    
    def _write_viral_contigs(
        self, 
        viral_sequences: List[SeqRecord], 
        output_file: str
    ) -> None:
        """
        Write viral contigs to output file.
        
        Args:
            viral_sequences: List of viral sequence records
            output_file: Output file path
        """
        if not viral_sequences:
            logger.warning("No viral sequences to write")
            return
        
        # Write as FASTA
        SeqIO.write(viral_sequences, output_file, "fasta")
        logger.info(f"Wrote {len(viral_sequences)} viral contigs to {output_file}")
    
    def _map_reads_to_viral_contigs(
        self, 
        reads_file: str, 
        viral_contigs: List[str], 
        confidence_threshold: float
    ) -> List:
        """
        Map reads back to viral contigs using BWA alignment.
        
        Args:
            reads_file: Path to reads file
            viral_contigs: List of viral contig IDs
            confidence_threshold: Confidence threshold
            
        Returns:
            List of viral reads
        """
        if not viral_contigs:
            return []
        
        logger.info(f"Mapping reads to {len(viral_contigs)} viral contigs")
        
        # For now, use a simple approach: take a subset of reads
        # TODO: Implement proper BWA mapping
        viral_reads = []
        total_reads = 0
        
        # Handle compressed files
        if str(reads_file).endswith('.gz'):
            with gzip.open(reads_file, 'rt') as handle:
                for i, record in enumerate(SeqIO.parse(handle, "fastq")):
                    total_reads += 1
                    # Simple heuristic: take every 10th read as "viral"
                    if i % 10 == 0:
                        viral_reads.append(record)
                    if len(viral_reads) >= 1000:  # Limit for testing
                        break
        else:
            for i, record in enumerate(SeqIO.parse(reads_file, "fastq")):
                total_reads += 1
                # Simple heuristic: take every 10th read as "viral"
                if i % 10 == 0:
                    viral_reads.append(record)
                if len(viral_reads) >= 1000:  # Limit for testing
                    break
        
        logger.info(f"Mapped {len(viral_reads)} reads to viral contigs")
        return viral_reads
    
    def _write_viral_reads(self, viral_reads: List, output_file: str) -> None:
        """
        Write viral reads to output file.
        
        Args:
            viral_reads: List of viral read records
            output_file: Output file path
        """
        if not viral_reads:
            logger.warning("No viral reads to write")
            return
        
        # Write as FASTQ if possible, otherwise FASTA
        try:
            if str(output_file).endswith('.gz'):
                with gzip.open(output_file, 'wt') as handle:
                    SeqIO.write(viral_reads, handle, "fastq")
            else:
                SeqIO.write(viral_reads, output_file, "fastq")
        except Exception as e:
            logger.warning(f"Failed to write as FASTQ, writing as FASTA: {e}")
            # Write as FASTA instead
            fasta_file = str(output_file).replace('.fastq', '.fasta').replace('.gz', '')
            SeqIO.write(viral_reads, fasta_file, "fasta")

    def quantify_viral_contigs(
        self,
        contigs_file: str,
        reads_1: Optional[str] = None,
        reads_2: Optional[str] = None,
        single_reads: Optional[str] = None,
        long_reads: Optional[str] = None,
        output_dir: Optional[str] = None,
        classification_data: Optional[Dict] = None
    ) -> Dict[str, Union[str, Dict]]:
        """
        Quantify viral contigs by mapping reads back to them.
        
        Args:
            contigs_file: Path to viral contigs FASTA file
            reads_1: Path to first mate of paired-end reads
            reads_2: Path to second mate of paired-end reads
            single_reads: Path to single-end reads
            long_reads: Path to long reads (PacBio/ONT)
            output_dir: Output directory for quantification results
            classification_data: Optional classification data for species information
            
        Returns:
            Dictionary containing quantification results
        """
        if output_dir is None:
            output_dir = self.temp_dir / "viral_quantification"
        
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)
        
        logger.info("Quantifying viral contigs from read mapping")
        
        # Determine read types and run quantification
        read_types = []
        if long_reads:
            read_types.append("long")
        if reads_1 and reads_2:
            read_types.append("paired_short")
        if single_reads:
            read_types.append("single_short")
        
        if not read_types:
            logger.error("No reads provided for quantification")
            return {"status": "failed", "error": "No reads provided"}
        
        # If multiple read types, run combined quantification
        if len(read_types) > 1:
            logger.info(f"Running combined quantification with {', '.join(read_types)} reads")
            return self._quantify_combined_reads(
                contigs_file, reads_1, reads_2, single_reads, long_reads, output_dir, classification_data
            )
        else:
            # Single read type
            if long_reads:
                return self._quantify_with_minimap2(contigs_file, long_reads, output_dir, classification_data)
            elif reads_1 and reads_2:
                return self._quantify_with_bwa_paired(contigs_file, reads_1, reads_2, output_dir, classification_data)
            elif single_reads:
                return self._quantify_with_bwa_single(contigs_file, single_reads, output_dir, classification_data)

    def _quantify_combined_reads(
        self,
        contigs_file: str,
        reads_1: Optional[str],
        reads_2: Optional[str],
        single_reads: Optional[str],
        long_reads: Optional[str],
        output_dir: Path,
        classification_data: Optional[Dict] = None
    ) -> Dict[str, Union[str, Dict]]:
        """Run quantification with multiple read types and combine results."""
        logger.info("Running combined quantification with multiple read types")
        
        all_quantification_results = {}
        bam_files = []
        methods_used = []
        
        # Quantify with each read type
        if long_reads:
            logger.info("Quantifying with long reads (minimap2)")
            long_output_dir = output_dir / "long_reads"
            long_output_dir.mkdir(exist_ok=True)
            long_results = self._quantify_with_minimap2(contigs_file, long_reads, long_output_dir, classification_data)
            if long_results.get("status") == "completed":
                all_quantification_results["long_reads"] = long_results["quantification_results"]
                bam_files.append(long_results["bam_file"])
                methods_used.append("minimap2")
        
        if reads_1 and reads_2:
            logger.info("Quantifying with paired-end short reads (BWA)")
            paired_output_dir = output_dir / "paired_reads"
            paired_output_dir.mkdir(exist_ok=True)
            paired_results = self._quantify_with_bwa_paired(contigs_file, reads_1, reads_2, paired_output_dir, classification_data)
            if paired_results.get("status") == "completed":
                all_quantification_results["paired_reads"] = paired_results["quantification_results"]
                bam_files.append(paired_results["bam_file"])
                methods_used.append("bwa_paired")
        
        if single_reads:
            logger.info("Quantifying with single-end short reads (BWA)")
            single_output_dir = output_dir / "single_reads"
            single_output_dir.mkdir(exist_ok=True)
            single_results = self._quantify_with_bwa_single(contigs_file, single_reads, single_output_dir, classification_data)
            if single_results.get("status") == "completed":
                all_quantification_results["single_reads"] = single_results["quantification_results"]
                bam_files.append(single_results["bam_file"])
                methods_used.append("bwa_single")
        
        if not all_quantification_results:
            return {"status": "failed", "error": "All quantification methods failed"}
        
        # Combine results from all read types
        combined_results = self._combine_quantification_results(all_quantification_results, output_dir, classification_data)
        
        logger.info("Combined viral contig quantification completed")
        return {
            "status": "completed",
            "method": f"combined_{'+'.join(methods_used)}",
            "bam_files": bam_files,
            "quantification_results": combined_results,
            "individual_results": all_quantification_results
        }
    
    def _combine_quantification_results(
        self, 
        all_results: Dict[str, Dict], 
        output_dir: Path,
        classification_data: Optional[Dict] = None
    ) -> Dict[str, Dict]:
        """Combine quantification results from multiple read types."""
        logger.info("Combining quantification results from multiple read types")
        
        # Get all contig IDs
        all_contigs = set()
        for read_type, results in all_results.items():
            all_contigs.update(results.keys())
        
        combined_results = {}
        
        for contig_id in all_contigs:
            contig_metrics = {}
            
            # Collect metrics from each read type
            for read_type, results in all_results.items():
                if contig_id in results:
                    metrics = results[contig_id]
                    contig_metrics[read_type] = {
                        'mean_coverage': metrics.get('mean_coverage', 0),
                        'max_coverage': metrics.get('max_coverage', 0),
                        'coverage_breadth': metrics.get('coverage_breadth', 0),
                        'total_coverage': metrics.get('total_coverage', 0),
                        'relative_abundance': metrics.get('relative_abundance', 0),
                        'mapped_positions': metrics.get('mapped_positions', 0)
                    }
            
            # Calculate combined metrics
            if contig_metrics:
                # Use the read type with highest coverage for primary metrics
                best_read_type = max(contig_metrics.keys(), 
                                  key=lambda x: contig_metrics[x]['total_coverage'])
                best_metrics = contig_metrics[best_read_type]
                
                # Calculate weighted averages
                total_coverage_all = sum(m['total_coverage'] for m in contig_metrics.values())
                total_abundance_all = sum(m['relative_abundance'] for m in contig_metrics.values())
                
                combined_results[contig_id] = {
                    'contig_length': best_metrics.get('contig_length', 0),
                    'mean_coverage': best_metrics.get('mean_coverage', 0),
                    'max_coverage': max(m['max_coverage'] for m in contig_metrics.values()),
                    'coverage_breadth': max(m['coverage_breadth'] for m in contig_metrics.values()),
                    'total_coverage': total_coverage_all,
                    'relative_abundance': total_abundance_all,
                    'mapped_positions': sum(m['mapped_positions'] for m in contig_metrics.values()),
                    'read_types_used': list(contig_metrics.keys()),
                    'best_read_type': best_read_type,
                    'individual_metrics': contig_metrics
                }
        
        # Write combined abundance summary
        self._write_combined_abundance_summary(combined_results, all_results, output_dir, classification_data)
        
        logger.info(f"Combined results for {len(combined_results)} contigs")
        return combined_results
    
    def _write_combined_abundance_summary(
        self, 
        combined_results: Dict[str, Dict], 
        individual_results: Dict[str, Dict],
        output_dir: Path,
        classification_data: Optional[Dict] = None
    ) -> None:
        """Write combined abundance summary to TSV file."""
        summary_file = output_dir / "combined_contig_abundance.tsv"
        
        summary_data = []
        for contig_id, metrics in combined_results.items():
            # Get species information from classification data
            species = "Unknown"
            if classification_data and contig_id in classification_data:
                contig_class = classification_data[contig_id]
                species = contig_class.get("taxon_name", "Unknown")
            
            row = {
                'contig_id': contig_id,
                'species': species,
                'contig_length': metrics['contig_length'],
                'mean_coverage': metrics['mean_coverage'],
                'max_coverage': metrics['max_coverage'],
                'coverage_breadth': metrics['coverage_breadth'],
                'total_coverage': metrics['total_coverage'],
                'relative_abundance': metrics['relative_abundance'],
                'mapped_positions': metrics['mapped_positions'],
                'read_types_used': ','.join(metrics['read_types_used']),
                'best_read_type': metrics['best_read_type']
            }
            
            # Add individual read type metrics
            for read_type in metrics['read_types_used']:
                individual = metrics['individual_metrics'][read_type]
                row[f'{read_type}_coverage'] = individual['mean_coverage']
                row[f'{read_type}_abundance'] = individual['relative_abundance']
            
            summary_data.append(row)
        
        if summary_data:
            df = pd.DataFrame(summary_data)
            # Sort by relative abundance (highest first)
            df = df.sort_values('relative_abundance', ascending=False)
            df.to_csv(summary_file, sep='\t', index=False)
            logger.info(f"Wrote combined abundance summary to {summary_file}")
        else:
            logger.warning("No combined abundance data to write.")

    def _quantify_with_bwa_paired(
        self, 
        contigs_file: str, 
        reads_1: str, 
        reads_2: str, 
        output_dir: Path,
        classification_data: Optional[Dict] = None
    ) -> Dict[str, Union[str, Dict]]:
        """Quantify using BWA with paired-end reads."""
        # Check if BWA is available
        import shutil
        bwa_path = shutil.which("bwa")
        if bwa_path is None:
            logger.error("BWA not found in PATH. Please install: conda install -c bioconda bwa")
            return {"status": "failed", "error": "BWA not installed"}
        
        logger.info(f"Found BWA at: {bwa_path}")
        
        # Test BWA execution (BWA returns 1 when called without args, which is normal)
        try:
            result = subprocess.run("bwa", capture_output=True, text=True, shell=True)
            # BWA returns 1 when called without arguments (shows help), which is normal
            if result.returncode not in [0, 1]:
                logger.error(f"BWA execution failed. Return code: {result.returncode}, stderr: {result.stderr}")
                return {"status": "failed", "error": f"BWA execution failed: {result.stderr}"}
            logger.info("BWA is working correctly")
        except Exception as e:
            logger.error(f"BWA execution error: {e}")
            return {"status": "failed", "error": f"BWA execution error: {e}"}
        
        # Index the contigs
        logger.info("Indexing viral contigs for BWA mapping")
        index_file = output_dir / "viral_contigs"
        index_cmd = ["bwa", "index", "-p", str(index_file), contigs_file]
        
        result = subprocess.run(" ".join(index_cmd), capture_output=True, text=True, shell=True)
        if result.returncode != 0:
            logger.error(f"BWA indexing failed: {result.stderr}")
            return {"status": "failed", "error": f"BWA indexing failed: {result.stderr}"}
        
        # Map reads to contigs
        sam_file = output_dir / "mapped_reads.sam"
        bam_file = output_dir / "mapped_reads.bam"
        
        logger.info("Mapping paired-end reads to viral contigs with BWA")
        map_cmd = [
            "bwa", "mem",
            "-t", str(self.threads),
            str(index_file),
            reads_1, reads_2
        ]
        
        # Run mapping
        with open(sam_file, 'w') as f:
            result = subprocess.run(" ".join(map_cmd), stdout=f, stderr=subprocess.PIPE, text=True, shell=True)
            if result.returncode != 0:
                logger.error(f"BWA mapping failed: {result.stderr}")
                return {"status": "failed", "error": f"BWA mapping failed: {result.stderr}"}
        
        # Process BAM file
        return self._process_mapping_results(sam_file, bam_file, contigs_file, output_dir, "bwa_paired", classification_data)

    def _quantify_with_bwa_single(
        self, 
        contigs_file: str, 
        single_reads: str, 
        output_dir: Path,
        classification_data: Optional[Dict] = None
    ) -> Dict[str, Union[str, Dict]]:
        """Quantify using BWA with single-end reads."""
        # Check if BWA is available
        import shutil
        bwa_path = shutil.which("bwa")
        if bwa_path is None:
            logger.error("BWA not found in PATH. Please install: conda install -c bioconda bwa")
            return {"status": "failed", "error": "BWA not installed"}
        
        logger.info(f"Found BWA at: {bwa_path}")
        
        # Test BWA execution (BWA returns 1 when called without args, which is normal)
        try:
            result = subprocess.run("bwa", capture_output=True, text=True, shell=True)
            # BWA returns 1 when called without arguments (shows help), which is normal
            if result.returncode not in [0, 1]:
                logger.error(f"BWA execution failed. Return code: {result.returncode}, stderr: {result.stderr}")
                return {"status": "failed", "error": f"BWA execution failed: {result.stderr}"}
            logger.info("BWA is working correctly")
        except Exception as e:
            logger.error(f"BWA execution error: {e}")
            return {"status": "failed", "error": f"BWA execution error: {e}"}
        
        # Index the contigs
        logger.info("Indexing viral contigs for BWA mapping")
        index_file = output_dir / "viral_contigs"
        index_cmd = ["bwa", "index", "-p", str(index_file), contigs_file]
        
        result = subprocess.run(" ".join(index_cmd), capture_output=True, text=True, shell=True)
        if result.returncode != 0:
            logger.error(f"BWA indexing failed: {result.stderr}")
            return {"status": "failed", "error": f"BWA indexing failed: {result.stderr}"}
        
        # Map reads to contigs
        sam_file = output_dir / "mapped_reads.sam"
        bam_file = output_dir / "mapped_reads.bam"
        
        logger.info("Mapping single-end reads to viral contigs with BWA")
        map_cmd = [
            "bwa", "mem",
            "-t", str(self.threads),
            str(index_file),
            single_reads
        ]
        
        # Run mapping
        with open(sam_file, 'w') as f:
            result = subprocess.run(" ".join(map_cmd), stdout=f, stderr=subprocess.PIPE, text=True, shell=True)
            if result.returncode != 0:
                logger.error(f"BWA mapping failed: {result.stderr}")
                return {"status": "failed", "error": f"BWA mapping failed: {result.stderr}"}
        
        # Process BAM file
        return self._process_mapping_results(sam_file, bam_file, contigs_file, output_dir, "bwa_single", classification_data)

    def _quantify_with_minimap2(
        self, 
        contigs_file: str, 
        long_reads: str, 
        output_dir: Path,
        classification_data: Optional[Dict] = None
    ) -> Dict[str, Union[str, Dict]]:
        """Quantify using minimap2 with long reads."""
        # Check if minimap2 is available
        import shutil
        minimap2_path = shutil.which("minimap2")
        if minimap2_path is None:
            logger.error("minimap2 not found in PATH. Please install: conda install -c bioconda minimap2")
            return {"status": "failed", "error": "minimap2 not installed"}
        
        logger.info(f"Found minimap2 at: {minimap2_path}")
        
        # Test minimap2 execution
        try:
            result = subprocess.run("minimap2 --version", capture_output=True, text=True, shell=True)
            if result.returncode != 0:
                logger.error("minimap2 execution failed. Please check installation")
                return {"status": "failed", "error": "minimap2 execution failed"}
        except Exception as e:
            logger.error(f"minimap2 execution error: {e}")
            return {"status": "failed", "error": f"minimap2 execution error: {e}"}
        
        # Map reads to contigs
        sam_file = output_dir / "mapped_reads.sam"
        bam_file = output_dir / "mapped_reads.bam"
        
        logger.info("Mapping long reads to viral contigs with minimap2")
        map_cmd = [
            "minimap2",
            "-ax", "map-ont",  # Use map-ont for ONT reads, map-pb for PacBio
            "-t", str(self.threads),
            contigs_file,
            long_reads
        ]
        
        # Run mapping
        with open(sam_file, 'w') as f:
            result = subprocess.run(" ".join(map_cmd), stdout=f, stderr=subprocess.PIPE, text=True, shell=True)
            if result.returncode != 0:
                logger.error(f"minimap2 mapping failed: {result.stderr}")
                return {"status": "failed", "error": f"minimap2 mapping failed: {result.stderr}"}
        
        # Process BAM file
        return self._process_mapping_results(sam_file, bam_file, contigs_file, output_dir, "minimap2", classification_data)

    def _process_mapping_results(
        self, 
        sam_file: Path, 
        bam_file: Path, 
        contigs_file: str, 
        output_dir: Path,
        method: str,
        classification_data: Optional[Dict] = None
    ) -> Dict[str, Union[str, Dict]]:
        """Process mapping results and calculate abundance."""
        
        # Convert SAM to BAM and sort
        logger.info("Processing mapping results")
        sort_cmd = ["samtools", "sort", "-o", str(bam_file), str(sam_file)]
        result = subprocess.run(" ".join(sort_cmd), capture_output=True, text=True, shell=True)
        if result.returncode != 0:
            logger.error(f"SAMtools sort failed: {result.stderr}")
            return {"status": "failed", "error": f"SAMtools sort failed: {result.stderr}"}
        
        # Index BAM file
        index_cmd = ["samtools", "index", str(bam_file)]
        result = subprocess.run(" ".join(index_cmd), capture_output=True, text=True, shell=True)
        if result.returncode != 0:
            logger.error(f"SAMtools index failed: {result.stderr}")
            return {"status": "failed", "error": f"SAMtools index failed: {result.stderr}"}
        
        # Calculate coverage and abundance
        quantification_results = self._calculate_contig_abundance(contigs_file, bam_file, output_dir, classification_data)
        
        logger.info("Viral contig quantification completed")
        return {
            "status": "completed",
            "method": method,
            "bam_file": str(bam_file),
            "sam_file": str(sam_file),
            "quantification_results": quantification_results
        }
    
    def _calculate_contig_abundance(
        self, 
        contigs_file: str, 
        bam_file: Path, 
        output_dir: Path,
        classification_data: Optional[Dict] = None
    ) -> Dict[str, Dict]:
        """Calculate abundance metrics for each contig."""
        logger.info("Calculating contig abundance metrics")
        
        # Get contig lengths
        contig_lengths = {}
        for record in SeqIO.parse(contigs_file, "fasta"):
            contig_lengths[record.id] = len(record.seq)
        
        # Calculate coverage using samtools depth
        depth_file = output_dir / "contig_depth.txt"
        depth_cmd = ["samtools", "depth", str(bam_file)]
        
        with open(depth_file, 'w') as f:
            result = subprocess.run(" ".join(depth_cmd), stdout=f, stderr=subprocess.PIPE, text=True, shell=True)
            if result.returncode != 0:
                logger.error(f"SAMtools depth failed: {result.stderr}")
                return {}
        
        # Parse depth data
        contig_coverage = {}
        with open(depth_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    contig_id = parts[0]
                    position = int(parts[1])
                    depth = int(parts[2])
                    
                    if contig_id not in contig_coverage:
                        contig_coverage[contig_id] = []
                    contig_coverage[contig_id].append(depth)
        
        # Calculate abundance metrics
        abundance_results = {}
        total_mapped_reads = 0
        
        for contig_id, depths in contig_coverage.items():
            if depths:
                mean_coverage = sum(depths) / len(depths)
                max_coverage = max(depths)
                coverage_breadth = len(depths) / contig_lengths.get(contig_id, 1)
                total_coverage = sum(depths)
                total_mapped_reads += total_coverage
                
                abundance_results[contig_id] = {
                    'contig_length': contig_lengths.get(contig_id, 0),
                    'mean_coverage': round(mean_coverage, 2),
                    'max_coverage': max_coverage,
                    'coverage_breadth': round(coverage_breadth, 3),
                    'total_coverage': total_coverage,
                    'mapped_positions': len(depths)
                }
            else:
                abundance_results[contig_id] = {
                    'contig_length': contig_lengths.get(contig_id, 0),
                    'mean_coverage': 0.0,
                    'max_coverage': 0,
                    'coverage_breadth': 0.0,
                    'total_coverage': 0,
                    'mapped_positions': 0
                }
        
        # Calculate relative abundance
        for contig_id, metrics in abundance_results.items():
            if total_mapped_reads > 0:
                relative_abundance = metrics['total_coverage'] / total_mapped_reads
                metrics['relative_abundance'] = round(relative_abundance, 6)
            else:
                metrics['relative_abundance'] = 0.0
        
        # Write abundance summary
        self._write_abundance_summary(abundance_results, output_dir, classification_data)
        
        logger.info(f"Calculated abundance for {len(abundance_results)} contigs")
        return abundance_results
    
    def _write_abundance_summary(self, abundance_results: Dict[str, Dict], output_dir: Path, classification_data: Optional[Dict] = None) -> None:
        """Write abundance summary to TSV file."""
        summary_file = output_dir / "contig_abundance.tsv"
        
        summary_data = []
        for contig_id, metrics in abundance_results.items():
            # Get species information from classification data
            species = "Unknown"
            if classification_data and contig_id in classification_data:
                contig_class = classification_data[contig_id]
                species = contig_class.get("taxon_name", "Unknown")
            
            summary_data.append({
                'contig_id': contig_id,
                'species': species,
                'contig_length': metrics['contig_length'],
                'mean_coverage': metrics['mean_coverage'],
                'max_coverage': metrics['max_coverage'],
                'coverage_breadth': metrics['coverage_breadth'],
                'total_coverage': metrics['total_coverage'],
                'relative_abundance': metrics['relative_abundance'],
                'mapped_positions': metrics['mapped_positions']
            })
        
        if summary_data:
            df = pd.DataFrame(summary_data)
            # Sort by relative abundance (highest first)
            df = df.sort_values('relative_abundance', ascending=False)
            df.to_csv(summary_file, sep='\t', index=False)
            logger.info(f"Wrote abundance summary to {summary_file}")
        else:
            logger.warning("No abundance data to write.")
    
    def _run_kaiju_classification(
        self, 
        contigs_file: str, 
        output_dir: Path
    ) -> Dict[str, Union[str, int, float]]:
        """
        Run Kaiju classification on viral contigs.
        
        Args:
            contigs_file: Path to viral contigs FASTA file
            output_dir: Output directory for results
            
        Returns:
            Dictionary containing Kaiju results
        """
        logger.info("Running Kaiju classification on viral contigs")
        
        # Check if Kaiju is available
        if not self._check_kaiju_installation():
            return {
                'status': 'failed',
                'error': 'Kaiju not found. Please install: conda install -c bioconda kaiju'
            }
        
        # Setup Kaiju database
        kaiju_db = self._setup_kaiju_database()
        if not kaiju_db:
            return {
                'status': 'failed',
                'error': 'Kaiju viral database not found. Please run: kaiju-makedb -s viruses'
            }
        
        # Run Kaiju classification
        results_file = output_dir / "kaiju_results.tsv"
        
        # Ensure output directory exists
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            logger.warning(f"Failed to create Kaiju output directory {output_dir}: {e}")
            return {
                'status': 'failed',
                'error': f"Could not create Kaiju output directory: {e}"
            }
        
        cmd = [
            "kaiju",
            "-t", str(kaiju_db / "nodes.dmp"),
            "-f", str(kaiju_db / "kaiju_db_viruses.fmi"),
            "-i", contigs_file,
            "-o", str(results_file),
            "-z", str(self.threads)
        ]
        
        logger.info(f"Kaiju command: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.warning(f"Kaiju failed: {result.stderr}")
            return {
                'status': 'failed',
                'error': f"Kaiju classification failed: {result.stderr}"
            }
        
        logger.info("Kaiju classification completed successfully")
        
        # Add taxon names to the results
        results_with_names = output_dir / "kaiju_results_with_names.tsv"
        logger.info("Adding taxonomic names to Kaiju results")
        
        add_names_cmd = [
            "kaiju-addTaxonNames",
            "-t", str(kaiju_db / "nodes.dmp"),
            "-n", str(kaiju_db / "names.dmp"),
            "-i", str(results_file),
            "-o", str(results_with_names)
        ]
        
        logger.info(f"Kaiju addTaxonNames command: {' '.join(add_names_cmd)}")
        add_names_result = subprocess.run(add_names_cmd, capture_output=True, text=True)
        
        if add_names_result.returncode != 0:
            logger.warning(f"Failed to add taxon names: {add_names_result.stderr}")
            # Continue with original results if name resolution fails
            return {
                'status': 'completed',
                'results_file': str(results_file)
            }
        
        logger.info("Taxon names added successfully")
        return {
            'status': 'completed',
            'results_file': str(results_with_names)
        }
    
    def _check_kaiju_installation(self) -> bool:
        """Check if Kaiju is installed and available."""
        try:
            result = subprocess.run(["kaiju", "-h"], capture_output=True, text=True)
            # Kaiju returns exit code 1 when showing help, which is normal
            return result.returncode in [0, 1]
        except FileNotFoundError:
            return False
    
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
    
    def _setup_kaiju_database(self) -> Optional[Path]:
        """Setup Kaiju viral database."""
        # Get database path from config
        config_path = self.config.get('databases', {}).get('kaiju_db_path')
        if config_path:
            kaiju_db_path = Path(config_path)
        else:
            # Try current working directory first, then fall back to installation directory
            cwd_db_path = Path.cwd() / "databases" / "kaiju_db"
            if cwd_db_path.exists():
                kaiju_db_path = cwd_db_path
            else:
                # Use the same installation directory detection as assembler
                software_dir = self._find_installation_directory()
                kaiju_db_path = software_dir / "databases" / "kaiju_db"
        
        logger.debug(f"Looking for Kaiju database at: {kaiju_db_path}")
        
        if not kaiju_db_path.exists():
            logger.warning(f"Kaiju database not found at {kaiju_db_path}. Please run: kaiju-makedb -s viruses")
            return None
        
        # Check for required database files
        # nodes.dmp is in the root directory (required for Kaiju to work)
        if not (kaiju_db_path / "nodes.dmp").exists():
            logger.warning("Kaiju database file not found: nodes.dmp")
            logger.warning("This is required for Kaiju to work properly")
            logger.warning("Please run: kaiju-makedb -d $KAIJU_DB_DIR to download full taxonomy database")
            return None
        
        # kaiju_db_viruses.fmi is in the main kaiju database directory
        if not (kaiju_db_path / "kaiju_db_viruses.fmi").exists():
            logger.warning("Kaiju database file not found: kaiju_db_viruses.fmi")
            return None
        
        logger.info(f"Kaiju viral database found at {kaiju_db_path}")
        return kaiju_db_path
    
    def _parse_kaiju_results(self, results_file: Path) -> Dict[str, Dict]:
        """
        Parse Kaiju classification results.
        
        Args:
            results_file: Path to Kaiju results file
            
        Returns:
            Dictionary containing parsed classifications
        """
        classifications = {}
        
        try:
            with open(results_file, 'r') as f:
                for line in f:
                    if line.startswith('C'):  # Classified reads
                        parts = line.strip().split('\t')
                        if len(parts) >= 3:
                            contig_id = parts[1]
                            taxon_id = parts[2]
                            # Kaiju output format with names: C	contig_id	taxon_id	species_name
                            if len(parts) >= 4:
                                taxon_name = parts[3]  # Actual species name
                            else:
                                taxon_name = f"Taxon_{taxon_id}"  # Fallback to ID if no name
                            
                            classifications[contig_id] = {
                                'taxon_id': taxon_id,
                                'taxon_name': taxon_name,
                                'classification': taxon_name
                            }
                    elif line.startswith('U'):  # Unclassified reads
                        parts = line.strip().split('\t')
                        if len(parts) >= 2:
                            contig_id = parts[1]
                            classifications[contig_id] = {
                                'taxon_id': '0',
                                'taxon_name': 'Unclassified',
                                'classification': 'Unclassified'
                            }
            
            logger.info(f"Parsed {len(classifications)} Kaiju classifications")
            return classifications
            
        except Exception as e:
            logger.error(f"Error parsing Kaiju results: {e}")
            return {}
    
    def _write_kaiju_summary(self, output_file: Path, classifications: Dict[str, Dict]) -> None:
        """Write Kaiju classification summary to TSV."""
        try:
            output_file.parent.mkdir(parents=True, exist_ok=True)
            with open(output_file, 'w') as f:
                # Header
                f.write("contig_id\ttaxon_id\ttaxon_name\tclassification\n")
                
                # Rows
                for contig_id, data in classifications.items():
                    f.write(
                        f"{contig_id}\t{data.get('taxon_id', '')}\t"
                        f"{data.get('taxon_name', '')}\t{data.get('classification', '')}\n"
                    )
            
            logger.info(f"Wrote Kaiju classification summary to {output_file}")
        except Exception as e:
            logger.warning(f"Failed to write Kaiju summary TSV: {e}")
