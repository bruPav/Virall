"""
Assembly validation and quality assessment module.
"""

import os
import subprocess
import tempfile
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Union
from loguru import logger

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class AssemblyValidator:
    """
    Validates and assesses the quality of assembled viral genomes using
    multiple tools and metrics. Also provides viral classification capabilities.
    """
    
    def __init__(self, threads: int = 8, config: Optional[Dict] = None):
        """Initialize the assembly validator."""
        self.threads = threads
        self.config = config or {}
        self.temp_dir = Path(tempfile.mkdtemp(prefix="viral_validate_"))
        # Setup CheckV database for viral contig quality assessment
        self.checkv_db_path = self._setup_checkv_database()
        logger.info(f"AssemblyValidator initialized with {threads} threads")
    
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
    
    def _setup_checkv_database(self) -> Optional[str]:
        """Setup CheckV database for viral contig quality assessment."""
        # Check config first
        config_path = self.config.get('databases', {}).get('checkv_db_path')
        if config_path:
            checkv_db_path = Path(config_path)
        else:
            # Check common container database locations first (where virall setup-db creates them)
            container_paths = [
                Path("/opt/virall/databases/checkv_db"),
                Path("/opt/virall/src/databases/checkv_db")
            ]
            checkv_db_path = None
            for container_path in container_paths:
                if container_path.exists():
                    checkv_db_path = container_path
                    break
            
            if not checkv_db_path:
                # Try current working directory first, then fall back to installation directory
                cwd_db_path = Path.cwd() / "databases" / "checkv_db"
                if cwd_db_path.exists():
                    checkv_db_path = cwd_db_path
                else:
                    # Use the same installation directory detection as assembler
                    software_dir = self._find_installation_directory()
                    checkv_db_path = software_dir / "databases" / "checkv_db"
            
            if not checkv_db_path.exists():
                # Fallback to home directory
                checkv_db_path = Path.home() / "checkv_db"
        
        if not checkv_db_path.exists():
            logger.info("CheckV database not found. Setting up CheckV database...")
            logger.info("This will download ~3GB of viral genome data")
            logger.info("Note: This step may take 10-30 minutes depending on your internet connection")
            
            try:
                # Create checkv database directory
                checkv_db_path.mkdir(exist_ok=True)
                
                # Download CheckV database
                cmd = ["checkv", "download_database", str(checkv_db_path)]
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                if result.returncode == 0:
                    logger.info("CheckV database downloaded successfully")
                    # CheckV downloads to a subdirectory, find the actual database path
                    actual_db_path = self._find_checkv_database_path(checkv_db_path)
                    return actual_db_path
                else:
                    logger.warning(f"CheckV database download failed: {result.stderr}")
                    return None
                    
            except Exception as e:
                logger.warning(f"Failed to setup CheckV database: {e}")
                return None
        else:
            # CheckV database exists, find the actual database path
            actual_db_path = self._find_checkv_database_path(checkv_db_path)
            if actual_db_path:
                logger.info(f"CheckV database found at {actual_db_path}")
                return actual_db_path
            else:
                logger.warning("CheckV database directory exists but database files not found")
                return None
    
    def _find_checkv_database_path(self, checkv_db_path: Path) -> Optional[str]:
        """Find the actual CheckV database path within the downloaded directory."""
        # Look for checkv-db-v* subdirectories
        for item in checkv_db_path.iterdir():
            if item.is_dir() and item.name.startswith("checkv-db-"):
                # Check if genome_db subdirectory exists
                genome_db_path = item / "genome_db"
                if genome_db_path.exists():
                    # CheckV automatically appends 'genome_db' to the path we provide
                    # So we need to return the parent directory (checkv-db-v*)
                    return str(item)
        
        # Check if the directory contains any CheckV database files
        # Look for common CheckV database files
        checkv_files = ["genome_db", "proteins", "taxonomy", "completeness"]
        has_checkv_files = any((checkv_db_path / file).exists() for file in checkv_files)
        
        if has_checkv_files:
            return str(checkv_db_path)
        
        # If no subdirectory found and no CheckV files, return None
        return None
    
    def run_checkv(
        self, 
        genome_file: str,
        output_dir: Optional[str] = None
    ) -> Dict[str, Union[str, float, int]]:
        """
        Run CheckV quality assessment on viral contigs.
        
        Args:
            genome_file: Path to viral contigs file
            output_dir: Output directory for CheckV results
            
        Returns:
            Dictionary with CheckV results
        """
        if not self.checkv_db_path:
            logger.warning("CheckV database not available, skipping CheckV analysis")
            return {
                "status": "skipped",
                "reason": "CheckV database not available"
            }
        
        if output_dir is None:
            output_dir = self.temp_dir / "checkv_output"
        
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Running CheckV on {genome_file}")
        
        # Try end_to_end first (best for complete DNA viral genomes)
        logger.info(f"Running CheckV end_to_end analysis on {genome_file}")
        
        # Convert to absolute paths to avoid subprocess issues
        genome_file_abs = os.path.abspath(genome_file)
        output_dir_abs = os.path.abspath(str(output_dir))
        checkv_db_abs = os.path.abspath(self.checkv_db_path)
        
        cmd = [
            "checkv", "end_to_end",
            genome_file_abs,
            output_dir_abs,
            "-d", checkv_db_abs,
            "-t", "1"  # Use single thread to avoid Prodigal race conditions
        ]
        
        # Debug: Log the exact command being run
        logger.info(f"CheckV command: {' '.join(cmd)}")
        logger.info(f"Working directory: {os.getcwd()}")
        logger.info(f"File exists: {os.path.exists(genome_file_abs)}")
        
        # Run with explicit working directory and environment
        result = subprocess.run(
            cmd, 
            capture_output=True, 
            text=True,
            cwd=os.getcwd(),  # Use current working directory
            env=os.environ.copy()  # Copy current environment
        )
        
        # Debug: Log the full output for troubleshooting
        logger.info(f"CheckV return code: {result.returncode}")
        if result.stdout:
            logger.info(f"CheckV stdout: {result.stdout[:500]}...")  # First 500 chars
        if result.stderr:
            logger.info(f"CheckV stderr: {result.stderr[:500]}...")  # First 500 chars
        if result.returncode == 0:
            logger.info("CheckV end_to_end completed successfully")
            return self._parse_checkv_results(output_dir, "end_to_end")
        
        # If end_to_end failed, try contamination mode (better for RNA viruses/fragments)
        logger.warning(f"CheckV end_to_end failed: {result.stderr}")
        logger.info("Trying CheckV contamination mode as fallback...")
        
        # Clean up failed output directory
        import shutil
        if output_dir.exists():
            shutil.rmtree(output_dir, ignore_errors=True)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        cmd = [
            "checkv", "contamination",
            genome_file_abs,
            output_dir_abs,
            "-d", checkv_db_abs,
            "-t", "1"  # Use single thread to avoid Prodigal race conditions
        ]
        
        # Debug: Log the exact command being run
        logger.info(f"CheckV contamination command: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd, 
            capture_output=True, 
            text=True,
            cwd=os.getcwd(),  # Use current working directory
            env=os.environ.copy()  # Copy current environment
        )
        
        # Debug: Log the full output for troubleshooting
        logger.info(f"CheckV contamination return code: {result.returncode}")
        if result.stdout:
            logger.info(f"CheckV contamination stdout: {result.stdout[:500]}...")
        if result.stderr:
            logger.info(f"CheckV contamination stderr: {result.stderr[:500]}...")
        if result.returncode == 0:
            logger.info("CheckV contamination completed successfully")
            return self._parse_checkv_results(output_dir, "contamination")
        
        # Both modes failed
        logger.error(f"CheckV failed in both end_to_end and contamination modes: {result.stderr}")
        return {
            "status": "failed",
            "error": f"CheckV failed in both modes. This may be due to Prodigal gene prediction issues with the sequence format.",
            "total_contigs": 0,
            "complete_contigs": 0,
            "high_quality_contigs": 0,
            "medium_quality_contigs": 0,
            "low_quality_contigs": 0,
            "undetermined_contigs": 0,
            "avg_completeness": 0.0,
            "avg_contamination": 0.0
        }
        
        return results
    
    def _parse_checkv_results(self, output_dir: Path, mode: str = "end_to_end") -> Dict[str, Union[str, float, int]]:
        """
        Parse CheckV results from output directory.
        
        Args:
            output_dir: CheckV output directory
            mode: CheckV mode used ("end_to_end" or "contamination")
            
        Returns:
            Dictionary with parsed CheckV results
        """
        try:
            if mode == "end_to_end":
                # Parse end_to_end results
                quality_file = output_dir / "quality_summary.tsv"
                if quality_file.exists():
                    with open(quality_file, 'r') as f:
                        lines = f.readlines()
                        if len(lines) > 1:
                            data = lines[1].strip().split('\t')
                            return {
                                "status": "completed",
                                "mode": "end_to_end",
                                "total_contigs": int(data[0]) if data[0].isdigit() else 0,
                                "complete_contigs": int(data[1]) if data[1].isdigit() else 0,
                                "high_quality_contigs": int(data[2]) if data[2].isdigit() else 0,
                                "medium_quality_contigs": int(data[3]) if data[3].isdigit() else 0,
                                "low_quality_contigs": int(data[4]) if data[4].isdigit() else 0,
                                "undetermined_contigs": int(data[5]) if data[5].isdigit() else 0,
                                "avg_completeness": float(data[6]) if data[6].replace('.', '').isdigit() else 0.0,
                                "avg_contamination": float(data[7]) if data[7].replace('.', '').isdigit() else 0.0
                            }
            
            elif mode == "contamination":
                # Parse contamination results
                contamination_file = output_dir / "contamination.tsv"
                if contamination_file.exists():
                    with open(contamination_file, 'r') as f:
                        lines = f.readlines()
                        if len(lines) > 1:
                            # Count contigs and calculate averages
                            total_contigs = len(lines) - 1
                            viral_genes_total = 0
                            host_genes_total = 0
                            
                            for line in lines[1:]:
                                data = line.strip().split('\t')
                                if len(data) >= 4:
                                    viral_genes_total += int(data[3]) if data[3].isdigit() else 0
                                    host_genes_total += int(data[4]) if data[4].isdigit() else 0
                            
                            return {
                                "status": "completed",
                                "mode": "contamination",
                                "total_contigs": total_contigs,
                                "viral_genes_total": viral_genes_total,
                                "host_genes_total": host_genes_total,
                                "avg_viral_genes": viral_genes_total / total_contigs if total_contigs > 0 else 0,
                                "avg_host_genes": host_genes_total / total_contigs if total_contigs > 0 else 0
                            }
            
            # Fallback if files don't exist
            return {
                "status": "completed",
                "mode": mode,
                "total_contigs": 0,
                "message": "CheckV completed but results could not be parsed"
            }
            
        except Exception as e:
            logger.warning(f"Error parsing CheckV results: {e}")
            return {
                "status": "completed",
                "mode": mode,
                "total_contigs": 0,
                "error": f"Results parsing failed: {e}"
            }
    
    # Assembly quality assessment methods removed
    
    def calculate_genome_statistics(self, genome_file: str) -> Dict[str, Union[int, float]]:
        """
        Calculate basic genome statistics.
        
        Args:
            genome_file: Path to genome file
            
        Returns:
            Dictionary containing genome statistics
        """
        logger.info(f"Calculating genome statistics for {genome_file}")
        
        sequences = list(SeqIO.parse(genome_file, "fasta"))
        
        if not sequences:
            return {
                "num_contigs": 0,
                "total_length": 0,
                "n50": 0,
                "longest_contig": 0,
                "shortest_contig": 0,
                "gc_content": 0.0
            }
        
        # Calculate basic statistics
        lengths = [len(seq) for seq in sequences]
        total_length = sum(lengths)
        
        # Calculate N50
        lengths_sorted = sorted(lengths, reverse=True)
        cumulative_length = 0
        n50 = 0
        for length in lengths_sorted:
            cumulative_length += length
            if cumulative_length >= total_length / 2:
                n50 = length
                break
        
        # Calculate GC content
        total_gc = 0
        total_bases = 0
        for seq in sequences:
            seq_str = str(seq.seq).upper()
            gc_count = seq_str.count('G') + seq_str.count('C')
            total_gc += gc_count
            total_bases += len(seq_str)
        
        gc_content = total_gc / total_bases if total_bases > 0 else 0.0
        
        stats = {
            "num_contigs": len(sequences),
            "total_contigs": len(sequences),  # Alias for CLI compatibility
            "total_length": total_length,
            "n50": n50,
            "longest_contig": max(lengths),
            "shortest_contig": min(lengths),
            "gc_content": gc_content,
            "average_length": total_length / len(sequences),
            "average_contig_length": total_length / len(sequences)  # Alias for CLI compatibility
        }
        
        logger.info(f"Genome statistics calculated: {stats['num_contigs']} contigs, {stats['total_length']} bp")
        
        return stats
    
    
    
    # Assembly quality assessment parsing removed
    
    def generate_quality_report(
        self, 
        genome_files: List[str],
        output_file: Optional[str] = None
    ) -> str:
        """
        Generate a comprehensive quality report for multiple genomes.
        
        Args:
            genome_files: List of genome file paths
            output_file: Output report file (optional)
            
        Returns:
            Path to generated report file
        """
        if output_file is None:
            output_file = self.temp_dir / "quality_report.html"
        
        logger.info(f"Generating quality report for {len(genome_files)} genomes")
        
        # Collect results for all genomes
        all_results = {}
        for genome_file in genome_files:
            logger.info(f"Processing {genome_file}")
            
            # Run validation tools
            stats = self.calculate_genome_statistics(genome_file)
            
            all_results[genome_file] = {
                "statistics": stats
            }
        
        # Generate HTML report
        html_content = self._generate_html_report(all_results)
        
        with open(output_file, 'w') as f:
            f.write(html_content)
        
        logger.info(f"Quality report generated: {output_file}")
        return str(output_file)
    
    def _generate_html_report(self, results: Dict) -> str:
        """Generate HTML quality report."""
        html = """
        <!DOCTYPE html>
        <html>
        <head>
            <title>Viral Genome Assembly Quality Report</title>
            <style>
                body { font-family: Arial, sans-serif; margin: 20px; }
                table { border-collapse: collapse; width: 100%; margin: 20px 0; }
                th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
                th { background-color: #f2f2f2; }
                .good { color: green; }
                .warning { color: orange; }
                .error { color: red; }
            </style>
        </head>
        <body>
            <h1>Viral Genome Assembly Quality Report</h1>
        """
        
        for genome_file, data in results.items():
            html += f"<h2>{Path(genome_file).name}</h2>"
            
            # Statistics table
            html += "<h3>Basic Statistics</h3><table>"
            stats = data.get("statistics", {})
            for key, value in stats.items():
                html += f"<tr><td>{key}</td><td>{value}</td></tr>"
            html += "</table>"
            
            # CheckV provides viral-specific quality assessment
        
        html += "</body></html>"
        return html
    
    def run_viral_classification(
        self, 
        genome_file: str,
        output_dir: Optional[str] = None,
        method: str = "blast"
    ) -> Dict[str, Union[str, float, int]]:
        """
        Run viral classification using BLAST.
        
        Args:
            genome_file: Path to genome file
            output_dir: Output directory (optional)
            method: Classification method ('blast')
            
        Returns:
            Dictionary containing classification results
        """
        if output_dir is None:
            output_dir = self.temp_dir / "viral_classification"
        
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)
        
        logger.info(f"Running viral classification on {genome_file} using BLAST")
        
        if method == "blast":
            return self._run_blast_classification(genome_file, output_dir)
        else:
            logger.warning(f"Unknown classification method: {method}, falling back to BLAST")
            return self._run_blast_classification(genome_file, output_dir)
    
    def _run_blast_classification(
        self, 
        genome_file: str, 
        output_dir: Path
    ) -> Dict[str, Union[str, float, int]]:
        """Run BLAST-based viral classification."""
        logger.info("Running BLAST-based viral classification")
        
        # Create BLAST database if it doesn't exist
        blast_db_path = self._setup_blast_database()
        if not blast_db_path:
            logger.warning("BLAST database not available, skipping classification")
            return {"status": "skipped", "reason": "No BLAST database"}
        
        # Run BLAST
        blast_output = output_dir / "blast_results.tsv"
        cmd = [
            "blastn",
            "-query", genome_file,
            "-db", blast_db_path,
            "-out", str(blast_output),
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids stitle",
            "-max_target_seqs", "10",
            "-evalue", "1e-5"
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.warning(f"BLAST failed: {result.stderr}")
            return {"status": "failed", "error": result.stderr}
        
        # Parse BLAST results
        classification_results = self._parse_blast_results(blast_output)
        classification_results["status"] = "completed"
        classification_results["output_file"] = str(blast_output)
        
        logger.info("BLAST classification completed")
        return classification_results
    
    def _setup_blast_database(self) -> Optional[str]:
        """Setup BLAST database for viral classification."""
        # Check for existing viral BLAST database
        possible_dbs = [
            "/home/ec2-user/blast_db/viral_genomes",
            "/home/ec2-user/ncbi/blast/db/viral_genomes",
            "~/blast_db/viral_genomes"
        ]
        
        for db_path in possible_dbs:
            expanded_path = str(Path(db_path).expanduser())
            if Path(expanded_path + ".nhr").exists():
                logger.info(f"Found BLAST database at {expanded_path}")
                return expanded_path
        
        logger.warning("No viral BLAST database found")
        logger.info("To create one, download viral genomes from NCBI and run:")
        logger.info("makeblastdb -in viral_genomes.fasta -dbtype nucl -out viral_genomes")
        return None
    
    def _parse_blast_results(self, blast_file: Path) -> Dict[str, Union[str, float, int]]:
        """Parse BLAST results for classification."""
        results = {
            "total_hits": 0,
            "best_hits": [],
            "taxonomic_groups": {},
            "average_identity": 0.0
        }
        
        if not blast_file.exists():
            return results
        
        try:
            df = pd.read_csv(blast_file, sep='\t', header=None, 
                           names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                                  'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                                  'evalue', 'bitscore', 'staxids', 'stitle'])
            
            if df.empty:
                return results
            
            results["total_hits"] = len(df)
            
            # Get best hits per query
            best_hits = df.groupby('qseqid').first().reset_index()
            results["best_hits"] = best_hits[['qseqid', 'sseqid', 'pident', 'evalue', 'stitle']].to_dict('records')
            
            # Calculate average identity
            results["average_identity"] = float(df['pident'].mean())
            
            # Extract taxonomic information from stitle
            for title in df['stitle'].dropna():
                # Simple taxonomic extraction (can be improved)
                if 'virus' in title.lower():
                    # Extract virus name (simplified)
                    virus_name = title.split()[0] if title.split() else "Unknown"
                    results["taxonomic_groups"][virus_name] = results["taxonomic_groups"].get(virus_name, 0) + 1
            
        except Exception as e:
            logger.warning(f"Failed to parse BLAST results: {e}")
        
        return results
    
    
    
    def cleanup(self) -> None:
        """Clean up temporary files."""
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
            logger.info("Cleaned up temporary files")
