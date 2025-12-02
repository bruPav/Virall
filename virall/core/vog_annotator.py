"""
VOG (Viral Orthologous Groups) annotation module for viral gene classification.

This module provides comprehensive viral gene annotation using the VOG database,
which contains curated viral protein families and functional annotations.
"""

import os
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import pandas as pd
from loguru import logger
from Bio import SeqIO
# Note: Using HMMER instead of BLAST for VOG annotation


class VOGAnnotator:
    """VOG database annotator for viral gene classification."""
    
    def __init__(self, vog_db_path: Optional[str] = None, config: Optional[Dict] = None, threads: int = 4):
        """
        Initialize VOG annotator.
        
        Args:
            vog_db_path: Path to VOG database directory
            config: Configuration dictionary
            threads: Number of threads for BLAST searches
        """
        self.threads = threads
        self.config = config or {}
        
        if vog_db_path:
            self.vog_db_path = Path(os.path.expanduser(vog_db_path))
        else:
            # ALWAYS check container bind mount paths first, even if config has a path
            # This ensures bind mounts work correctly in Singularity containers
            container_paths = [
                Path("/opt/virall/databases/vog_db"),
                Path("/opt/virall/src/databases/vog_db")
            ]
            found_path = None
            for container_path in container_paths:
                if container_path.exists() and (container_path / "vog.hmm").exists():
                    found_path = container_path
                    logger.debug(f"Found VOG database at container path: {container_path}")
                    break
            
            if found_path:
                self.vog_db_path = found_path
            else:
                # Try to get from config, but only if it exists
                config_path = self.config.get('databases', {}).get('vog_db_path')
                if config_path:
                    config_vog_path = Path(config_path)
                    if config_vog_path.exists() and (config_vog_path / "vog.hmm").exists():
                        self.vog_db_path = config_vog_path
                        logger.debug(f"Using VOG database path from config: {self.vog_db_path}")
                    else:
                        logger.debug(f"Config path {config_vog_path} does not exist, using _find_vog_database()")
                        self.vog_db_path = self._find_vog_database()
                else:
                    self.vog_db_path = self._find_vog_database()
        self.vog_fasta = None
        self.vog_annotations = None
        self.vog_hmm = None
        
        if self.vog_db_path and self.vog_db_path.exists():
            self._load_vog_database()
        else:
            # Don't show warning immediately - try to find database first
            if not self.vog_db_path:
                # Try to find database if no path was provided
                self.vog_db_path = self._find_vog_database()
                if self.vog_db_path and self.vog_db_path.exists():
                    self._load_vog_database()
                else:
                    logger.warning("VOG database not found. Run setup_vog_database() first.")
            else:
                logger.warning("VOG database not found. Run setup_vog_database() first.")
    
    def _run_command(self, cmd: List[str], log_file: Optional[Path] = None, **kwargs) -> subprocess.CompletedProcess:
        """
        Run a command and stream output to a log file or logger to avoid memory issues.
        
        Args:
            cmd: Command to run
            log_file: Path to log file (optional)
            **kwargs: Additional arguments for subprocess.run
            
        Returns:
            CompletedProcess object
        """
        cmd_str = " ".join(cmd)
        logger.info(f"Running command: {cmd_str}")
        
        if log_file:
            # Ensure parent directory exists
            if isinstance(log_file, str):
                log_file = Path(log_file)
            log_file.parent.mkdir(parents=True, exist_ok=True)
            
            with open(log_file, "w") as f:
                # Stream stdout and stderr to the file
                # Remove capture_output if present in kwargs to avoid conflict
                kwargs.pop('capture_output', None)
                kwargs.pop('stdout', None)
                kwargs.pop('stderr', None)
                kwargs.pop('text', None)
                
                return subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, text=True, **kwargs)
        else:
            return subprocess.run(cmd, **kwargs)
    
    def _find_vog_database(self) -> Path:
        """Find VOG database in common locations."""
        # Check for VIRALL_DATABASE_DIR environment variable first
        env_db_dir = os.environ.get("VIRALL_DATABASE_DIR")
        if env_db_dir:
            env_vog_path = Path(env_db_dir) / "vog_db"
            if env_vog_path.exists() and (env_vog_path / "vog.hmm").exists():
                return env_vog_path

        # Check common container database locations first (where virall setup-db creates them)
        container_paths = [
            Path("/opt/virall/databases/vog_db"),
            Path("/opt/virall/src/databases/vog_db")
        ]
        for container_path in container_paths:
            if container_path.exists() and (container_path / "vog.hmm").exists():
                return container_path
        
        # Try current working directory first, then fall back to installation directory
        cwd_db_path = Path.cwd() / "databases" / "vog_db"
        if cwd_db_path.exists() and (cwd_db_path / "vog.hmm").exists():
            return cwd_db_path
        
        # Use the same installation directory detection as assembler
        software_dir = self._find_installation_directory()
        
        possible_paths = [
            software_dir / "databases" / "vog_db",
            Path.home() / "vog_db",
            Path.home() / ".vog_db",
            Path("/opt/vog_db"),
            Path("/usr/local/vog_db"),
            Path.cwd() / "vog_db"
        ]
        
        for path in possible_paths:
            if path.exists() and (path / "vog.hmm").exists():
                # Also check for annotations file in hmm subdirectory
                annotations_file = path / "hmm" / "vog.annotations.tsv"
                if annotations_file.exists():
                    return path
                # If no annotations file, still return the path as it has vog.hmm
                return path
        
        return software_dir / "databases" / "vog_db"  # Default to installation directory
    
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
    
    def _load_vog_database(self):
        """Load VOG database files."""
        try:
            self.vog_hmm = self.vog_db_path / "vog.hmm"
            
            if not self.vog_hmm.exists():
                logger.info(f"VOG HMM database not found in {self.vog_db_path} - will be downloaded during setup")
                self.vog_db_path = None
                return
            
            # Load VOG annotations for function names
            self.vog_annotations = self._load_vog_annotations()
            
            logger.info(f"VOG HMM database loaded from {self.vog_db_path}")
            
        except Exception as e:
            logger.error(f"Error loading VOG database: {e}")
            self.vog_db_path = None
    
    def _load_vog_annotations(self):
        """Load VOG annotations file for function names."""
        try:
            annotations_file = self.vog_db_path / "hmm" / "vog.annotations.tsv"
            if not annotations_file.exists():
                logger.warning(f"VOG annotations file not found: {annotations_file}")
                return {}
            
            vog_dict = {}
            with open(annotations_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 5:
                        vog_id = parts[0]
                        protein_count = parts[1]
                        species_count = parts[2]
                        functional_category = parts[3]
                        consensus_function = parts[4]
                        
                        vog_dict[vog_id] = {
                            'name': vog_id,
                            'function': consensus_function,
                            'category': functional_category,
                            'consensus_function': consensus_function,
                            'consensus_annotation': consensus_function
                        }
            
            logger.info(f"Loaded {len(vog_dict)} VOG annotations")
            return vog_dict
            
        except Exception as e:
            logger.warning(f"Error loading VOG annotations: {e}")
            return {}
    
    def setup_vog_database(self, download_path: Optional[str] = None) -> bool:
        """
        Download and setup VOG database.
        
        Args:
            download_path: Path to download VOG database
            
        Returns:
            True if setup successful, False otherwise
        """
        if download_path:
            self.vog_db_path = Path(os.path.expanduser(download_path))
        else:
            # Default to software directory
            software_dir = Path(__file__).parent.parent.parent.parent
            self.vog_db_path = software_dir / "databases" / "vog_db"
        
        self.vog_db_path.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Setting up VOG database in {self.vog_db_path}")
        
        try:
            # Download VOG database files from alternative source (ViralRecall)
            vog_urls = {
                "hmm.tar.gz": "https://zenodo.org/records/12666277/files/hmm.tar.gz?download=1"
            }
            
            for filename, url in vog_urls.items():
                file_path = self.vog_db_path / filename
                if not file_path.exists():
                    logger.info(f"Downloading {filename}...")
                    # Use curl for better macOS compatibility
                    subprocess.run([
                        "curl", "-L", "-o", str(file_path), url
                    ], check=True)
                    logger.info(f"Downloaded {filename}")
                else:
                    logger.info(f"{filename} already exists, skipping download")
            
            # Extract the tar.gz file
            logger.info("Extracting VOG database files...")
            subprocess.run([
                "tar", "-xzf", str(self.vog_db_path / "hmm.tar.gz"), 
                "-C", str(self.vog_db_path)
            ], check=True)
            
            # Move HMM files to the correct location
            hmm_dir = self.vog_db_path / "hmm"
            if hmm_dir.exists():
                # Copy vogdb.hmm to the main directory
                vog_hmm_source = hmm_dir / "vogdb.hmm"
                if vog_hmm_source.exists():
                    vog_hmm_target = self.vog_db_path / "vog.hmm"
                    subprocess.run([
                        "cp", str(vog_hmm_source), str(vog_hmm_target)
                    ], check=True)
                    logger.info("VOG HMM database extracted successfully")
                    
                    # Press the HMM database for HMMER (only if not already pressed)
                    hmm_pressed_file = vog_hmm_target.with_suffix('.hmm.h3m')
                    if hmm_pressed_file.exists():
                        logger.info("HMM database already pressed, skipping")
                    else:
                        logger.info("Pressing HMM database for HMMER...")
                        subprocess.run([
                            "hmmpress", str(vog_hmm_target)
                        ], check=True)
                        logger.info("HMM database pressed successfully")
                else:
                    logger.warning("vogdb.hmm not found in extracted files")
                    return False
            else:
                logger.warning("HMM directory not found in extracted files")
                return False
            
            # Note: This version only provides HMM files, not FASTA sequences
            # We'll use HMMER for annotation instead of BLAST
            logger.info("VOG database setup completed (HMM-only version)")
            
            # Load database
            self._load_vog_database()
            
            logger.info("VOG database setup completed successfully")
            return True
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Error setting up VOG database: {e}")
            return False
        except Exception as e:
            logger.error(f"Unexpected error during VOG setup: {e}")
            return False
    
    def annotate_proteins(self, protein_file: str, output_dir: str) -> Dict[str, Dict]:
        """
        Annotate proteins using VOG HMM database.
        
        Args:
            protein_file: Path to FASTA file with protein sequences
            output_dir: Output directory for results
            
        Returns:
            Dictionary with annotation results
        """
        if not self.vog_db_path or not self.vog_db_path.exists():
            logger.error("VOG database not available")
            return {}
        
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Annotating proteins from {protein_file} using VOG HMM database")
        
        try:
            # Run HMMER search against VOG HMM database
            hmm_output = output_dir / "vog_hmm.txt"
            domtbl_output = output_dir / "vog_domains.txt"
            hmm_cmd = [
                "hmmscan",
                "--cpu", str(self.threads),
                "--tblout", str(hmm_output),
                "--domtblout", str(domtbl_output),
                str(self.vog_hmm),
                protein_file
            ]
            
            logger.info("Running HMMER search against VOG HMM database...")
            log_file = output_dir / "hmmer.log"
            result = self._run_command(hmm_cmd, log_file=log_file)
            
            if result.returncode != 0:
                logger.error(f"HMMER search failed. See log at {log_file}")
                return {}
            
            # Parse HMMER results
            annotations = self._parse_hmmer_results(hmm_output)
            # Parse domain table to compute coverage and merge
            try:
                coverage_map = self._parse_hmmer_domtbl(domtbl_output)
                for qid, hit in annotations.items():
                    key = (qid, hit.get('vog_id', ''))
                    if key in coverage_map:
                        hit['coverage'] = round(coverage_map[key] * 100.0, 1)
            except Exception as _:
                # Leave coverage as N/A if domtbl parsing fails
                pass
            
            # Write results
            self._write_annotation_results(annotations, output_dir)
            
            logger.info(f"VOG HMM annotation completed. Found {len(annotations)} hits.")
            return annotations
            
        except Exception as e:
            logger.error(f"Error during VOG annotation: {e}")
            return {}
    
    def _parse_hmmer_results(self, hmm_output: Path) -> Dict[str, Dict]:
        """Parse HMMER tabular output results."""
        annotations = {}
        
        try:
            with open(hmm_output, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    
                    parts = line.strip().split()
                    if len(parts) < 15:
                        continue
                    
                    query_id = parts[2]
                    vog_id = parts[0]
                    evalue = float(parts[4])
                    score = float(parts[5])
                    bias = float(parts[6])
                    hmm_start = int(float(parts[7]))  # Convert to int via float first
                    hmm_end = int(float(parts[8]))
                    query_start = int(float(parts[9]))
                    query_end = int(float(parts[10]))
                    env_start = int(float(parts[11]))
                    env_end = int(float(parts[12]))
                    acc = float(parts[13])
                    description = ' '.join(parts[15:]) if len(parts) > 15 else 'Unknown'
                    
                    # Only keep significant hits
                    if evalue < 1e-5:
                        if query_id not in annotations or evalue < annotations[query_id]['evalue']:
                            annotations[query_id] = {
                                'vog_id': vog_id,
                                'description': description,
                                'evalue': evalue,
                                'score': score,
                                'bias': bias,
                                'hmm_start': hmm_start,
                                'hmm_end': hmm_end,
                                'query_start': query_start,
                                'query_end': query_end,
                                'env_start': env_start,
                                'env_end': env_end,
                                'accuracy': acc,
                                'vog_name': vog_id,
                                'vog_function': description,
                                'vog_category': 'VOG',
                                'vog_consensus_function': description,
                                'vog_consensus_annotation': description
                            }
                            
        except Exception as e:
            logger.error(f"Error parsing HMMER results: {e}")
        
        return annotations

    def _parse_hmmer_domtbl(self, domtbl_file: Path) -> Dict[Tuple[str, str], float]:
        """Parse HMMER domtblout to compute query coverage per (query_id, vog_id).
        Returns a mapping (query_id, vog_id) -> coverage_fraction (0..1).
        """
        coverage: Dict[Tuple[str, str], float] = {}
        # For multiple domains per (query, target), keep the maximum aligned span proportion
        try:
            with open(domtbl_file, 'r') as f:
                for line in f:
                    if not line or line.startswith('#'):
                        continue
                    parts = line.strip().split()
                    # Expect at least 23 columns in domtblout
                    if len(parts) < 23:
                        continue
                    target = parts[0]  # VOG ID
                    query = parts[3]   # protein (query) id
                    try:
                        qlen = int(parts[5])
                        ali_from = int(parts[17])
                        ali_to = int(parts[18])
                    except ValueError:
                        continue
                    if qlen <= 0:
                        continue
                    aligned_len = max(0, ali_to - ali_from + 1)
                    frac = min(1.0, aligned_len / float(qlen))
                    key = (query, target)
                    # Keep the maximum coverage observed across domains
                    if key not in coverage or frac > coverage[key]:
                        coverage[key] = frac
        except FileNotFoundError:
            return {}
        return coverage

    
    
    def _write_annotation_results(self, annotations: Dict[str, Dict], output_dir: Path):
        """Write annotation results to files."""
        # Write detailed annotations (resolve VOG IDs to human-readable function/category when available)
        annotations_file = output_dir / "vog_annotations.tsv"
        with open(annotations_file, 'w') as f:
            f.write("query_id\tvog_id\tvog_name\tfunction\tcategory\tconsensus_function\tconsensus_annotation\tevalue\tidentity\tcoverage\tscore\n")
            
            for query_id, hit in annotations.items():
                vog_id_write = hit.get('vog_id', 'N/A')
                vog_name_write = hit.get('vog_name', 'Unknown')
                function_write = hit.get('vog_function', 'Unknown')
                category_write = hit.get('vog_category', 'Unknown')
                consensus_function_write = hit.get('vog_consensus_function', function_write)
                consensus_annotation_write = hit.get('vog_consensus_annotation', function_write)

                # Prefer resolved annotations from loaded VOG mapping
                if isinstance(self.vog_annotations, dict) and vog_id_write in self.vog_annotations:
                    vog_info = self.vog_annotations[vog_id_write]
                    function_write = vog_info.get('function', function_write)
                    category_write = vog_info.get('category', category_write)
                    consensus_function_write = vog_info.get('consensus_function', function_write)
                    consensus_annotation_write = vog_info.get('consensus_annotation', function_write)

                row = [
                    str(query_id),
                    str(vog_id_write),
                    str(vog_name_write),
                    str(function_write),
                    str(category_write),
                    str(consensus_function_write),
                    str(consensus_annotation_write),
                    str(hit.get('evalue', 'N/A')),
                    str(hit.get('identity', 'N/A')),
                    str(hit.get('coverage', 'N/A')),
                    str(hit.get('score', 'N/A')),
                ]
                f.write("\t".join(row) + "\n")
        
        # Write summary (resolve function/category when available)
        summary_file = output_dir / "vog_summary.tsv"
        with open(summary_file, 'w') as f:
            f.write("vog_id\tvog_name\tfunction\tcategory\tcount\n")
            
            # Count occurrences of each VOG
            vog_counts = {}
            for hit in annotations.values():
                vog_id = hit.get('vog_id', 'Unknown')
                if vog_id not in vog_counts:
                    vog_counts[vog_id] = {
                        'name': hit.get('vog_name', 'Unknown'),
                        'function': hit.get('vog_function', 'Unknown'),
                        'category': hit.get('vog_category', 'Unknown'),
                        'count': 0
                    }
                vog_counts[vog_id]['count'] += 1
            
            for vog_id, info in vog_counts.items():
                function_write = info['function']
                category_write = info['category']
                if isinstance(self.vog_annotations, dict) and vog_id in self.vog_annotations:
                    vog_info = self.vog_annotations[vog_id]
                    function_write = vog_info.get('function', function_write)
                    category_write = vog_info.get('category', category_write)
                f.write("\t".join([
                    str(vog_id),
                    str(info['name']),
                    str(function_write),
                    str(category_write),
                    str(info['count'])
                ]) + "\n")
        
        logger.info(f"VOG annotation results written to {output_dir}")
    
    def get_viral_classification(self, annotations: Dict[str, Dict]) -> Dict[str, str]:
        """
        Infer viral classification based on VOG annotations.
        
        Args:
            annotations: VOG annotation results (ignored, we'll read from TSV file)
            
        Returns:
            Dictionary mapping contig IDs to viral classifications
        """
        # Group annotations by contig ID
        contig_classifications = {}
        
        # Read the processed VOG annotations from the TSV file
        # The raw HMMER results don't have the processed function information
        try:
            import pandas as pd
            from pathlib import Path
            
            # Find the VOG annotations TSV file
            vog_tsv_file = None
            for parent_dir in [Path.cwd(), Path.cwd().parent]:
                for tsv_file in parent_dir.rglob("vog_annotations.tsv"):
                    vog_tsv_file = tsv_file
                    break
                if vog_tsv_file:
                    break
            
            if vog_tsv_file and vog_tsv_file.exists():
                df = pd.read_csv(vog_tsv_file, sep='\t')
                
                for _, row in df.iterrows():
                    query_id = row['query_id']
                    function = row['function']
                    
                    # Extract contig ID from gene ID (remove gene suffix like _1, _2, etc.)
                    contig_id = query_id.rsplit('_', 1)[0] if '_' in query_id else query_id
                    
                    # Extract human-readable function from VOG annotation
                    if function and function != '1 1 1 -':  # Skip generic placeholders
                        # Parse the function to get the human-readable name
                        # Format: "sp|P0C6T8|R1A_CVBEN Replicase polyprotein 1a"
                        # We want: "Replicase polyprotein 1a"
                        if '|' in function:
                            # Split by '|' and take the last part (after the third |)
                            parts = function.split('|')
                            if len(parts) >= 4:
                                human_readable = parts[-1].strip()
                                contig_classifications[contig_id] = human_readable
                            else:
                                contig_classifications[contig_id] = function
                        else:
                            contig_classifications[contig_id] = function
                    else:
                        contig_classifications[contig_id] = 'Viral_unknown'
            else:
                # Fallback to raw annotations if TSV file not found
                for query_id, hit in annotations.items():
                    contig_id = query_id.rsplit('_', 1)[0] if '_' in query_id else query_id
                    contig_classifications[contig_id] = 'Viral_unknown'
                    
        except Exception as e:
            # Fallback to raw annotations if TSV parsing fails
            for query_id, hit in annotations.items():
                contig_id = query_id.rsplit('_', 1)[0] if '_' in query_id else query_id
                contig_classifications[contig_id] = 'Viral_unknown'
        
        return contig_classifications
    
    def cleanup(self):
        """Clean up temporary files."""
        # Remove any temporary files if needed
        pass
