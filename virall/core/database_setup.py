"""
Database setup module for Virall.

This module handles downloading and formatting all required databases:
- VOG: Viral gene annotation database
- Kaiju: Viral taxonomy classification database
- CheckV: Viral contig quality assessment database
"""

import os
import subprocess
import time
import tarfile
import urllib.request
from pathlib import Path
from typing import Dict, Optional, Tuple
from loguru import logger

from .vog_annotator import VOGAnnotator


class DatabaseSetup:
    """Handles setup and formatting of all Virall databases."""
    
    def __init__(self, base_dir: Optional[str] = None, progress_callback=None):
        """
        Initialize database setup.
        
        Args:
            base_dir: Base directory for databases. If None, uses installation directory or current directory.
            progress_callback: Optional callback function(message: str) to display progress messages.
        """
        self.base_dir = self._determine_base_dir(base_dir)
        self.progress_callback = progress_callback or (lambda msg: None)
        # Try to create directory, but handle read-only filesystem gracefully
        try:
            self.base_dir.mkdir(parents=True, exist_ok=True)
        except OSError as e:
            if e.errno == 30:  # Read-only file system
                # In container with read-only filesystem, try common bind mount locations
                logger.warning(f"Cannot create directory at {self.base_dir} (read-only filesystem)")
                # Try common bind mount locations
                for bind_path in [Path("/opt/virall/databases"), Path("/opt/virall/src/databases")]:
                    if bind_path.exists() and os.access(bind_path, os.W_OK):
                        logger.info(f"Using bind-mounted directory: {bind_path}")
                        self.base_dir = bind_path
                        break
                else:
                    # Fall back to current working directory
                    logger.warning("No writable bind mount found, using current directory")
                    self.base_dir = Path.cwd() / "databases"
                    try:
                        self.base_dir.mkdir(parents=True, exist_ok=True)
                    except OSError:
                        raise OSError(f"Cannot create database directory. Please specify --base-dir with a writable path.")
            else:
                raise
        
        # Set up LD_LIBRARY_PATH for HMMER (needs OpenMPI libraries)
        self._setup_mpi_library_path()
    
    def _determine_base_dir(self, base_dir: Optional[str]) -> Path:
        """Determine base directory for databases."""
        if base_dir:
            return Path(base_dir)
        
        # Check for VIRALL_DATABASE_DIR environment variable first
        env_db_dir = os.environ.get("VIRALL_DATABASE_DIR")
        if env_db_dir:
            env_path = Path(env_db_dir)
            # If it exists and is writable, use it
            if env_path.exists() and os.access(env_path, os.W_OK):
                return env_path
            # If it doesn't exist but parent is writable, use it (we'll create it)
            if not env_path.exists() and env_path.parent.exists() and os.access(env_path.parent, os.W_OK):
                return env_path
        
        # Check for common bind mount locations first (for containers)
        for bind_path in [Path("/opt/virall/databases"), Path("/opt/virall/src/databases")]:
            if bind_path.exists() and os.access(bind_path, os.W_OK):
                logger.info(f"Found writable bind mount: {bind_path}")
                return bind_path
        
        # Try to find installation directory
        try:
            from .validator import AssemblyValidator
            validator = AssemblyValidator()
            installation_dir = validator._find_installation_directory()
            return installation_dir / "databases"
        except:
            return Path.cwd() / "databases"
    
    def _setup_mpi_library_path(self):
        """Set up LD_LIBRARY_PATH for HMMER (needs OpenMPI libraries)."""
        conda_prefix = os.environ.get('CONDA_PREFIX', '')
        if conda_prefix:
            lib_path = Path(conda_prefix) / "lib"
            if lib_path.exists():
                current_ld_path = os.environ.get('LD_LIBRARY_PATH', '')
                new_ld_path = f"{lib_path}:{current_ld_path}" if current_ld_path else str(lib_path)
                os.environ['LD_LIBRARY_PATH'] = new_ld_path
    
    def setup_vog_database(self, vog_db_path: Optional[str] = None) -> Tuple[bool, str]:
        """
        Set up VOG database.
        
        Args:
            vog_db_path: Custom path for VOG database. If None, uses base_dir/vog_db.
        
        Returns:
            Tuple of (success: bool, message: str)
        """
        vog_path = Path(vog_db_path) if vog_db_path else self.base_dir / "vog_db"
        vog_path.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Setting up VOG database at {vog_path}")
        self.progress_callback("Starting VOG database download...")
        
        try:
            # Set LD_LIBRARY_PATH for subprocess calls
            self._setup_mpi_library_path()
            
            annotator = VOGAnnotator()
            self.progress_callback("Downloading VOG database (~1GB, this may take 5-10 minutes)...")
            if annotator.setup_vog_database(str(vog_path)):
                self.progress_callback("Extracting VOG database...")
                self.progress_callback("Pressing HMM database for HMMER...")
                # Verify the database was pressed
                if (vog_path / "vog.hmm.h3m").exists():
                    return True, f"VOG database setup completed successfully at {vog_path}"
                else:
                    return False, f"VOG database setup reported success but .h3m file not found. You may need to press it manually: hmmpress {vog_path / 'vog.hmm'}"
            else:
                return False, "VOG database setup failed"
        except Exception as e:
            logger.error(f"VOG database setup error: {e}")
            return False, f"VOG database setup failed: {e}"
    
    def setup_kaiju_database(self, kaiju_db_path: Optional[str] = None) -> Tuple[bool, str]:
        """
        Set up Kaiju database.
        
        Args:
            kaiju_db_path: Custom path for Kaiju database. If None, uses base_dir/kaiju_db.
        
        Returns:
            Tuple of (success: bool, message: str)
        """
        kaiju_path = Path(kaiju_db_path) if kaiju_db_path else self.base_dir / "kaiju_db"
        kaiju_path.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Setting up Kaiju database at {kaiju_path}")
        self.progress_callback("Starting Kaiju database download...")
        
        try:
            # Check if kaiju-makedb is available
            if subprocess.run(["which", "kaiju-makedb"], capture_output=True).returncode != 0:
                return False, "kaiju-makedb not found. Please install Kaiju first."
            
            # Download Kaiju viral database
            # Note: kaiju-makedb doesn't have a -d option, it creates files in the current working directory
            logger.info("Downloading Kaiju viral database...")
            self.progress_callback("Downloading Kaiju viral database (this may take 10-30 minutes)...")
            self.progress_callback("Note: This downloads the full NCBI taxonomy database first, then creates viral subset")
            result = subprocess.run(
                ["kaiju-makedb", "-s", "viruses"],
                cwd=str(kaiju_path),  # Change to target directory before running
                capture_output=False,  # Show output in real-time
                text=True
            )
            
            if result.returncode == 0:
                # Clean up temporary files
                # Note: kaiju-makedb creates files in a 'viruses/' subdirectory
                logger.info("Cleaning up temporary database files...")
                self.progress_callback("Cleaning up temporary database files...")
                viruses_dir = kaiju_path / "viruses"
                if viruses_dir.exists():
                    for temp_file in ["kaiju_db_viruses.bwt", "kaiju_db_viruses.sa"]:
                        temp_path = viruses_dir / temp_file
                        if temp_path.exists():
                            temp_path.unlink()
                # Also check root directory
                for temp_file in ["kaiju_db_viruses.bwt", "kaiju_db_viruses.sa"]:
                    temp_path = kaiju_path / temp_file
                    if temp_path.exists():
                        temp_path.unlink()
                
                # Verify nodes.dmp exists (required for Kaiju)
                # Note: kaiju-makedb -s viruses should already download taxonomy files
                # But if nodes.dmp is missing, we need to download it manually
                if not (kaiju_path / "nodes.dmp").exists():
                    logger.info("nodes.dmp not found. Downloading taxonomy database manually...")
                    self.progress_callback("Downloading full taxonomy database for Kaiju...")
                    # Download taxonomy manually since kaiju-makedb doesn't have a -d option
                    taxdump_url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
                    taxdump_file = kaiju_path / "taxdump.tar.gz"
                    try:
                        import urllib.request
                        urllib.request.urlretrieve(taxdump_url, taxdump_file)
                        import tarfile
                        with tarfile.open(taxdump_file, "r:gz") as tar:
                            tar.extractall(kaiju_path)
                        taxdump_file.unlink()  # Remove tar file
                        logger.info("Taxonomy database downloaded successfully")
                    except Exception as e:
                        logger.error(f"Taxonomy database download failed: {e}")
                        return False, f"Taxonomy database download failed: {e}"
                
                # Verify database was created
                # kaiju-makedb creates files in a 'viruses/' subdirectory
                viruses_dir = kaiju_path / "viruses"
                db_file = viruses_dir / "kaiju_db_viruses.fmi" if viruses_dir.exists() else kaiju_path / "kaiju_db_viruses.fmi"
                
                if db_file.exists():
                    msg = f"Kaiju database setup completed successfully at {kaiju_path}"
                    if (kaiju_path / "nodes.dmp").exists():
                        msg += " (taxonomy database ready)"
                    return True, msg
                else:
                    return False, f"Kaiju database setup may be incomplete: {db_file} not found"
            else:
                return False, f"Kaiju database setup failed: {result.stderr}"
        except Exception as e:
            logger.error(f"Kaiju database setup error: {e}")
            return False, f"Kaiju database setup failed: {e}"
    
    def setup_checkv_database(self, checkv_db_path: Optional[str] = None) -> Tuple[bool, str]:
        """
        Set up CheckV database.
        
        Args:
            checkv_db_path: Custom path for CheckV database. If None, uses base_dir/checkv_db.
        
        Returns:
            Tuple of (success: bool, message: str)
        """
        checkv_path = Path(checkv_db_path) if checkv_db_path else self.base_dir / "checkv_db"
        checkv_path.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Setting up CheckV database at {checkv_path}")
        self.progress_callback("Starting CheckV database download...")
        
        try:
            # Check if checkv command is available
            if subprocess.run(["which", "checkv"], capture_output=True).returncode != 0:
                return False, "checkv command not found. Please install CheckV first."
            
            # Try downloading with retry mechanism
            MAX_RETRIES = 3
            RETRY_COUNT = 0
            success = False
            
            while RETRY_COUNT < MAX_RETRIES and not success:
                if RETRY_COUNT > 0:
                    logger.info(f"Retry attempt {RETRY_COUNT + 1} of {MAX_RETRIES}...")
                    self.progress_callback(f"Retry attempt {RETRY_COUNT + 1} of {MAX_RETRIES}...")
                    time.sleep(60)  # Wait 60 seconds before retry
                
                logger.info("Downloading CheckV database...")
                self.progress_callback(f"Downloading CheckV database (~3GB, attempt {RETRY_COUNT + 1}, this may take 10-30 minutes)...")
                result = subprocess.run(
                    ["checkv", "download_database", str(checkv_path)],
                    capture_output=False,  # Show output in real-time
                    text=True
                )
                
                if result.returncode == 0:
                    # Check if database directory has content
                    if checkv_path.exists() and any(checkv_path.iterdir()):
                        return True, f"CheckV database setup completed successfully at {checkv_path}"
                    else:
                        logger.warning("CheckV database directory is empty")
                        RETRY_COUNT += 1
                else:
                    RETRY_COUNT += 1
                    if RETRY_COUNT < MAX_RETRIES:
                        logger.warning(f"Download failed: {result.stderr}")
                    else:
                        # Try manual download as fallback
                        return self._setup_checkv_manual(checkv_path)
            
            return False, "CheckV database setup failed after all attempts"
        except Exception as e:
            logger.error(f"CheckV database setup error: {e}")
            return False, f"CheckV database setup failed: {e}"
    
    def _setup_checkv_manual(self, checkv_path: Path) -> Tuple[bool, str]:
        """Manual CheckV database download as fallback."""
        try:
            logger.info("Trying manual CheckV database download...")
            self.progress_callback("Trying manual CheckV database download...")
            checkv_url = "https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz"
            checkv_tar = checkv_path / "checkv-db.tar.gz"
            
            logger.info(f"Downloading from {checkv_url}...")
            self.progress_callback(f"Downloading from {checkv_url}...")
            urllib.request.urlretrieve(checkv_url, checkv_tar)
            
            # Verify and extract
            if checkv_tar.exists() and checkv_tar.stat().st_size > 0:
                logger.info("Extracting CheckV database...")
                self.progress_callback("Extracting CheckV database...")
                with tarfile.open(checkv_tar, 'r:gz') as tar:
                    tar.extractall(checkv_path)
                checkv_tar.unlink()
                
                if any(checkv_path.iterdir()):
                    return True, f"CheckV database downloaded and extracted successfully at {checkv_path}"
                else:
                    return False, "CheckV database extraction failed"
            else:
                return False, "CheckV database download failed"
        except Exception as e:
            logger.error(f"Manual CheckV download error: {e}")
            return False, f"Manual CheckV download failed: {e}"
    
    def setup_all_databases(self, 
                           vog_db: Optional[str] = None,
                           kaiju_db: Optional[str] = None,
                           checkv_db: Optional[str] = None) -> Dict[str, Tuple[bool, str]]:
        """
        Set up all databases.
        
        Args:
            vog_db: Custom path for VOG database
            kaiju_db: Custom path for Kaiju database
            checkv_db: Custom path for CheckV database
        
        Returns:
            Dictionary mapping database name to (success, message) tuple
        """
        results = {}
        
        results['vog'] = self.setup_vog_database(vog_db)
        results['kaiju'] = self.setup_kaiju_database(kaiju_db)
        results['checkv'] = self.setup_checkv_database(checkv_db)
        
        return results

