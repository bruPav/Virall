"""
RNA-Bloom assembler for transcriptome assembly from various RNA-seq data types.
Supports single-cell, bulk RNA-seq, and long-read RNA-seq data.
"""

import os
import subprocess
import tempfile
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Union, Tuple
from loguru import logger

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class RNABloomAssembler:
    """
    RNA-Bloom assembler for comprehensive transcriptome assembly.
    
    Supports:
    - Single-cell RNA-seq (10X Genomics, etc.)
    - Bulk RNA-seq (paired-end, single-end)
    - Long-read RNA-seq (ONT, PacBio)
    - Hybrid assembly (short + long reads)
    - Reference-guided assembly
    """
    
    def __init__(self, threads: int = 8, memory: str = "16G", temp_dir: Optional[Path] = None):
        """
        Initialize RNA-Bloom assembler.
        
        Args:
            threads: Number of threads to use
            memory: Memory allocation (e.g., "16G")
            temp_dir: Temporary directory for intermediate files
        """
        self.threads = threads
        self.memory = memory
        self.temp_dir = temp_dir or Path(tempfile.mkdtemp(prefix="rnabloom_"))
        self.temp_dir.mkdir(parents=True, exist_ok=True)
        
        # Check if RNA-Bloom is available
        if not self._check_rnabloom_installation():
            raise RuntimeError("RNA-Bloom not found. Please install with: conda install -c bioconda rnabloom")
        
        logger.info(f"RNA-Bloom assembler initialized with {threads} threads, {memory} memory")
    
    def _check_rnabloom_installation(self) -> bool:
        """Check if RNA-Bloom is installed and available."""
        try:
            result = subprocess.run(["rnabloom", "-help"], capture_output=True, text=True)
            return result.returncode == 0
        except FileNotFoundError:
            return False
    
    def assemble_single_cell(
        self,
        r1_file: str,
        r2_file: str,
        output_dir: Path,
        sample_name: str = "sample",
        min_length: int = 200,
        stranded: bool = False,
        reference: Optional[str] = None,
        use_pooled_mode: bool = True
    ) -> Dict[str, str]:
        """
        Assemble single-cell RNA-seq data using RNA-Bloom.
        
        Args:
            r1_file: Path to R1 reads (cell barcodes + UMIs)
            r2_file: Path to R2 reads (cDNA sequences)
            output_dir: Output directory for assembly results
            sample_name: Sample name for output files
            min_length: Minimum transcript length
            stranded: Whether reads are strand-specific
            reference: Optional reference transcriptome for guided assembly
            use_pooled_mode: Whether to use pooled assembly mode for single-cell data
            
        Returns:
            Dictionary with assembly output file paths
        """
        logger.info(f"Assembling single-cell RNA-seq data for {sample_name}")
        
        # Create output directory
        output_dir.mkdir(parents=True, exist_ok=True)
        
        if use_pooled_mode:
            # Use pooled assembly mode for better single-cell support
            logger.info("Using RNA-Bloom pooled assembly mode for single-cell data")
            return self._assemble_single_cell_pooled(
                r1_file, r2_file, output_dir, sample_name, min_length, stranded, reference
            )
        else:
            # Use standard assembly mode
            logger.info("Using RNA-Bloom standard assembly mode for single-cell data")
            return self._assemble_single_cell_standard(
                r1_file, r2_file, output_dir, sample_name, min_length, stranded, reference
            )
    
    def _assemble_single_cell_standard(
        self,
        r1_file: str,
        r2_file: str,
        output_dir: Path,
        sample_name: str,
        min_length: int,
        stranded: bool,
        reference: Optional[str]
    ) -> Dict[str, str]:
        """Standard single-cell assembly (current approach)."""
        # Prepare RNA-Bloom command for single-cell data
        cmd = [
            "rnabloom",
            "-left", r1_file,
            "-right", r2_file,
            "-t", str(self.threads),
            "-outdir", str(output_dir),
            "-name", sample_name,
            "-length", str(min_length)
        ]
        
        # Add stranded option if specified
        if stranded:
            cmd.append("-stranded")
            logger.info("Using stranded assembly for single-cell data")
        
        # Add reference if provided
        if reference:
            cmd.extend(["-ref", reference])
            logger.info(f"Using reference-guided assembly with {reference}")
        
        # Add memory settings (convert from "16G" to "16")
        if self.memory:
            memory_value = self.memory.replace('G', '').replace('g', '').replace('M', '').replace('m', '')
            cmd.extend(["-mem", memory_value])
        
        logger.info(f"RNA-Bloom command: {' '.join(cmd)}")
        
        # Run RNA-Bloom
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(output_dir))
            
            if result.returncode != 0:
                logger.error(f"RNA-Bloom failed: {result.stderr}")
                raise RuntimeError(f"RNA-Bloom assembly failed: {result.stderr}")
            
            logger.info("RNA-Bloom single-cell assembly completed successfully")
            
            # Find and move output files to clean structure
            transcripts_file = self._organize_output_files(output_dir, sample_name, "transcripts.fa")
            transcripts_short_file = self._organize_output_files(output_dir, sample_name, "transcripts.short.fa")
            transcripts_nr_file = self._organize_output_files(output_dir, sample_name, "transcripts.nr.fa")
            
            # Return output file paths
            return {
                "transcripts": transcripts_file,
                "transcripts_short": transcripts_short_file,
                "transcripts_nr": transcripts_nr_file,
                "assembly_log": str(output_dir / "rnabloom.log")
            }
            
        except Exception as e:
            logger.error(f"RNA-Bloom assembly failed: {e}")
            raise
    
    def _assemble_single_cell_pooled(
        self,
        r1_file: str,
        r2_file: str,
        output_dir: Path,
        sample_name: str,
        min_length: int,
        stranded: bool,
        reference: Optional[str]
    ) -> Dict[str, str]:
        """Pooled assembly mode for single-cell data (preserves cell information)."""
        logger.info("Creating cell reads list for pooled assembly")
        
        # Create a reads list file for pooled assembly
        reads_list_file = output_dir / "cell_reads_list.txt"
        self._create_cell_reads_list(r1_file, r2_file, reads_list_file)
        
        # Prepare RNA-Bloom command for pooled assembly
        cmd = [
            "rnabloom",
            "-pool", str(reads_list_file),
            "-t", str(self.threads),
            "-outdir", str(output_dir),
            "-length", str(min_length)
        ]
        
        # Add stranded option if specified
        if stranded:
            cmd.append("-stranded")
            logger.info("Using stranded assembly for pooled single-cell data")
        
        # Add reference if provided
        if reference:
            cmd.extend(["-ref", reference])
            logger.info(f"Using reference-guided assembly with {reference}")
        
        # Add memory settings (convert from "16G" to "16")
        if self.memory:
            memory_value = self.memory.replace('G', '').replace('g', '').replace('M', '').replace('m', '')
            cmd.extend(["-mem", memory_value])
        
        logger.info(f"RNA-Bloom pooled command: {' '.join(cmd)}")
        
        # Run RNA-Bloom
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(output_dir))
            
            if result.returncode != 0:
                logger.error(f"RNA-Bloom failed with return code {result.returncode}")
                logger.error(f"RNA-Bloom stderr: {result.stderr}")
                logger.error(f"RNA-Bloom stdout: {result.stdout}")
                raise RuntimeError(f"RNA-Bloom assembly failed: {result.stderr}")
            
            logger.info("RNA-Bloom pooled single-cell assembly completed successfully")
            
            # Find and move output files to clean structure
            # For pooled mode, files are in subdirectories like cell_1/
            transcripts_file = self._organize_pooled_output_files(output_dir, "transcripts.fa")
            transcripts_short_file = self._organize_pooled_output_files(output_dir, "transcripts.short.fa")
            transcripts_nr_file = self._organize_pooled_output_files(output_dir, "transcripts.nr.fa")
            
            # Return output file paths
            return {
                "transcripts": transcripts_file,
                "transcripts_short": transcripts_short_file,
                "transcripts_nr": transcripts_nr_file,
                "assembly_log": str(output_dir / "rnabloom.log")
            }
            
        except Exception as e:
            logger.error(f"RNA-Bloom pooled assembly failed: {e}")
            raise
    
    def _create_cell_reads_list(self, r1_file: str, r2_file: str, output_file: Path) -> None:
        """
        Create a reads list file for RNA-Bloom pooled assembly.
        
        Args:
            r1_file: Path to R1 reads file
            r2_file: Path to R2 reads file
            output_file: Path to output reads list file
        """
        logger.info("Creating cell reads list for pooled assembly")
        
        # For now, create a simple reads list with the provided files
        # In a full implementation, this would parse cell barcodes and create separate entries
        with open(output_file, 'w') as f:
            f.write("#name left right\n")
            f.write(f"cell_1 {r1_file} {r2_file}\n")
        
        logger.info(f"Created reads list: {output_file}")
    
    def assemble_bulk_rna(
        self,
        r1_file: str,
        r2_file: str,
        output_dir: Path,
        sample_name: str = "sample",
        min_length: int = 200,
        stranded: bool = False,
        reference: Optional[str] = None,
        single_end: Optional[str] = None
    ) -> Dict[str, str]:
        """
        Assemble bulk RNA-seq data using RNA-Bloom.
        
        Args:
            r1_file: Path to R1 reads
            r2_file: Path to R2 reads
            output_dir: Output directory for assembly results
            sample_name: Sample name for output files
            min_length: Minimum transcript length
            stranded: Whether reads are strand-specific
            reference: Optional reference transcriptome for guided assembly
            single_end: Optional single-end reads file
            
        Returns:
            Dictionary with assembly output file paths
        """
        logger.info(f"Assembling bulk RNA-seq data for {sample_name}")
        
        # Create output directory
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Prepare RNA-Bloom command for bulk RNA-seq
        cmd = [
            "rnabloom",
            "-left", r1_file,
            "-right", r2_file,
            "-t", str(self.threads),
            "-outdir", str(output_dir),
            "-name", sample_name,
            "-length", str(min_length)
        ]
        
        # Add single-end reads if provided
        if single_end:
            cmd.extend(["-sef", single_end])
            logger.info("Including single-end reads in assembly")
        
        # Add stranded option if specified
        if stranded:
            cmd.append("-stranded")
            logger.info("Using stranded assembly for bulk RNA-seq data")
        
        # Add reference if provided
        if reference:
            cmd.extend(["-ref", reference])
            logger.info(f"Using reference-guided assembly with {reference}")
        
        # Add memory settings (convert from "16G" to "16")
        if self.memory:
            memory_value = self.memory.replace('G', '').replace('g', '').replace('M', '').replace('m', '')
            cmd.extend(["-mem", memory_value])
        
        logger.info(f"RNA-Bloom command: {' '.join(cmd)}")
        
        # Run RNA-Bloom
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(output_dir))
            
            if result.returncode != 0:
                logger.error(f"RNA-Bloom failed: {result.stderr}")
                raise RuntimeError(f"RNA-Bloom assembly failed: {result.stderr}")
            
            logger.info("RNA-Bloom bulk RNA-seq assembly completed successfully")
            
            # Find and move output files to clean structure
            transcripts_file = self._organize_output_files(output_dir, sample_name, "transcripts.fa")
            transcripts_short_file = self._organize_output_files(output_dir, sample_name, "transcripts.short.fa")
            transcripts_nr_file = self._organize_output_files(output_dir, sample_name, "transcripts.nr.fa")
            
            # Return output file paths
            return {
                "transcripts": transcripts_file,
                "transcripts_short": transcripts_short_file,
                "transcripts_nr": transcripts_nr_file,
                "assembly_log": str(output_dir / "rnabloom.log")
            }
            
        except Exception as e:
            logger.error(f"RNA-Bloom assembly failed: {e}")
            raise
    
    def assemble_long_reads(
        self,
        long_reads: str,
        output_dir: Path,
        sample_name: str = "sample",
        min_length: int = 200,
        stranded: bool = False,
        pacbio: bool = False,
        short_reads_forward: Optional[str] = None,
        short_reads_reverse: Optional[str] = None
    ) -> Dict[str, str]:
        """
        Assemble long-read RNA-seq data using RNA-Bloom.
        
        Args:
            long_reads: Path to long reads (ONT or PacBio)
            output_dir: Output directory for assembly results
            sample_name: Sample name for output files
            min_length: Minimum transcript length
            stranded: Whether reads are strand-specific (for direct RNA)
            pacbio: Whether reads are PacBio (vs ONT)
            short_reads_forward: Optional short reads for polishing
            short_reads_reverse: Optional reverse short reads for polishing
            
        Returns:
            Dictionary with assembly output file paths
        """
        logger.info(f"Assembling long-read RNA-seq data for {sample_name}")
        
        # Create output directory
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Prepare RNA-Bloom command for long reads
        cmd = [
            "rnabloom",
            "-long", long_reads,
            "-t", str(self.threads),
            "-outdir", str(output_dir),
            "-name", sample_name,
            "-length", str(min_length)
        ]
        
        # Add PacBio flag if specified
        if pacbio:
            cmd.append("-lrpb")
            logger.info("Using PacBio-specific parameters")
        
        # Add stranded option for direct RNA sequencing
        if stranded:
            cmd.append("-stranded")
            logger.info("Using stranded assembly for direct RNA sequencing")
        
        # Add short reads for polishing if provided
        if short_reads_forward:
            cmd.extend(["-sef", short_reads_forward])
            if short_reads_reverse:
                cmd.extend(["-ser", short_reads_reverse])
            logger.info("Including short reads for polishing")
        
        # Add memory settings (convert from "16G" to "16")
        if self.memory:
            memory_value = self.memory.replace('G', '').replace('g', '').replace('M', '').replace('m', '')
            cmd.extend(["-mem", memory_value])
        
        logger.info(f"RNA-Bloom command: {' '.join(cmd)}")
        
        # Run RNA-Bloom
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(output_dir))
            
            if result.returncode != 0:
                logger.error(f"RNA-Bloom failed: {result.stderr}")
                raise RuntimeError(f"RNA-Bloom assembly failed: {result.stderr}")
            
            logger.info("RNA-Bloom long-read assembly completed successfully")
            
            # Find and move output files to clean structure
            transcripts_file = self._organize_output_files(output_dir, sample_name, "transcripts.fa")
            transcripts_short_file = self._organize_output_files(output_dir, sample_name, "transcripts.short.fa")
            
            # Return output file paths
            return {
                "transcripts": transcripts_file,
                "transcripts_short": transcripts_short_file,
                "assembly_log": str(output_dir / "rnabloom.log")
            }
            
        except Exception as e:
            logger.error(f"RNA-Bloom assembly failed: {e}")
            raise
    
    def assemble_hybrid(
        self,
        short_r1: str,
        short_r2: str,
        long_reads: str,
        output_dir: Path,
        sample_name: str = "sample",
        min_length: int = 200,
        stranded: bool = False,
        pacbio: bool = False
    ) -> Dict[str, str]:
        """
        Assemble hybrid RNA-seq data (short + long reads) using RNA-Bloom.
        
        Args:
            short_r1: Path to short R1 reads
            short_r2: Path to short R2 reads
            long_reads: Path to long reads
            output_dir: Output directory for assembly results
            sample_name: Sample name for output files
            min_length: Minimum transcript length
            stranded: Whether reads are strand-specific
            pacbio: Whether long reads are PacBio (vs ONT)
            
        Returns:
            Dictionary with assembly output file paths
        """
        logger.info(f"Assembling hybrid RNA-seq data for {sample_name}")
        
        # Create output directory
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Prepare RNA-Bloom command for hybrid assembly
        cmd = [
            "rnabloom",
            "-left", short_r1,
            "-right", short_r2,
            "-long", long_reads,
            "-t", str(self.threads),
            "-outdir", str(output_dir),
            "-name", sample_name,
            "-length", str(min_length)
        ]
        
        # Add PacBio flag if specified
        if pacbio:
            cmd.append("-lrpb")
            logger.info("Using PacBio-specific parameters for long reads")
        
        # Add stranded option if specified
        if stranded:
            cmd.append("-stranded")
            logger.info("Using stranded assembly for hybrid data")
        
        # Add memory settings (convert from "16G" to "16")
        if self.memory:
            memory_value = self.memory.replace('G', '').replace('g', '').replace('M', '').replace('m', '')
            cmd.extend(["-mem", memory_value])
        
        logger.info(f"RNA-Bloom command: {' '.join(cmd)}")
        
        # Run RNA-Bloom
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(output_dir))
            
            if result.returncode != 0:
                logger.error(f"RNA-Bloom failed: {result.stderr}")
                raise RuntimeError(f"RNA-Bloom assembly failed: {result.stderr}")
            
            logger.info("RNA-Bloom hybrid assembly completed successfully")
            
            # Find and move output files to clean structure
            transcripts_file = self._organize_output_files(output_dir, sample_name, "transcripts.fa")
            transcripts_short_file = self._organize_output_files(output_dir, sample_name, "transcripts.short.fa")
            transcripts_nr_file = self._organize_output_files(output_dir, sample_name, "transcripts.nr.fa")
            
            # Return output file paths
            return {
                "transcripts": transcripts_file,
                "transcripts_short": transcripts_short_file,
                "transcripts_nr": transcripts_nr_file,
                "assembly_log": str(output_dir / "rnabloom.log")
            }
            
        except Exception as e:
            logger.error(f"RNA-Bloom assembly failed: {e}")
            raise
    
    def create_cell_sample_sheet(
        self,
        cell_data: Dict[str, Dict[str, str]],
        output_file: Path
    ) -> str:
        """
        Create a sample sheet for RNA-Bloom cell-level assembly.
        
        Args:
            cell_data: Dictionary mapping cell IDs to read file paths
            output_file: Path to output sample sheet
            
        Returns:
            Path to created sample sheet
        """
        logger.info(f"Creating RNA-Bloom sample sheet for {len(cell_data)} cells")
        
        with open(output_file, 'w') as f:
            # Write header
            f.write("#name left right\n")
            
            # Write cell data
            for cell_id, paths in cell_data.items():
                r1_path = paths.get('r1')
                r2_path = paths.get('r2')
                if r1_path and r2_path:
                    f.write(f"{cell_id}\t{r1_path}\t{r2_path}\n")
                else:
                    logger.warning(f"Skipping cell {cell_id}: missing R1 or R2 files")
        
        logger.info(f"Sample sheet created: {output_file}")
        return str(output_file)
    
    def assemble_multiple_cells(
        self,
        sample_sheet: str,
        output_dir: Path,
        min_length: int = 200,
        stranded: bool = False,
        reference: Optional[str] = None
    ) -> Dict[str, str]:
        """
        Assemble multiple cells using RNA-Bloom with sample sheet.
        
        Args:
            sample_sheet: Path to sample sheet file
            output_dir: Output directory for assembly results
            min_length: Minimum transcript length
            stranded: Whether reads are strand-specific
            reference: Optional reference transcriptome for guided assembly
            
        Returns:
            Dictionary with assembly output file paths
        """
        logger.info("Assembling multiple cells with RNA-Bloom")
        
        # Create output directory
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Prepare RNA-Bloom command for multiple cells
        cmd = [
            "rnabloom",
            "-pool", sample_sheet,
            "-t", str(self.threads),
            "-outdir", str(output_dir),
            "-length", str(min_length)
        ]
        
        # Add stranded option if specified
        if stranded:
            cmd.append("-stranded")
            logger.info("Using stranded assembly for multiple cells")
        
        # Add reference if provided
        if reference:
            cmd.extend(["-ref", reference])
            logger.info(f"Using reference-guided assembly with {reference}")
        
        # Add memory settings (convert from "16G" to "16")
        if self.memory:
            memory_value = self.memory.replace('G', '').replace('g', '').replace('M', '').replace('m', '')
            cmd.extend(["-mem", memory_value])
        
        logger.info(f"RNA-Bloom command: {' '.join(cmd)}")
        
        # Run RNA-Bloom
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(output_dir))
            
            if result.returncode != 0:
                logger.error(f"RNA-Bloom failed: {result.stderr}")
                raise RuntimeError(f"RNA-Bloom assembly failed: {result.stderr}")
            
            logger.info("RNA-Bloom multi-cell assembly completed successfully")
            
            # Return output file paths
            return {
                "transcripts": str(output_dir / "rnabloom.transcripts.fa"),
                "transcripts_short": str(output_dir / "rnabloom.transcripts.short.fa"),
                "transcripts_nr": str(output_dir / "rnabloom.transcripts.nr.fa"),
                "assembly_log": str(output_dir / "rnabloom.log")
            }
            
        except Exception as e:
            logger.error(f"RNA-Bloom assembly failed: {e}")
            raise
    
    def _organize_output_files(self, output_dir: Path, sample_name: str, filename: str) -> str:
        """
        Find RNA-Bloom output files and organize them in a clean structure.
        
        Args:
            output_dir: Base output directory
            sample_name: Sample name
            filename: Name of the file to find (e.g., "transcripts.fa")
            
        Returns:
            Path to the organized file
        """
        full_filename = f"{sample_name}.{filename}"
        target_path = output_dir / full_filename
        
        # First try the direct path
        if target_path.exists():
            return str(target_path)
        
        # Search recursively for the file
        found_file = None
        for file_path in output_dir.rglob(full_filename):
            if file_path.is_file():
                found_file = file_path
                break
        
        if found_file and found_file != target_path:
            # Move file to clean location
            logger.info(f"Moving {found_file} to {target_path}")
            shutil.move(str(found_file), str(target_path))
        
        return str(target_path)
    
    def _organize_pooled_output_files(self, output_dir: Path, filename: str) -> str:
        """
        Organize output files from RNA-Bloom pooled mode.
        
        Args:
            output_dir: Base output directory
            filename: Name of the file to find (e.g., "transcripts.fa")
            
        Returns:
            Path to the organized file
        """
        target_path = output_dir / filename
        
        # Look for files in subdirectories (cell_1/, cell_2/, etc.)
        found_files = []
        for subdir in output_dir.iterdir():
            if subdir.is_dir() and not subdir.name.startswith('.'):
                file_path = subdir / filename
                if file_path.exists():
                    found_files.append(file_path)
        
        if found_files:
            # Combine all files from different cells
            logger.info(f"Found {len(found_files)} {filename} files from pooled assembly")
            self._combine_pooled_files(found_files, target_path)
        else:
            logger.warning(f"Could not find {filename} in pooled assembly output")
        
        return str(target_path)
    
    def _combine_pooled_files(self, source_files: List[Path], target_file: Path) -> None:
        """
        Combine files from multiple cells into a single file.
        
        Args:
            source_files: List of source files to combine
            target_file: Target file to write combined content
        """
        from Bio import SeqIO
        
        all_sequences = []
        for source_file in source_files:
            for record in SeqIO.parse(source_file, "fasta"):
                # Add cell information to sequence ID
                cell_name = source_file.parent.name
                record.id = f"{cell_name}_{record.id}"
                all_sequences.append(record)
        
        # Write combined file
        SeqIO.write(all_sequences, target_file, "fasta")
        logger.info(f"Combined {len(all_sequences)} sequences from {len(source_files)} cells into {target_file}")
    
    def cleanup(self) -> None:
        """Clean up temporary files."""
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
            logger.info("Cleaned up RNA-Bloom temporary files")
