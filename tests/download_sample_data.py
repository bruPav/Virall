#!/usr/bin/env python3
"""
Helper script to download sample single-cell data for testing.
This script provides easy access to commonly used test datasets.
"""

import os
import sys
import subprocess
import requests
from pathlib import Path

def download_10x_pbmc_1k():
    """Download 10X Genomics PBMC 1k dataset (smallest for testing)."""
    print("📥 Downloading 10X Genomics PBMC 1k dataset...")
    
    # This is a placeholder - actual download would depend on the specific dataset
    # For now, we'll provide instructions
    print("\n🔗 Manual Download Instructions:")
    print("1. Go to: https://www.10xgenomics.com/resources/datasets")
    print("2. Search for 'PBMC' or 'peripheral blood mononuclear cells'")
    print("3. Look for the smallest dataset (usually 1k cells)")
    print("4. Download the FASTQ files")
    print("5. Place them in: tests/data/real_sc_data/10x_genomics/")
    
    print("\n📁 Expected files:")
    print("   - Sample1_S1_L001_R1_001.fastq.gz")
    print("   - Sample1_S1_L001_R2_001.fastq.gz") 
    print("   - Sample1_S1_L001_I1_001.fastq.gz")

def create_sample_data_structure():
    """Create the proper directory structure for sample data."""
    base_dir = Path(__file__).parent / "data" / "real_sc_data"
    
    directories = [
        "10x_genomics",
        "other_formats",
        "10x_genomics/sample1",
        "10x_genomics/sample2"
    ]
    
    for dir_path in directories:
        full_path = base_dir / dir_path
        full_path.mkdir(parents=True, exist_ok=True)
        print(f"✅ Created directory: {full_path}")

def check_disk_space():
    """Check available disk space."""
    import shutil
    
    total, used, free = shutil.disk_usage("/")
    
    print(f"💾 Disk Space Check:")
    print(f"   Total: {total // (1024**3)} GB")
    print(f"   Used:  {used // (1024**3)} GB") 
    print(f"   Free:  {free // (1024**3)} GB")
    
    # Recommend at least 5GB free for testing
    if free < 5 * (1024**3):
        print("⚠️  Warning: Less than 5GB free space. Consider freeing up space.")
    else:
        print("✅ Sufficient disk space available")

def main():
    print("=" * 60)
    print("SINGLE-CELL DATA DOWNLOAD HELPER")
    print("=" * 60)
    
    # Check disk space
    check_disk_space()
    
    # Create directory structure
    print("\n📁 Creating directory structure...")
    create_sample_data_structure()
    
    # Provide download instructions
    print("\n📥 Download Instructions:")
    download_10x_pbmc_1k()
    
    print("\n" + "=" * 60)
    print("NEXT STEPS:")
    print("=" * 60)
    print("1. Download the sample data as instructed above")
    print("2. Place files in: tests/data/real_sc_data/10x_genomics/")
    print("3. Run: python tests/test_real_single_cell.py")
    print("4. Check the results in: tests/real_sc_test_output/")

if __name__ == "__main__":
    main()
