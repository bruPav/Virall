#!/usr/bin/env python3
"""
Test script for real single-cell data processing.
This script helps you test the single-cell functionality with real datasets.
"""

import os
import sys
import subprocess
from pathlib import Path

# Add the virall package to the path
sys.path.insert(0, str(Path(__file__).parent.parent))

def find_test_data():
    """Find available test data in the real_sc_data directory."""
    test_data_dir = Path(__file__).parent / "data" / "real_sc_data"
    
    print("Searching for test data...")
    print(f"Looking in: {test_data_dir}")
    
    if not test_data_dir.exists():
        print("❌ Test data directory not found")
        return None
    
    # Look for 10X Genomics data
    formats = {
        "10x_genomics": test_data_dir / "10x_genomics",
        "other_formats": test_data_dir / "other_formats"
    }
    
    found_data = {}
    
    for format_name, format_dir in formats.items():
        if format_dir.exists():
            print(f"\n📁 {format_name.upper()} data:")
            
            # Look for common file patterns
            patterns = {
                "R1": "*R1*.fastq*",
                "R2": "*R2*.fastq*", 
                "I1": "*I1*.fastq*",
                "I2": "*I2*.fastq*"
            }
            
            format_files = {}
            for file_type, pattern in patterns.items():
                files = list(format_dir.glob(pattern))
                if files:
                    format_files[file_type] = files[0]  # Take first match
                    print(f"  ✅ {file_type}: {files[0].name}")
                else:
                    print(f"  ❌ {file_type}: Not found")
            
            if format_files:
                found_data[format_name] = format_files
    
    return found_data

def test_with_real_data(data_files, min_cells=100):
    """Test single-cell functionality with real data."""
    print(f"\n🧪 Testing with real data (min_cells={min_cells})...")
    
    # Check if we have the required files
    required_files = ["R1", "R2"]
    missing_files = [f for f in required_files if f not in data_files]
    
    if missing_files:
        print(f"❌ Missing required files: {missing_files}")
        return False
    
    # Create output directory
    output_dir = Path(__file__).parent / "real_sc_test_output"
    output_dir.mkdir(exist_ok=True)
    
    # Build command
    cmd = [
        "python", "-m", "virall.cli", "assemble",
        "--single-cell",
        "--short-reads-1", str(data_files["R1"]),
        "--short-reads-2", str(data_files["R2"]),
        "--output-dir", str(output_dir),
        "--min-cells", str(min_cells),
        "--threads", "4",
        "--memory", "8G"
    ]
    
    # Add index file if available
    if "I1" in data_files:
        cmd.extend(["--index-reads", str(data_files["I1"])])
    
    print(f"Running command: {' '.join(cmd)}")
    
    try:
        # Run the command
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)  # 5 min timeout
        
        print(f"\n📊 Results:")
        print(f"Return code: {result.returncode}")
        
        if result.stdout:
            print(f"STDOUT:\n{result.stdout}")
        
        if result.stderr:
            print(f"STDERR:\n{result.stderr}")
        
        # Check if output files were created
        expected_files = [
            "single_cell_processed/processed_R1.fastq",
            "single_cell_processed/processed_R2.fastq", 
            "single_cell_processed/cell_assignments.tsv",
            "single_cell_processed/cell_metadata.tsv",
            "single_cell_processed/single_cell_qc_report.txt"
        ]
        
        print(f"\n📁 Checking output files:")
        for expected_file in expected_files:
            file_path = output_dir / expected_file
            if file_path.exists():
                size = file_path.stat().st_size
                print(f"  ✅ {expected_file} ({size:,} bytes)")
            else:
                print(f"  ❌ {expected_file} (not found)")
        
        return result.returncode == 0
        
    except subprocess.TimeoutExpired:
        print("❌ Test timed out after 5 minutes")
        return False
    except Exception as e:
        print(f"❌ Test failed with error: {e}")
        return False

def download_sample_data():
    """Provide instructions for downloading sample data."""
    print("\n📥 To download sample single-cell data:")
    print("\n1. 10X Genomics Public Datasets:")
    print("   https://www.10xgenomics.com/resources/datasets")
    print("   - Look for 'PBMC' datasets")
    print("   - Download the smallest available dataset first")
    
    print("\n2. Single Cell Portal:")
    print("   https://singlecell.broadinstitute.org/")
    print("   - Search for 'PBMC' or 'peripheral blood'")
    print("   - Download FASTQ files")
    
    print("\n3. GEO (Gene Expression Omnibus):")
    print("   https://www.ncbi.nlm.nih.gov/geo/")
    print("   - Search for 'single cell RNA-seq'")
    print("   - Look for datasets with FASTQ files")
    
    print("\n4. Recommended for testing:")
    print("   - 10X Genomics PBMC 1k cells (smallest)")
    print("   - File naming: Sample1_S1_L001_R1_001.fastq.gz")
    print("   - Place in: tests/data/real_sc_data/10x_genomics/")

def main():
    print("=" * 60)
    print("REAL SINGLE-CELL DATA TESTING")
    print("=" * 60)
    
    # Find test data
    data_files = find_test_data()
    
    if not data_files:
        print("\n❌ No test data found!")
        download_sample_data()
        return
    
    # Test with each found dataset
    for format_name, files in data_files.items():
        print(f"\n{'='*40}")
        print(f"TESTING {format_name.upper()} DATA")
        print(f"{'='*40}")
        
        success = test_with_real_data(files, min_cells=10)  # Lower threshold for testing
        
        if success:
            print(f"✅ {format_name} test PASSED")
        else:
            print(f"❌ {format_name} test FAILED")
    
    print(f"\n{'='*60}")
    print("TESTING COMPLETE")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()
