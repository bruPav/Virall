#!/usr/bin/env python3
"""
Test script for single-cell preprocessing functionality.
"""

import os
import sys
import tempfile
from pathlib import Path

# Add the virall package to the path
sys.path.insert(0, str(Path(__file__).parent.parent))

def test_single_cell_preprocessor():
    """Test the single-cell preprocessor functionality."""
    print("Testing SingleCellPreprocessor...")
    
    try:
        from virall.core.single_cell_preprocessor import SingleCellPreprocessor
        
        # Initialize preprocessor
        preprocessor = SingleCellPreprocessor(threads=2)
        print("✓ SingleCellPreprocessor initializes correctly")
        
        # Test with dummy data
        test_data_dir = Path(__file__).parent / "data" / "sample_sc_data"
        
        if not test_data_dir.exists():
            print("✗ Test data directory not found")
            return False
        
        r1_file = test_data_dir / "sample_R1.fastq"
        r2_file = test_data_dir / "sample_R2.fastq"
        
        if not (r1_file.exists() and r2_file.exists()):
            print("✗ Test data files not found")
            return False
        
        # Test barcode extraction
        print("Testing barcode extraction...")
        barcodes = preprocessor.extract_cell_barcodes(str(r1_file))
        print(f"✓ Extracted {len(barcodes)} barcodes")
        
        # Test processing with temporary output
        with tempfile.TemporaryDirectory() as temp_dir:
            print("Testing 10X data processing...")
            results = preprocessor.process_10x_data(
                r1_file=str(r1_file),
                r2_file=str(r2_file),
                sample_id="test_sample",
                output_dir=temp_dir,
                min_cells=1  # Lower threshold for test data
            )
            
            print(f"✓ Processed {results['num_cells']} cells")
            print(f"✓ Output directory: {results['output_dir']}")
        
        print("✓ SingleCellPreprocessor tests passed")
        return True
        
    except Exception as e:
        print(f"✗ SingleCellPreprocessor test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    print("=" * 50)
    print("SINGLE-CELL PREPROCESSOR TEST")
    print("=" * 50)
    
    success = test_single_cell_preprocessor()
    
    print("\n" + "=" * 50)
    print("TEST RESULT")
    print("=" * 50)
    print(f"Status: {'✓ PASSED' if success else '✗ FAILED'}")
    
    if not success:
        sys.exit(1)
