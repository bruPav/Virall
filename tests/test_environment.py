#!/usr/bin/env python3
"""
Test environment setup for single-cell feature development.
This script helps ensure the current functionality still works before adding new features.
"""

import os
import sys
import tempfile
import shutil
from pathlib import Path

# Add the virall package to the path
sys.path.insert(0, str(Path(__file__).parent.parent))

def test_current_functionality():
    """Test that current functionality still works."""
    print("Testing current functionality...")
    
    try:
        # Test imports
        from virall.core.assembler import ViralAssembler
        from virall.core.preprocessor import Preprocessor
        from virall.core.viral_identifier import ViralIdentifier
        from virall.core.validator import AssemblyValidator
        from virall.core.gene_predictor import ViralGenePredictor
        print("✓ All core modules import successfully")
        
        # Test basic initialization
        with tempfile.TemporaryDirectory() as temp_dir:
            assembler = ViralAssembler(output_dir=temp_dir, threads=2, memory="4G")
            print("✓ ViralAssembler initializes correctly")
            
            preprocessor = Preprocessor(threads=2)
            print("✓ Preprocessor initializes correctly")
            
            validator = AssemblyValidator(threads=2)
            print("✓ AssemblyValidator initializes correctly")
            
        print("✓ All basic functionality tests passed")
        return True
        
    except Exception as e:
        print(f"✗ Test failed: {e}")
        return False

def create_test_data():
    """Create minimal test data for single-cell testing."""
    test_data_dir = Path(__file__).parent / "data" / "sample_sc_data"
    test_data_dir.mkdir(parents=True, exist_ok=True)
    
    # Create dummy FASTQ files for testing
    dummy_r1 = test_data_dir / "sample_R1.fastq"
    dummy_r2 = test_data_dir / "sample_R2.fastq"
    dummy_i1 = test_data_dir / "sample_I1.fastq"
    
    # Create minimal FASTQ content
    fastq_content = """@read1
ATCGATCGATCG
+
IIIIIIIIIIII
@read2
GCTAGCTAGCTA
+
IIIIIIIIIIII
"""
    
    with open(dummy_r1, 'w') as f:
        f.write(fastq_content)
    
    with open(dummy_r2, 'w') as f:
        f.write(fastq_content)
    
    with open(dummy_i1, 'w') as f:
        f.write(fastq_content)
    
    print(f"✓ Created test data in {test_data_dir}")
    return test_data_dir

def check_dependencies():
    """Check if required dependencies are available."""
    print("Checking dependencies...")
    
    required_tools = [
        'fastqc', 'trimmomatic', 'spades.py', 'flye', 
        'bwa', 'samtools', 'minimap2'
    ]
    
    missing_tools = []
    for tool in required_tools:
        if shutil.which(tool) is None:
            missing_tools.append(tool)
    
    if missing_tools:
        print(f"⚠ Missing tools: {', '.join(missing_tools)}")
        print("Some features may not work properly")
    else:
        print("✓ All required tools are available")
    
    return len(missing_tools) == 0

if __name__ == "__main__":
    print("=" * 50)
    print("VIRALL TEST ENVIRONMENT SETUP")
    print("=" * 50)
    
    # Check dependencies
    deps_ok = check_dependencies()
    
    # Test current functionality
    func_ok = test_current_functionality()
    
    # Create test data
    test_data_dir = create_test_data()
    
    print("\n" + "=" * 50)
    print("TEST SUMMARY")
    print("=" * 50)
    print(f"Dependencies: {'✓ OK' if deps_ok else '⚠ Some missing'}")
    print(f"Functionality: {'✓ OK' if func_ok else '✗ FAILED'}")
    print(f"Test data: ✓ Created in {test_data_dir}")
    
    if func_ok:
        print("\n✓ Environment is ready for single-cell feature development!")
    else:
        print("\n✗ Please fix issues before proceeding with development")
        sys.exit(1)
