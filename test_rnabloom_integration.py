#!/usr/bin/env python3
"""
Test script to verify RNA-Bloom integration with Virall.
"""

import sys
from pathlib import Path

# Add the virall package to the path
sys.path.insert(0, str(Path(__file__).parent))

def test_rnabloom_import():
    """Test if RNA-Bloom assembler can be imported."""
    try:
        from virall.core.rnabloom_assembler import RNABloomAssembler
        print("✅ RNA-Bloom assembler import successful")
        return True
    except ImportError as e:
        print(f"❌ RNA-Bloom assembler import failed: {e}")
        return False

def test_rnabloom_initialization():
    """Test if RNA-Bloom assembler can be initialized."""
    try:
        from virall.core.rnabloom_assembler import RNABloomAssembler
        assembler = RNABloomAssembler(threads=4, memory="8G")
        print("✅ RNA-Bloom assembler initialization successful")
        return True
    except Exception as e:
        print(f"❌ RNA-Bloom assembler initialization failed: {e}")
        return False

def test_viral_assembler_integration():
    """Test if ViralAssembler can use RNA-Bloom."""
    try:
        from virall.core.assembler import ViralAssembler
        
        # Test with RNA mode enabled
        assembler = ViralAssembler(
            output_dir="test_output",
            threads=4,
            memory="8G",
            rna_mode=True
        )
        
        if hasattr(assembler, 'rnabloom_assembler'):
            if assembler.rnabloom_assembler is not None:
                print("✅ ViralAssembler RNA-Bloom integration successful")
                return True
            else:
                print("⚠️  RNA-Bloom not available in ViralAssembler")
                return False
        else:
            print("❌ ViralAssembler does not have RNA-Bloom integration")
            return False
            
    except Exception as e:
        print(f"❌ ViralAssembler RNA-Bloom integration failed: {e}")
        return False

def main():
    """Run all tests."""
    print("Testing RNA-Bloom integration with Virall...")
    print("=" * 50)
    
    tests = [
        test_rnabloom_import,
        test_rnabloom_initialization,
        test_viral_assembler_integration
    ]
    
    passed = 0
    total = len(tests)
    
    for test in tests:
        if test():
            passed += 1
        print()
    
    print("=" * 50)
    print(f"Tests passed: {passed}/{total}")
    
    if passed == total:
        print("🎉 All tests passed! RNA-Bloom integration is working.")
        return 0
    else:
        print("⚠️  Some tests failed. Check the output above.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
