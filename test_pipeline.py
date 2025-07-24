#!/usr/bin/env python3
"""
Test script for the MMGBSA pipeline.

This script tests the pipeline components without running full calculations.
"""

import os
import sys
from pathlib import Path

# Add the current directory to Python path
sys.path.insert(0, str(Path(__file__).parent))

def test_imports():
    """Test that all required modules can be imported."""
    print("Testing imports...")
    
    try:
        import openmm
        print(f"✓ OpenMM version: {openmm.version.version}")
    except ImportError as e:
        print(f"✗ OpenMM import failed: {e}")
        return False
    
    try:
        import numpy
        print(f"✓ NumPy version: {numpy.__version__}")
    except ImportError as e:
        print(f"✗ NumPy import failed: {e}")
        return False
    
    try:
        import yaml
        print(f"✓ PyYAML imported successfully")
    except ImportError as e:
        print(f"✗ PyYAML import failed: {e}")
        return False
    
    try:
        from mmgbsa_pipeline import MMGBSAPipeline
        print("✓ MMGBSAPipeline imported successfully")
    except ImportError as e:
        print(f"✗ MMGBSAPipeline import failed: {e}")
        return False
    
    return True

def test_external_tools():
    """Test that external tools are available."""
    print("\nTesting external tools...")
    
    tools = ['tleap', 'cpptraj', 'MMPBSA.py']
    all_available = True
    
    for tool in tools:
        try:
            result = os.system(f"which {tool} > /dev/null 2>&1")
            if result == 0:
                print(f"✓ {tool} found")
            else:
                print(f"✗ {tool} not found in PATH")
                all_available = False
        except Exception as e:
            print(f"✗ Error checking {tool}: {e}")
            all_available = False
    
    # Test optional tools
    optional_tools = ['pdb4amber', 'reduce', 'obabel']
    print("\nTesting optional tools (for protein/ligand processing):")
    for tool in optional_tools:
        try:
            result = os.system(f"which {tool} > /dev/null 2>&1")
            if result == 0:
                print(f"✓ {tool} found")
            else:
                print(f"⚠ {tool} not found (will use fallback methods)")
        except Exception as e:
            print(f"⚠ Error checking {tool}: {e}")
    
    return all_available

def test_input_files():
    """Test that example input files exist."""
    print("\nTesting input files...")
    
    expected_files = [
        "3htb/3htb_cleaned_protonated_fixed.pdb",
        "3htb/JZ4.mol2",
        "3htb/JZ4.pdb"
    ]
    
    all_exist = True
    for file_path in expected_files:
        if os.path.exists(file_path):
            print(f"✓ {file_path} found")
        else:
            print(f"✗ {file_path} not found")
            all_exist = False
    
    return all_exist

def test_pipeline_initialization():
    """Test pipeline initialization with existing files."""
    print("\nTesting pipeline initialization...")
    
    try:
        from mmgbsa_pipeline import MMGBSAPipeline
        
        # Test with existing files
        if os.path.exists("3htb/3htb_cleaned_protonated_fixed.pdb") and os.path.exists("3htb/JZ4.mol2"):
            pipeline = MMGBSAPipeline(
                protein_file="3htb/3htb_cleaned_protonated_fixed.pdb",
                ligand_file="3htb/JZ4.mol2",
                output_dir="test_output",
                method="mmgbsa"
            )
            print("✓ Pipeline initialized successfully")
            return True
        else:
            print("✗ Required input files not found for initialization test")
            return False
            
    except Exception as e:
        print(f"✗ Pipeline initialization failed: {e}")
        return False

def test_simplified_interface():
    """Test that the simplified interface works correctly."""
    print("\nTesting simplified interface...")
    
    try:
        from mmgbsa_pipeline import MMGBSAPipeline
        
        # Test that we can create a pipeline with just protein and ligand
        if os.path.exists("3htb/3htb_cleaned_protonated_fixed.pdb") and os.path.exists("3htb/JZ4.mol2"):
            pipeline = MMGBSAPipeline(
                protein_file="3htb/3htb_cleaned_protonated_fixed.pdb",
                ligand_file="3htb/JZ4.mol2",
                output_dir="test_output"
            )
            print("✓ Simplified interface works (protein + ligand only)")
            print("✓ Complex will be created automatically")
            return True
        else:
            print("✗ Required input files not found for interface test")
            return False
            
    except Exception as e:
        print(f"✗ Simplified interface test failed: {e}")
        return False

def main():
    """Run all tests."""
    print("MMGBSA Pipeline Test Suite")
    print("=" * 40)
    print("Testing simplified interface (protein + ligand only)")
    print("Complex creation, charging, and protonation handled automatically")
    print()
    
    tests = [
        ("Import Test", test_imports),
        ("External Tools Test", test_external_tools),
        ("Input Files Test", test_input_files),
        ("Pipeline Initialization Test", test_pipeline_initialization),
        ("Simplified Interface Test", test_simplified_interface)
    ]
    
    results = []
    for test_name, test_func in tests:
        print(f"\n{test_name}:")
        print("-" * len(test_name))
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"✗ Test failed with exception: {e}")
            results.append((test_name, False))
    
    # Summary
    print("\n" + "=" * 40)
    print("TEST SUMMARY")
    print("=" * 40)
    
    passed = 0
    total = len(results)
    
    for test_name, result in results:
        status = "PASS" if result else "FAIL"
        print(f"{test_name}: {status}")
        if result:
            passed += 1
    
    print(f"\nOverall: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All tests passed! The pipeline should work correctly.")
        print("\nYou can now run the pipeline with just:")
        print("python mmgbsa_pipeline.py --protein protein.pdb --ligand ligand.mol2 --output_dir results/")
    else:
        print("⚠️  Some tests failed. Please check the issues above.")
        print("\nCommon solutions:")
        print("1. Install missing Python packages: pip install -r requirements.txt")
        print("2. Install AmberTools and ensure it's in your PATH")
        print("3. Ensure input files are present in the 3htb/ directory")
        print("4. Optional: Install pdb4amber, reduce, and obabel for better file processing")
    
    return passed == total

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1) 