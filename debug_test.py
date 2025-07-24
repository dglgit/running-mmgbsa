#!/usr/bin/env python3
"""
Debug script to test file operations and directory creation.
"""

import os
import shutil
from pathlib import Path

def test_file_operations():
    """Test basic file operations."""
    print("Testing file operations...")
    
    # Test with the actual file paths
    protein_file = Path("3htb/3htb_cleaned_protonated_fixed.pdb")
    ligand_file = Path("3htb/JZ4.mol2")
    output_dir = Path("debug_test_results")
    
    print(f"Current working directory: {os.getcwd()}")
    print(f"Protein file: {protein_file} (exists: {protein_file.exists()})")
    print(f"Ligand file: {ligand_file} (exists: {ligand_file.exists()})")
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    work_dir = output_dir / "work"
    work_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Created work directory: {work_dir}")
    
    # Test copying files
    if protein_file.exists():
        protein_work = work_dir / "protein_input.pdb"
        try:
            shutil.copy2(protein_file, protein_work)
            print(f"✓ Successfully copied protein file to: {protein_work}")
        except Exception as e:
            print(f"✗ Failed to copy protein file: {e}")
    else:
        print("✗ Protein file not found")
    
    if ligand_file.exists():
        ligand_work = work_dir / f"ligand_input{ligand_file.suffix}"
        try:
            shutil.copy2(ligand_file, ligand_work)
            print(f"✓ Successfully copied ligand file to: {ligand_work}")
        except Exception as e:
            print(f"✗ Failed to copy ligand file: {e}")
    else:
        print("✗ Ligand file not found")
    
    # List contents of work directory
    print(f"\nContents of work directory:")
    for item in work_dir.iterdir():
        print(f"  - {item.name}")

if __name__ == "__main__":
    test_file_operations() 