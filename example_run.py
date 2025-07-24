#!/usr/bin/env python3
"""
Example script demonstrating how to use the MMGBSA pipeline programmatically.

This script shows how to run MMGBSA calculations on a protein-ligand complex
using the refactored pipeline.
"""

import os
import sys
from pathlib import Path

# Add the current directory to Python path to import the pipeline
sys.path.insert(0, str(Path(__file__).parent))

from mmgbsa_pipeline import MMGBSAPipeline

def run_example():
    """Run an example MMGBSA calculation."""
    
    # Example file paths (modify these to match your actual files)
    protein_file = "3htb/3htb_cleaned_protonated_fixed.pdb"  # Your protein PDB
    ligand_file = "3htb/JZ4.mol2"  # Your ligand file
    output_dir = "3htb_results"
    
    # Check if input files exist
    required_files = [protein_file, ligand_file]
    for file_path in required_files:
        if not os.path.exists(file_path):
            print(f"Warning: Input file {file_path} not found.")
            print("Please ensure all input files are present before running.")
            return
    
    print("Starting MMGBSA calculation...")
    print(f"Protein file: {protein_file}")
    print(f"Ligand file: {ligand_file}")
    print(f"Output directory: {output_dir}")
    print("Note: Complex will be created automatically with proper charging and protonation")
    
    # Create pipeline with default parameters
    pipeline = MMGBSAPipeline(
        protein_file=protein_file,
        ligand_file=ligand_file,
        output_dir=output_dir,
        method="mmgbsa",  # or "mmpbsa"
        temperature=300,
        equilibration_steps=5000,  # Reduced for example
        production_steps=100000,   # Reduced for example
        trajectory_interval=500,
        salt_concentration=0.15,
        start_frame=1,
        end_frame=100,  # Reduced for example
        frame_interval=1
    )
    
    try:
        # Run the pipeline
        pipeline.run()
        
        print("\n" + "="*50)
        print("CALCULATION COMPLETED SUCCESSFULLY!")
        print("="*50)
        
        # Display results summary
        results_file = os.path.join(output_dir, "mmgbsa_results.dat")
        if os.path.exists(results_file):
            print(f"\nResults saved to: {results_file}")
            print("\nKey output files:")
            print(f"  - Main results: {output_dir}/mmgbsa_results.dat")
            print(f"  - Decomposition: {output_dir}/mmgbsa_decomp.dat")
            print(f"  - Energy components: {output_dir}/mmgbsa_energy.dat")
            print(f"  - Summary: {output_dir}/calculation_summary.txt")
            print(f"  - Complex structure: {output_dir}/complex.pdb")
        
    except Exception as e:
        print(f"\nERROR: Calculation failed: {e}")
        print("Check the log output above for details.")
        return False
    
    return True

def run_mmpbsa_example():
    """Run an example MMPBSA calculation."""
    
    # Example file paths
    protein_file = "3htb/3htb_cleaned_protonated_fixed.pdb"
    ligand_file = "3htb/JZ4.mol2"
    output_dir = "3htb_mmpbsa_results"
    
    print("\n" + "="*50)
    print("RUNNING MMPBSA EXAMPLE")
    print("="*50)
    
    # Create pipeline for MMPBSA
    pipeline = MMGBSAPipeline(
        protein_file=protein_file,
        ligand_file=ligand_file,
        output_dir=output_dir,
        method="mmpbsa",  # Using MMPBSA instead of MMGBSA
        temperature=300,
        equilibration_steps=5000,
        production_steps=100000,
        trajectory_interval=500,
        salt_concentration=0.15,
        start_frame=1,
        end_frame=100,
        frame_interval=1
    )
    
    try:
        pipeline.run()
        print(f"\nMMPBSA calculation completed! Results in: {output_dir}")
    except Exception as e:
        print(f"MMPBSA calculation failed: {e}")

def run_simple_example():
    """Run a very simple example with minimal parameters."""
    
    protein_file = "3htb/3htb_cleaned_protonated_fixed.pdb"
    ligand_file = "3htb/JZ4.mol2"
    output_dir = "simple_results"
    
    print("\n" + "="*50)
    print("RUNNING SIMPLE EXAMPLE")
    print("="*50)
    
    # Create pipeline with minimal parameters
    pipeline = MMGBSAPipeline(
        protein_file=protein_file,
        ligand_file=ligand_file,
        output_dir=output_dir,
        method="mmgbsa"
    )
    
    try:
        pipeline.run()
        print(f"\nSimple calculation completed! Results in: {output_dir}")
    except Exception as e:
        print(f"Simple calculation failed: {e}")

if __name__ == "__main__":
    print("MMGBSA/MMPBSA Pipeline Example")
    print("="*40)
    print("This example demonstrates the simplified pipeline that only requires")
    print("protein and ligand files. The complex will be created automatically")
    print("with proper charging and protonation.")
    print()
    
    # Run MMGBSA example
    success = run_example()
    
    if success:
        # Ask user what they want to run next
        print("\nWhat would you like to run next?")
        print("1. MMPBSA example")
        print("2. Simple example (minimal parameters)")
        print("3. Exit")
        
        choice = input("Enter your choice (1-3): ").strip()
        
        if choice == "1":
            run_mmpbsa_example()
        elif choice == "2":
            run_simple_example()
        elif choice == "3":
            print("Exiting...")
        else:
            print("Invalid choice, exiting...")
    
    print("\nExample completed!") 