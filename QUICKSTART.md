# Quick Start Guide

This guide will help you get the MMGBSA pipeline running quickly with the new simplified interface.

## Prerequisites

1. **Install AmberTools**
   ```bash
   # Download from https://ambermd.org/GetAmber.php
   # Ensure tleap, cpptraj, and MMPBSA.py are in your PATH
   ```

2. **Install Python dependencies**
   ```bash
   pip install -r requirements.txt
   ```

3. **Optional: Install enhanced processing tools**
   ```bash
   # For better protein processing
   conda install -c conda-forge pdb4amber reduce
   
   # For ligand processing
   conda install -c conda-forge openbabel
   ```

## Quick Test

Run the test script to verify everything is set up correctly:

```bash
python test_pipeline.py
```

## Basic Usage

### Method 1: Command Line (Simplest)

```bash
python mmgbsa_pipeline.py \
    --protein protein.pdb \
    --ligand ligand.mol2 \
    --output_dir results/
```

**That's it!** The pipeline will automatically:
- Clean and protonate your protein
- Process your ligand for proper charging
- Create the protein-ligand complex
- Run the MD simulation
- Perform MMGBSA analysis

### Method 2: Configuration File (Recommended)

1. Copy the template:
   ```bash
   cp config_template.yaml my_config.yaml
   ```

2. Edit `my_config.yaml` with your file paths and parameters

3. Run with configuration:
   ```bash
   python run_from_config.py my_config.yaml
   ```

### Method 3: Python Script

```python
from mmgbsa_pipeline import MMGBSAPipeline

pipeline = MMGBSAPipeline(
    protein_file="protein.pdb",
    ligand_file="ligand.mol2",
    output_dir="results/",
    method="mmgbsa"
)

pipeline.run()
```

## Example with Your Files

If you have the example files in this repository:

```bash
# Test run with reduced parameters
python mmgbsa_pipeline.py \
    --protein 3htb-example/3htb_protein.pdb \
    --ligand 3htb-example/3htb_ligand.pdb \
    --output_dir 3htb_results \
    --production_steps 1000 \
    --end_frame 100
```

## Expected Output

After running, you'll find:

```
3htb_results/
├── mmgbsa_results.dat      # Main binding energy results
├── mmgbsa_decomp.dat       # Per-residue contributions  
├── mmgbsa_energy.dat       # Detailed energy components
├── calculation_summary.txt # Summary of parameters used
├── complex.pdb            # Created protein-ligand complex
└── work/                  # All intermediate files
```

## What the Pipeline Does Automatically

1. **Protein Processing**:
   - Cleans the protein structure
   - Adds missing hydrogens
   - Assigns proper protonation states

2. **Ligand Processing**:
   - Converts formats if needed
   - Adds hydrogens
   - Generates 3D coordinates if needed

3. **Complex Creation**:
   - Combines protein and ligand
   - Generates proper topology files
   - Creates coordinate files

4. **Simulation & Analysis**:
   - Runs MD simulation
   - Performs MMGBSA/MMPBSA analysis
   - Organizes results

## Common Issues

1. **"tleap not found"**: Install AmberTools and add to PATH
2. **"OpenMM import failed"**: Run `pip install openmm`
3. **"Input file not found"**: Check file paths are correct
4. **"Permission denied"**: Ensure write permissions for output directory
5. **"Protein processing failed"**: Install optional tools (pdb4amber, reduce) or use cleaner input files

## Next Steps

- Read the full [README.md](README.md) for detailed documentation
- Check [example_run.py](example_run.py) for more examples
- Use `python test_pipeline.py` to diagnose issues

## Getting Help

If you encounter issues:

1. Run `python test_pipeline.py` to check your setup
2. Check the log output for error messages
3. Verify all input files exist and are properly formatted
4. Ensure AmberTools and OpenMM are properly installed
5. Install optional tools for better file processing 