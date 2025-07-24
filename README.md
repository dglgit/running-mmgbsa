# MMGBSA/MMPBSA Pipeline

A simplified pipeline for performing MMGBSA (Molecular Mechanics Generalized Born Surface Area) and MMPBSA (Molecular Mechanics Poisson-Boltzmann Surface Area) calculations using OpenMM and AmberTools.

## 🎯 **Key Features**

- **Ultra-Simple Interface**: Only requires PDB files for protein and ligand
- **Automatic Processing**: Handles complex creation, charging, and protonation
- **Consistent Charge Assignment**: Uses antechamber for AM1-BCC charge assignment
- **Robust Error Handling**: Graceful fallbacks and comprehensive logging
- **Multiple Usage Modes**: Command line, configuration file, and programmatic API

## 📋 **Requirements**

### System Dependencies
- Python 3.7+
- OpenMM
- AmberTools (tleap, cpptraj, MMPBSA.py, antechamber, parmchk2)
- Open Babel (obabel) - optional, for fallback processing

### Python Dependencies
```
numpy
openmm
pyyaml
```

## 🚀 **Quick Start**

### Basic Usage (PDB Files Only)

```bash
python mmgbsa_pipeline.py \
    --protein protein.pdb \
    --ligand ligand.pdb \
    --output_dir results/
```

### Advanced Usage

```bash
python mmgbsa_pipeline.py \
    --protein protein.pdb \
    --ligand ligand.pdb \
    --output_dir results/ \
    --method mmgbsa \
    --production_steps 10000 \
    --end_frame 100 \
    --frame_interval 10
```

## 📁 **Input Files**

### Required Files
- **Protein PDB**: Clean protein structure in PDB format
- **Ligand PDB**: Ligand structure in PDB format (will be automatically processed)

### File Requirements
- Both files must be in **PDB format** (`.pdb` extension)
- Protein should be reasonably clean (missing atoms will be added)
- Ligand will be automatically processed for hydrogen addition

## 🔧 **Pipeline Steps**

1. **File Validation**: Ensures PDB files exist and are valid
2. **Protein Processing**: Cleaning and protonation using pdb4amber and reduce
3. **Ligand Processing**: Hydrogen addition and AM1-BCC charge assignment using antechamber
4. **Complex Creation**: Automatic protein-ligand combination using tleap
5. **Topology Generation**: Amber topology files for protein, ligand, and complex
6. **MD Simulation**: Molecular dynamics using OpenMM
7. **Trajectory Processing**: Format conversion for analysis
8. **MMGBSA/MMPBSA Analysis**: Binding free energy calculation

## 📊 **Output Files**

### Main Results
- `{method}_results.dat`: Binding free energy results
- `{method}_decomp.dat`: Per-residue decomposition
- `{method}_energy.dat`: Detailed energy components

### Intermediate Files (in `work/` directory)
- `complex.pdb`: Protein-ligand complex structure
- `complex_prmtop`, `complex.inpcrd`: Complex topology files
- `protein_prmtop`, `protein.inpcrd`: Protein topology files
- `ligand_prmtop`, `ligand.inpcrd`: Ligand topology files
- `complex_traj.dcd`: MD trajectory
- `{method}.in`: Analysis input file

## ⚙️ **Configuration**

### Command Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `--protein` | Protein PDB file | Required |
| `--ligand` | Ligand PDB file | Required |
| `--output_dir` | Output directory | Required |
| `--method` | Analysis method (mmgbsa/mmpbsa) | mmgbsa |
| `--production_steps` | MD production steps | 10000 |
| `--end_frame` | Number of frames for analysis | 100 |
| `--frame_interval` | Frame interval for analysis | 10 |

### Configuration File

Create a YAML configuration file:

```yaml
# config.yaml
protein_file: "protein.pdb"
ligand_file: "ligand.pdb"
output_dir: "results"
method: "mmgbsa"
simulation_params:
  production_steps: 10000
  end_frame: 100
  frame_interval: 10
```

Run with configuration:
```bash
python run_from_config.py config.yaml
```

## 🔍 **Troubleshooting**

### Common Issues

1. **Antechamber Library Errors**: Missing BLAS/LAPACK libraries
   ```
   dyld: Library not loaded: @rpath/libblas.3.dylib
   ```
   - **Solution**: Install missing system libraries
     ```bash
     # macOS with Homebrew
     brew install openblas lapack
     
     # Or reinstall AmberTools with proper dependencies
     conda install -c conda-forge ambertools
     ```
   - **Pipeline Behavior**: Falls back to original PDB ligand if antechamber fails

2. **cpptraj Library Issues**: Missing system libraries
   - **Solution**: Pipeline automatically falls back to DCD files
   - **Analysis**: Continues with available trajectory format

3. **Atom Type Issues**: Incompatible ligand atom types
   - **Solution**: Pipeline detects and falls back to PDB format
   - **Recommendation**: Use properly formatted ligand files

### System Requirements

- **macOS**: May require additional library installations (BLAS/LAPACK)
- **Linux**: Generally works out-of-the-box
- **Windows**: Use WSL or similar environment

### Charge Assignment Method

The pipeline uses **AM1-BCC** (Austin Model 1 with Bond Charge Correction) for ligand charge assignment:
- **Consistent**: Same method used for both MOL2 and frcmod files
- **Standard**: Widely accepted for biomolecular simulations
- **Reliable**: Implemented in antechamber with proper parameterization

## 📚 **Examples**

### Example 1: Basic MMGBSA
```bash
python mmgbsa_pipeline.py \
    --protein 3htb/3htb_cleaned_protonated_fixed.pdb \
    --ligand 3htb/JZ4.pdb \
    --output_dir 3htb_results
```

### Example 2: MMPBSA with Custom Parameters
```bash
python mmgbsa_pipeline.py \
    --protein protein.pdb \
    --ligand ligand.pdb \
    --output_dir results/ \
    --method mmpbsa \
    --production_steps 50000 \
    --end_frame 500
```

### Example 3: Programmatic Usage
```python
from mmgbsa_pipeline import MMGBSAPipeline

pipeline = MMGBSAPipeline(
    protein_file="protein.pdb",
    ligand_file="ligand.pdb", 
    output_dir="results/",
    method="mmgbsa",
    production_steps=10000,
    end_frame=100
)

pipeline.run()
```

## 🔬 **Technical Details**

### Force Fields
- **Protein**: ff14SB
- **Ligand**: GAFF (General Amber Force Field)
- **Solvent**: Implicit solvent (GB/PB)

### Analysis Methods
- **MMGBSA**: Generalized Born surface area
- **MMPBSA**: Poisson-Boltzmann surface area

### Trajectory Processing
- **Input**: DCD format from OpenMM
- **Output**: mdcrd format for MMPBSA.py
- **Fallback**: Direct DCD usage if conversion fails

## 📄 **License**

This project is open source and available under the MIT License.

## 🤝 **Contributing**

Contributions are welcome! Please feel free to submit issues and pull requests.

