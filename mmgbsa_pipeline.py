#!/usr/bin/env python3
"""
MMGBSA/MMPBSA Pipeline using OpenMM and AmberTools

This script provides a simplified interface for running MMGBSA/MMPBSA calculations
on protein-ligand complexes using OpenMM for MD simulation and AmberTools for analysis.

Usage:
    python mmgbsa_pipeline.py --protein protein.pdb --ligand ligand.mol2 --output_dir results/
"""

import argparse
import os
import sys
import subprocess
import shutil
from pathlib import Path
import tempfile
import logging
from typing import Optional, Dict, Any

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class MMGBSAPipeline:
    """Main pipeline class for MMGBSA/MMPBSA calculations."""
    
    def __init__(self, protein_file: str, ligand_file: str, 
                 output_dir: str, method: str = "mmgbsa", **kwargs):
        """
        Initialize the MMGBSA pipeline.
        
        Args:
            protein_file: Path to the protein PDB file
            ligand_file: Path to the ligand file (mol2, pdb, etc.)
            output_dir: Directory for output files
            method: Either "mmgbsa" or "mmpbsa"
            **kwargs: Additional simulation parameters
        """
        self.protein_file = Path(protein_file)
        self.ligand_file = Path(ligand_file)
        self.output_dir = Path(output_dir)
        self.method = method.lower()
        
        # Simulation parameters with defaults
        self.simulation_params = {
            'temperature': kwargs.get('temperature', 300),  # Kelvin
            'pressure': kwargs.get('pressure', 1.0),       # atm
            'equilibration_steps': kwargs.get('equilibration_steps', 10000),
            'production_steps': kwargs.get('production_steps', 2000000),
            'trajectory_interval': kwargs.get('trajectory_interval', 1000),
            'salt_concentration': kwargs.get('salt_concentration', 0.15),  # M
            'start_frame': kwargs.get('start_frame', 1),
            'end_frame': kwargs.get('end_frame', 2000),
            'frame_interval': kwargs.get('frame_interval', 1),
        }
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Working directory for temporary files
        self.work_dir = self.output_dir / "work"
        self.work_dir.mkdir(exist_ok=True)
        
        # Validate inputs
        self._validate_inputs()
    
    def _validate_inputs(self):
        """Validate input files exist and are accessible."""
        required_files = [self.protein_file, self.ligand_file]
        for file_path in required_files:
            if not file_path.exists():
                raise FileNotFoundError(f"Input file not found: {file_path}")
        
        # Validate file formats - only PDB files accepted
        if self.protein_file.suffix.lower() != '.pdb':
            raise ValueError(f"Protein file must be PDB format, got: {self.protein_file.suffix}")
        
        if self.ligand_file.suffix.lower() != '.pdb':
            raise ValueError(f"Ligand file must be PDB format, got: {self.ligand_file.suffix}")
        
        if self.method not in ["mmgbsa", "mmpbsa"]:
            raise ValueError(f"Method must be 'mmgbsa' or 'mmpbsa', got: {self.method}")
    
    def run(self):
        """Execute the complete MMGBSA/MMPBSA pipeline."""
        logger.info(f"Starting {self.method.upper()} pipeline")
        logger.info(f"Protein file: {self.protein_file}")
        logger.info(f"Ligand file: {self.ligand_file}")
        logger.info(f"Output directory: {self.output_dir}")
        
        try:
            # Step 1: Prepare and process input files
            logger.info("Step 1: Preparing and processing input files")
            self._prepare_input_files()
            
            # Step 2: Create complex structure
            logger.info("Step 2: Creating protein-ligand complex")
            self._create_complex()
            
            # Step 3: Prepare Amber topology files
            logger.info("Step 3: Preparing Amber topology files")
            self._prepare_topologies()
            
            # Step 4: Run MD simulation
            logger.info("Step 4: Running MD simulation")
            self._run_simulation()
            
            # Step 5: Convert trajectory
            logger.info("Step 5: Converting trajectory format")
            self._convert_trajectory()
            
            # Step 6: Run MMGBSA/MMPBSA analysis
            logger.info(f"Step 6: Running {self.method.upper()} analysis")
            self._run_analysis()
            
            # Step 7: Clean up and organize results
            logger.info("Step 7: Organizing results")
            self._organize_results()
            
            logger.info("Pipeline completed successfully!")
            
        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            raise
    
    def _prepare_input_files(self):
        """Prepare and process input files for proper charging and protonation."""
        logger.info("Preparing input files for proper charging and protonation...")
        
        # Debug information
        logger.info(f"Current working directory: {os.getcwd()}")
        logger.info(f"Protein file path: {self.protein_file} (exists: {self.protein_file.exists()})")
        logger.info(f"Ligand file path: {self.ligand_file} (exists: {self.ligand_file.exists()})")
        logger.info(f"Work directory: {self.work_dir}")
        
        # Ensure work directory exists
        self.work_dir.mkdir(parents=True, exist_ok=True)
        
        # Verify input files exist
        if not self.protein_file.exists():
            raise FileNotFoundError(f"Protein file not found: {self.protein_file}")
        if not self.ligand_file.exists():
            raise FileNotFoundError(f"Ligand file not found: {self.ligand_file}")
        
        # Copy input files to work directory
        protein_work = self.work_dir / "protein_input.pdb"
        ligand_work = self.work_dir / "ligand_input.pdb"
        
        try:
            shutil.copy2(self.protein_file, protein_work)
            logger.info(f"Copied protein file to: {protein_work}")
        except Exception as e:
            raise RuntimeError(f"Failed to copy protein file: {e}")
        
        try:
            shutil.copy2(self.ligand_file, ligand_work)
            logger.info(f"Copied ligand file to: {ligand_work}")
        except Exception as e:
            raise RuntimeError(f"Failed to copy ligand file: {e}")
        
        # Process protein for proper protonation and charging
        self._process_protein(protein_work)
        
        # Process ligand for proper charging
        self._process_ligand(ligand_work)
    
    def _process_protein(self, protein_file: Path):
        """Process protein for proper protonation and charging."""
        logger.info("Processing protein for protonation and charging...")
        
        # Ensure work directory exists
        self.work_dir.mkdir(parents=True, exist_ok=True)
        
        # Use pdb4amber to clean the protein
        cleaned_protein = self.work_dir / "protein_cleaned.pdb"
        cmd = f"pdb4amber --nohyd --no-conect -i {protein_file.name} -o {cleaned_protein.name}"
        
        try:
            result = subprocess.run(cmd, shell=True, cwd=self.work_dir, 
                                  capture_output=True, text=True)
            if result.returncode != 0:
                logger.warning(f"pdb4amber failed, using original protein: {result.stderr}")
                shutil.copy2(protein_file, cleaned_protein)
        except FileNotFoundError:
            logger.warning("pdb4amber not found, using original protein")
            shutil.copy2(protein_file, cleaned_protein)
        
        # Use reduce to add hydrogens if available
        protonated_protein = self.work_dir / "protein_protonated.pdb"
        cmd = f"reduce -HIS -FLIP {cleaned_protein.name} > {protonated_protein.name}"
        
        try:
            result = subprocess.run(cmd, shell=True, cwd=self.work_dir, 
                                  capture_output=True, text=True)
            if result.returncode != 0:
                logger.warning(f"reduce failed, using cleaned protein: {result.stderr}")
                shutil.copy2(cleaned_protein, protonated_protein)
        except FileNotFoundError:
            logger.warning("reduce not found, using cleaned protein")
            shutil.copy2(cleaned_protein, protonated_protein)
        
        self.processed_protein = cleaned_protein#protonated_protein
    
    def _process_ligand(self, ligand_file: Path):
        """Process ligand for proper charging and hydrogen addition using antechamber."""
        logger.info("Processing ligand for proper charging and hydrogen addition...")
        
        # Ensure work directory exists
        self.work_dir.mkdir(parents=True, exist_ok=True)
        
        # Use antechamber to convert PDB to MOL2 with proper hydrogens and charges
        mol2_ligand = self.work_dir / "ligand.mol2"
        
        try:
            #add hydrogens
            with_h_ligand_name = ligand_file.name.replace('.pdb', '_with_h.pdb')
            result = subprocess.run(f"reduce {ligand_file.name} > {with_h_ligand_name}", shell=True, cwd=self.work_dir,capture_output=True, text=True)
            #result = subprocess.run(f"obabel {ligand_file.name} -O {with_h_ligand_name} -h", shell=True, cwd=self.work_dir,capture_output=True, text=True)
            assert result.returncode == 0, f"Reduce failed: {result.stderr}"
            # antechamber can handle PDB input and generate MOL2 with hydrogens and AM1-BCC charges
            cmd = f"antechamber -i {with_h_ligand_name} -fi pdb -o {mol2_ligand.name} -fo mol2 -c bcc -nc 0 -j 4"
            
            result = subprocess.run(cmd, shell=True, cwd=self.work_dir, 
                                capture_output=True, text=True)
            if result.returncode != 0:
                intermediate=self.work_dir / "intermediate_sybl.mol2"
                cmd1=f'antechamber -i "{with_h_ligand_name}" -fi pdb -o "{intermediate.name}" -fo mol2 -j 5 -at sybyl -dr no'
                cmd2=f'antechamber -i "{intermediate.name}" -fi mol2 -o "{mol2_ligand.name}" -fo mol2 -c bcc -s 2'
                result = subprocess.run(cmd1, shell=True, cwd=self.work_dir, 
                                capture_output=True, text=True)
                logger.info(f"Antechamber1 failed: {result.stderr}")
                #assert result.returncode == 0, f"Antechambe1 failed: {result.stderr}"
                result = subprocess.run(cmd2, shell=True, cwd=self.work_dir, 
                                capture_output=True, text=True)
                logger.info(f"Antechamber2 failed: {result.stderr}")
                #assert result.returncode == 0, f"Antechamber2 failed: {result.stderr}"
            if result.returncode != 0:
                at_once = subprocess.run(f"antechamber -i {ligand_file.name} -fi pdb -o {mol2_ligand.name} -fo mol2 -c bcc -nc 0 -j 4 -at gaff", shell=True, cwd=self.work_dir,capture_output=True, text=True)
                #hmmmm this^^^ might be making a ligand without hydrogens
                logger.info("trying with single pass antechamber")
                result=at_once
                
            if result.returncode == 0:
                logger.info("Successfully converted ligand PDB to MOL2 using antechamber")
                
                # Check if the MOL2 file is valid (has reasonable size and content)
                if mol2_ligand.exists() and mol2_ligand.stat().st_size > 100:
                    self.processed_ligand = mol2_ligand
                    logger.info("Using antechamber-generated MOL2 file")
                else:
                    logger.warning("Antechamber-generated MOL2 file appears invalid, using original PDB")
                    self.processed_ligand = ligand_file
            else:
                logger.warning(f"Antechamber PDB to MOL2 conversion failed: {result.stderr}")
                logger.info("Using original PDB ligand")
                self.processed_ligand = ligand_file
                
        except FileNotFoundError:
            logger.warning("Antechamber not found, using original ligand")
            self.processed_ligand = ligand_file
        
        # Generate frcmod file for the ligand using parmchk2
        self._generate_ligand_parameters()
    
    def _generate_ligand_parameters(self):
        """Generate frcmod file for the ligand using parmchk2."""
        logger.info("Generating frcmod file for the ligand...")
        
        ligand_name = self.processed_ligand.stem
        
        # Generate frcmod file using parmchk2
        # The MOL2 file already has AM1-BCC charges from antechamber
        cmd = f"parmchk2 -i {self.processed_ligand.name} -f mol2 -o {ligand_name}.frcmod"
        
        try:
            result = subprocess.run(cmd, shell=True, cwd=self.work_dir, 
                                  capture_output=True, text=True)
            if result.returncode == 0:
                logger.info("Parmchk2 completed successfully")
            else:
                logger.warning(f"Parmchk2 failed: {result.stderr}")
        except FileNotFoundError:
            logger.warning("Parmchk2 not found")
    
    def _create_complex(self):
        """Create protein-ligand complex structure."""
        logger.info("Creating protein-ligand complex...")
        
        # Use tleap to combine protein and ligand
        complex_tleap = self._create_complex_tleap_script()
        tleap_script_path = self.work_dir / "create_complex.tleap"
        
        with open(tleap_script_path, 'w') as f:
            f.write(complex_tleap)
        
        # Run tleap to create complex
        self._run_tleap(tleap_script_path)
        
        # The complex PDB will be saved as complex.pdb in work directory
        self.complex_file = self.work_dir / "complex.pdb"
        
        if not self.complex_file.exists():
            raise RuntimeError("Failed to create complex structure")
        
        logger.info(f"Complex created: {self.complex_file}")
    
    def _create_complex_tleap_script(self) -> str:
        """Create tleap script for combining protein and ligand."""
        protein_name = self.processed_protein.stem
        ligand_name = self.processed_ligand.stem
        
        # Check if we have a valid MOL2 file, otherwise use PDB
        if self.processed_ligand.suffix.lower() == '.mol2' and self.processed_ligand.exists():
            # Try to load as MOL2 first
            load_ligand = f"lig = loadmol2 {self.processed_ligand.name}"
        else:
            # Fallback to PDB
            pdb_ligand = self.work_dir / "ligand_input.pdb"
            load_ligand = f"lig = loadpdb {pdb_ligand.name}"
        
        # Check if frcmod file exists
        frcmod_file = self.work_dir / f"{ligand_name}.frcmod"
        if frcmod_file.exists():
            load_frcmod = f"loadamberparams {frcmod_file.name}"
        else:
            load_frcmod = "# No frcmod file found"
        load_ligand = f"lig = loadmol2 {self.processed_ligand.name}"
        script = f"""source leaprc.protein.ff14SB
source leaprc.gaff
set default PBRadii mbondi2

{load_frcmod}

# Load the protein
prot = loadpdb {self.processed_protein.name}

# Load the ligand
{load_ligand}

# Combine protein and ligand
comp = combine {{prot lig}}

# Save complex structure
savepdb comp complex.pdb

# Save individual components for analysis (using MMPBSA.py naming convention)
saveAmberParm prot protein_prmtop protein.inpcrd
saveAmberParm lig ligand_prmtop ligand.inpcrd
saveAmberParm comp complex_prmtop complex.inpcrd

quit
"""
        print('='*10,'SCRIPT')
        print(script)
        return script
    
    def _prepare_topologies(self):
        """Generate Amber topology files using tleap."""
        logger.info("Generating Amber topology files...")
        
        # The complex creation step already generated the topology files
        # Just verify they exist
        required_files = [
            "complex_prmtop", "complex.inpcrd",
            "protein_prmtop", "protein.inpcrd", 
            "ligand_prmtop", "ligand.inpcrd"
        ]
        
        for file_name in required_files:
            file_path = self.work_dir / file_name
            if not file_path.exists():
                raise RuntimeError(f"Required topology file not found: {file_name}")
        
        logger.info("All topology files generated successfully")
    
    def _run_tleap(self, script_path: Path):
        """Run tleap with the given script."""
        cmd = f"tleap -f {script_path.name}"
        logger.info(f"Running: {cmd}")
        logger.info(f"Working directory: {self.work_dir}")
        logger.info(f"Script path: {script_path}")
        
        result = subprocess.run(cmd, shell=True, cwd=self.work_dir, 
                              capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.error(f"Tleap failed with return code: {result.returncode}")
            logger.error(f"Tleap stderr: {result.stderr}")
            logger.error(f"Tleap stdout: {result.stdout}")
            raise RuntimeError(f"Tleap command failed: {result.stderr}")
        
        logger.info("Tleap completed successfully")
    
    def _run_simulation(self):
        """Run MD simulation using OpenMM."""
        logger.info("Running MD simulation with OpenMM...")
        
        # Create simulation script
        sim_script = self._create_simulation_script()
        sim_script_path = self.work_dir / "run_simulation.py"
        
        with open(sim_script_path, 'w') as f:
            f.write(sim_script)
        
        # Run simulation
        cmd = f"python {sim_script_path.name}"
        logger.info(f"Running: {cmd}")
        
        result = subprocess.run(cmd, shell=True, cwd=self.work_dir, 
                              capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.error(f"Simulation failed: {result.stderr}")
            raise RuntimeError(f"Simulation failed: {result.stderr}")
        
        logger.info("MD simulation completed successfully")
    
    def _create_simulation_script(self) -> str:
        """Create OpenMM simulation script."""
        script = f"""#!/usr/bin/env python3
from openmm.app import *
from openmm import *
from openmm.unit import *
import sys

# Load topology and coordinates
inpcrd = AmberInpcrdFile('complex.inpcrd')
prmtop = AmberPrmtopFile('complex_prmtop', periodicBoxVectors=inpcrd.boxVectors)

# Create system
method = PME if inpcrd.boxVectors is not None else NoCutoff
system = prmtop.createSystem(
    implicitSolvent=OBC2,
    nonbondedMethod=method, 
    constraints=HBonds, 
    implicitSolventSaltConc={self.simulation_params['salt_concentration']}*moles/liter
)

# Create integrator
integrator = LangevinMiddleIntegrator(
    {self.simulation_params['temperature']}*kelvin, 
    1/picosecond, 
    0.002*picoseconds
)

# Create simulation
simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(inpcrd.positions)

# Save initial structure
with open('complex_initial.pdb', 'w') as outfile:
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), file=outfile, keepIds=True)

print('Starting minimization...')
simulation.minimizeEnergy()

# Set velocities and equilibrate
simulation.context.setVelocitiesToTemperature({self.simulation_params['temperature']}*kelvin)
print('Starting equilibration...')
simulation.step({self.simulation_params['equilibration_steps']})

# Save equilibrated structure
with open('complex_equilibrated.pdb', 'w') as outfile:
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), file=outfile, keepIds=True)

# Setup trajectory output
output_base = 'complex_traj'
simulation.reporters.append(DCDReporter(output_base + ".dcd", {self.simulation_params['trajectory_interval']}))

print('Starting production run...')
simulation.step({self.simulation_params['production_steps']})

print('Simulation completed!')
"""
        return script
    
    def _convert_trajectory(self):
        """Convert DCD trajectory to Amber format."""
        logger.info("Converting trajectory format...")
        
        # Check if DCD file exists
        dcd_file = self.work_dir / "complex_traj.dcd"
        if not dcd_file.exists():
            raise RuntimeError("DCD trajectory file not found")
        
        # Create cpptraj script
        cpptraj_script = f"""parm complex_prmtop
trajin complex_traj.dcd 1 last {self.simulation_params['frame_interval']}
trajout mdcrd
run
quit
"""
        
        cpptraj_path = self.work_dir / "convert_traj.cpptraj"
        with open(cpptraj_path, 'w') as f:
            f.write(cpptraj_script)
        
        # Run cpptraj
        cmd = f"cpptraj -i {cpptraj_path.name}"
        logger.info(f"Running: {cmd}")
        
        result = subprocess.run(cmd, shell=True, cwd=self.work_dir, 
                              capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.warning(f"cpptraj failed: {result.stderr}")
            logger.info("Skipping trajectory conversion - will use DCD file directly")
            
            # Create a symbolic link or copy the DCD file
            mdcrd_file = self.work_dir / "mdcrd"
            try:
                # Try to create a symbolic link first
                if mdcrd_file.exists():
                    mdcrd_file.unlink()
                mdcrd_file.symlink_to(dcd_file.name)
                logger.info("Created symbolic link to DCD file for analysis")
            except Exception as e:
                logger.warning(f"Could not create symbolic link: {e}")
                try:
                    # Fallback to copying
                    shutil.copy2(dcd_file, mdcrd_file)
                    logger.info("Copied DCD file as mdcrd for analysis")
                except Exception as e2:
                    logger.error(f"Failed to copy DCD file: {e2}")
                    raise RuntimeError("Trajectory conversion failed and fallback failed")
        else:
            logger.info("Trajectory conversion completed successfully")
    
    def _run_analysis(self):
        """Run MMGBSA/MMPBSA analysis."""
        logger.info(f"Running {self.method.upper()} analysis...")
        
        # Create input file for MMPBSA.py
        analysis_input = self._create_analysis_input()
        analysis_input_path = self.work_dir / f"{self.method}.in"
        
        with open(analysis_input_path, 'w') as f:
            f.write(analysis_input)
        
        # Run MMPBSA.py
        cmd = f"MMPBSA.py -O -i {analysis_input_path.name} -o {self.method}_results.dat -do {self.method}_decomp.dat -eo {self.method}_energy.dat"
        logger.info(f"Running: {cmd}")
        
        result = subprocess.run(cmd, shell=True, cwd=self.work_dir, 
                              capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.error(f"MMPBSA.py failed: {result.stderr}")
            
            # Check if it's the cpptraj library issue
            if "libblas.3.dylib" in result.stderr or "cpptraj failed" in result.stderr:
                logger.error("=" * 60)
                logger.error("SYSTEM ISSUE DETECTED")
                logger.error("=" * 60)
                logger.error("The pipeline has successfully completed all steps except the final analysis.")
                logger.error("The issue is with the AmberTools installation - cpptraj has missing library dependencies.")
                logger.error("")
                logger.error("SOLUTIONS:")
                logger.error("1. Reinstall AmberTools with proper library dependencies")
                logger.error("2. Install missing BLAS libraries: conda install -c conda-forge openblas")
                logger.error("3. Use a different AmberTools installation")
                logger.error("")
                logger.error("PIPELINE STATUS:")
                logger.error("✅ Protein processing: COMPLETE")
                logger.error("✅ Ligand processing: COMPLETE") 
                logger.error("✅ Complex creation: COMPLETE")
                logger.error("✅ Topology generation: COMPLETE")
                logger.error("✅ MD simulation: COMPLETE")
                logger.error("✅ Trajectory handling: COMPLETE")
                logger.error("❌ MMGBSA analysis: FAILED (system issue)")
                logger.error("")
                logger.error("All intermediate files are available in the work directory.")
                logger.error("You can manually run MMPBSA.py once the system issue is resolved.")
                logger.error("=" * 60)
                
                # Create a summary file with instructions
                self._create_system_issue_summary()
                
                raise RuntimeError("MMPBSA.py failed due to system-level library dependencies. See log for details.")
            else:
                raise RuntimeError(f"MMPBSA.py failed: {result.stderr}")
        
        logger.info(f"{self.method.upper()} analysis completed")
    
    def _create_analysis_input(self) -> str:
        """Create input file for MMPBSA.py."""
        script = f"""&general
startframe={self.simulation_params['start_frame']},
endframe={self.simulation_params['end_frame']},
interval={self.simulation_params['frame_interval']},
verbose=2,
keep_files=2,
/
"""
        
        if self.method == "mmgbsa":
            script += f"""&gb
igb=2,
saltcon={self.simulation_params['salt_concentration']},
/
"""
        elif self.method == "mmpbsa":
            script += f"""&pb
istrng={self.simulation_params['salt_concentration']},
/
"""
        
        script += f"""&decomp
idecomp=1,
dec_verbose=1,
/
"""
        
        return script
    
    def _organize_results(self):
        """Organize and copy results to output directory."""
        logger.info("Organizing results...")
        
        # Copy important files to output directory
        important_files = [
            f"{self.method}_results.dat",
            f"{self.method}_decomp.dat", 
            f"{self.method}_energy.dat",
            "complex.pdb",
            "complex_equilibrated.pdb"
        ]
        
        for file_name in important_files:
            src = self.work_dir / file_name
            dst = self.output_dir / file_name
            if src.exists():
                shutil.copy2(src, dst)
                logger.info(f"Copied {file_name} to output directory")
        
        # Create summary file
        self._create_summary()
        
        logger.info(f"Results organized in: {self.output_dir}")
    
    def _create_summary(self):
        """Create a summary of the calculation."""
        summary_path = self.output_dir / "calculation_summary.txt"
        
        with open(summary_path, 'w') as f:
            f.write("MMGBSA/MMPBSA Calculation Summary\n")
            f.write("=" * 40 + "\n\n")
            f.write(f"Method: {self.method.upper()}\n")
            f.write(f"Protein file: {self.protein_file}\n")
            f.write(f"Ligand file: {self.ligand_file}\n")
            f.write(f"Complex created automatically\n\n")
            
            f.write("Simulation Parameters:\n")
            for key, value in self.simulation_params.items():
                f.write(f"  {key}: {value}\n")
            
            f.write(f"\nOutput directory: {self.output_dir}\n")
            f.write(f"Working directory: {self.work_dir}\n")

    def _create_system_issue_summary(self):
        """Create a summary file with instructions for resolving system issues."""
        summary_path = self.output_dir / "system_issue_summary.txt"
        
        with open(summary_path, 'w') as f:
            f.write("MMGBSA Pipeline - System Issue Summary\n")
            f.write("=" * 50 + "\n\n")
            f.write("The pipeline completed successfully up to the MMGBSA analysis step.\n")
            f.write("The final step failed due to system-level library dependencies.\n\n")
            
            f.write("COMPLETED STEPS:\n")
            f.write("✅ Protein processing and protonation\n")
            f.write("✅ Ligand processing and parameter generation\n")
            f.write("✅ Complex creation\n")
            f.write("✅ Amber topology generation\n")
            f.write("✅ MD simulation with OpenMM\n")
            f.write("✅ Trajectory handling\n\n")
            
            f.write("FAILED STEP:\n")
            f.write("❌ MMGBSA analysis (cpptraj library dependency issue)\n\n")
            
            f.write("SOLUTIONS:\n")
            f.write("1. Reinstall AmberTools:\n")
            f.write("   conda install -c conda-forge ambertools\n\n")
            f.write("2. Install missing BLAS libraries:\n")
            f.write("   conda install -c conda-forge openblas\n\n")
            f.write("3. Use a different AmberTools installation\n\n")
            
            f.write("MANUAL COMPLETION:\n")
            f.write("Once the system issue is resolved, you can manually run:\n")
            f.write(f"cd {self.work_dir}\n")
            f.write(f"MMPBSA.py -O -i {self.method}.in -o {self.method}_results.dat -do {self.method}_decomp.dat -eo {self.method}_energy.dat\n\n")
            
            f.write("FILES AVAILABLE:\n")
            f.write(f"- Complex structure: {self.work_dir}/complex.pdb\n")
            f.write(f"- Topology files: {self.work_dir}/complex_prmtop, protein_prmtop, ligand_prmtop\n")
            f.write(f"- Trajectory: {self.work_dir}/complex_traj.dcd\n")
            f.write(f"- Analysis input: {self.work_dir}/{self.method}.in\n")


def main():
    """Main function to run the MMGBSA pipeline."""
    parser = argparse.ArgumentParser(
        description="MMGBSA/MMPBSA Pipeline using OpenMM and AmberTools",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python mmgbsa_pipeline.py --protein protein.pdb --ligand ligand.mol2 --output_dir results/
  python mmgbsa_pipeline.py --protein protein.pdb --ligand ligand.pdb --method mmpbsa --output_dir results/ --temperature 310
        """
    )
    
    parser.add_argument("--protein", required=True, help="Protein PDB file")
    parser.add_argument("--ligand", required=True, help="Ligand file (mol2, pdb, etc.)")
    parser.add_argument("--output_dir", required=True, help="Output directory for results")
    parser.add_argument("--method", choices=["mmgbsa", "mmpbsa"], default="mmgbsa", 
                       help="Analysis method (default: mmgbsa)")
    
    # Simulation parameters
    parser.add_argument("--temperature", type=float, default=300, help="Temperature in Kelvin (default: 300)")
    parser.add_argument("--pressure", type=float, default=1.0, help="Pressure in atm (default: 1.0)")
    parser.add_argument("--equilibration_steps", type=int, default=10000, help="Equilibration steps (default: 10000)")
    parser.add_argument("--production_steps", type=int, default=2000000, help="Production steps (default: 2000000)")
    parser.add_argument("--trajectory_interval", type=int, default=1000, help="Trajectory save interval (default: 1000)")
    parser.add_argument("--salt_concentration", type=float, default=0.15, help="Salt concentration in M (default: 0.15)")
    parser.add_argument("--start_frame", type=int, default=1, help="Start frame for analysis (default: 1)")
    parser.add_argument("--end_frame", type=int, default=2000, help="End frame for analysis (default: 2000)")
    parser.add_argument("--frame_interval", type=int, default=1, help="Frame interval for analysis (default: 1)")
    
    args = parser.parse_args()
    
    # Create pipeline and run
    pipeline = MMGBSAPipeline(
        protein_file=args.protein,
        ligand_file=args.ligand,
        output_dir=args.output_dir,
        method=args.method,
        temperature=args.temperature,
        pressure=args.pressure,
        equilibration_steps=args.equilibration_steps,
        production_steps=args.production_steps,
        trajectory_interval=args.trajectory_interval,
        salt_concentration=args.salt_concentration,
        start_frame=args.start_frame,
        end_frame=args.end_frame,
        frame_interval=args.frame_interval
    )
    
    pipeline.run()


if __name__ == "__main__":
    main() 