#!/usr/bin/env python3
"""
Script to run MMGBSA/MMPBSA pipeline from a YAML configuration file.

Usage:
    python run_from_config.py config.yaml
"""

import argparse
import yaml
import sys
from pathlib import Path

# Add the current directory to Python path
sys.path.insert(0, str(Path(__file__).parent))

from mmgbsa_pipeline import MMGBSAPipeline

def load_config(config_file):
    """Load configuration from YAML file."""
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    return config

def validate_config(config):
    """Validate configuration parameters."""
    required_sections = ['input', 'output', 'simulation', 'analysis']
    
    for section in required_sections:
        if section not in config:
            raise ValueError(f"Missing required section: {section}")
    
    # Check input files
    input_files = config['input']
    for file_type, file_path in input_files.items():
        if not Path(file_path).exists():
            print(f"Warning: Input file {file_path} not found")
    
    # Validate method
    method = config['output']['method']
    if method not in ['mmgbsa', 'mmpbsa']:
        raise ValueError(f"Invalid method: {method}. Must be 'mmgbsa' or 'mmpbsa'")
    
    return True

def run_from_config(config_file):
    """Run the pipeline from a configuration file."""
    print(f"Loading configuration from: {config_file}")
    
    # Load and validate configuration
    config = load_config(config_file)
    validate_config(config)
    
    # Extract parameters
    input_params = config['input']
    output_params = config['output']
    simulation_params = config['simulation']
    analysis_params = config['analysis']
    
    # Create pipeline
    pipeline = MMGBSAPipeline(
        protein_file=input_params['protein_file'],
        ligand_file=input_params['ligand_file'],
        output_dir=output_params['output_dir'],
        method=output_params['method'],
        temperature=simulation_params['temperature'],
        pressure=simulation_params['pressure'],
        equilibration_steps=simulation_params['equilibration_steps'],
        production_steps=simulation_params['production_steps'],
        trajectory_interval=simulation_params['trajectory_interval'],
        salt_concentration=simulation_params['salt_concentration'],
        start_frame=analysis_params['start_frame'],
        end_frame=analysis_params['end_frame'],
        frame_interval=analysis_params['frame_interval']
    )
    
    # Run pipeline
    pipeline.run()
    
    print(f"\nPipeline completed! Results saved to: {output_params['output_dir']}")

def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Run MMGBSA/MMPBSA pipeline from YAML configuration",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example configuration file (config.yaml):
  input:
    protein_file: "protein.pdb"
    ligand_file: "ligand.mol2"
  output:
    output_dir: "results/"
    method: "mmgbsa"
  simulation:
    temperature: 300
    production_steps: 2000000
  analysis:
    start_frame: 1
    end_frame: 2000

Note: The complex will be created automatically from protein and ligand files.
        """
    )
    
    parser.add_argument("config_file", help="YAML configuration file")
    parser.add_argument("--check", action="store_true", 
                       help="Only validate configuration, don't run")
    
    args = parser.parse_args()
    
    if not Path(args.config_file).exists():
        print(f"Error: Configuration file {args.config_file} not found")
        sys.exit(1)
    
    try:
        if args.check:
            # Just validate configuration
            config = load_config(args.config_file)
            validate_config(config)
            print("Configuration is valid!")
        else:
            # Run the pipeline
            run_from_config(args.config_file)
            
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 