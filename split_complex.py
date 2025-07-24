#!/usr/bin/env python3
"""
Split a protein-ligand complex PDB file into separate protein and ligand files.

This script takes a PDB file containing a protein-ligand complex and splits it into
separate files for the protein and ligand, preserving their positions and handling
CONECT records correctly.

Usage:
    python split_complex.py complex.pdb LIGAND_NAME output_prefix
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple


class PDBComplexSplitter:
    """Split a protein-ligand complex PDB file into separate components."""
    
    def __init__(self, pdb_file: str, ligand_name: str, output_prefix: str):
        self.pdb_file = Path(pdb_file)
        self.ligand_name = ligand_name.upper()
        self.output_prefix = output_prefix
        
        # Store atoms and connections
        self.protein_atoms = []
        self.ligand_atoms = []
        self.protein_connections = []
        self.ligand_connections = []
        self.other_records = []
        
        # Track atom numbers for CONECT records
        self.protein_atom_numbers = set()
        self.ligand_atom_numbers = set()
        
    def parse_pdb(self):
        """Parse the PDB file and separate protein and ligand atoms."""
        if not self.pdb_file.exists():
            raise FileNotFoundError(f"PDB file not found: {self.pdb_file}")
        
        print(f"Parsing PDB file: {self.pdb_file}")
        print(f"Looking for ligand: {self.ligand_name}")
        
        with open(self.pdb_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.rstrip()
                
                if not line:
                    continue
                
                record_type = line[:6].strip()
                
                if record_type == "ATOM" or record_type == "HETATM":
                    self._process_atom_record(line, line_num)
                elif record_type == "CONECT":
                    self._process_conect_record(line, line_num)
                elif record_type == "TER":
                    # Keep TER records in protein file
                    self.protein_atoms.append(line)
                else:
                    # Keep other records (TITLE, REMARK, etc.) in protein file
                    self.other_records.append(line)
    
    def _process_atom_record(self, line: str, line_num: int):
        """Process an ATOM or HETATM record."""
        try:
            # Parse atom record
            atom_num = int(line[6:11])
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain_id = line[21]
            res_num = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            
            # Determine if this is part of the ligand
            if res_name == self.ligand_name:
                self.ligand_atoms.append(line)
                self.ligand_atom_numbers.add(atom_num)
                print(f"  Line {line_num}: Found ligand atom {atom_name} in {res_name}")
            else:
                self.protein_atoms.append(line)
                self.protein_atom_numbers.add(atom_num)
                
        except (ValueError, IndexError) as e:
            print(f"Warning: Could not parse atom record at line {line_num}: {e}")
            # Keep unparseable records in protein file
            self.protein_atoms.append(line)
    
    def _process_conect_record(self, line: str, line_num: int):
        """Process a CONECT record and assign to appropriate component."""
        try:
            # Parse CONECT record
            parts = line.split()
            if len(parts) < 2:
                print(f"Warning: Invalid CONECT record at line {line_num}")
                return
            
            # Get the main atom number
            main_atom = int(parts[1])
            
            # Check if this CONECT involves ligand atoms
            connected_atoms = [int(parts[i]) for i in range(2, len(parts)) if parts[i].isdigit()]
            all_atoms = [main_atom] + connected_atoms
            
            # Determine if this CONECT belongs to ligand or protein
            ligand_atoms_in_conect = [a for a in all_atoms if a in self.ligand_atom_numbers]
            protein_atoms_in_conect = [a for a in all_atoms if a in self.protein_atom_numbers]
            
            if ligand_atoms_in_conect and protein_atoms_in_conect:
                # Mixed CONECT - this shouldn't happen in a properly split complex
                print(f"Warning: Mixed CONECT record at line {line_num} - skipping")
                return
            elif ligand_atoms_in_conect:
                # This CONECT belongs to the ligand
                self.ligand_connections.append(line)
                print(f"  Line {line_num}: Found ligand CONECT for atom {main_atom}")
            else:
                # This CONECT belongs to the protein
                self.protein_connections.append(line)
                
        except (ValueError, IndexError) as e:
            print(f"Warning: Could not parse CONECT record at line {line_num}: {e}")
    
    def write_output_files(self):
        """Write the separated protein and ligand files."""
        # Write protein file
        protein_file = f"{self.output_prefix}_protein.pdb"
        with open(protein_file, 'w') as f:
            # Write other records first
            for record in self.other_records:
                f.write(record + '\n')
            
            # Write protein atoms
            for atom in self.protein_atoms:
                f.write(atom + '\n')
            
            # Write protein CONECT records
            for conect in self.protein_connections:
                f.write(conect + '\n')
        
        print(f"Protein file written: {protein_file}")
        print(f"  Atoms: {len(self.protein_atoms)}")
        print(f"  CONECT records: {len(self.protein_connections)}")
        
        # Write ligand file
        ligand_file = f"{self.output_prefix}_ligand.pdb"
        with open(ligand_file, 'w') as f:
            # Write ligand atoms
            for atom in self.ligand_atoms:
                f.write(atom + '\n')
            
            # Write ligand CONECT records
            for conect in self.ligand_connections:
                f.write(conect + '\n')
        
        print(f"Ligand file written: {ligand_file}")
        print(f"  Atoms: {len(self.ligand_atoms)}")
        print(f"  CONECT records: {len(self.ligand_connections)}")
        
        return protein_file, ligand_file
    
    def run(self) -> Tuple[str, str]:
        """Run the complete splitting process."""
        self.parse_pdb()
        
        if not self.ligand_atoms:
            raise ValueError(f"No atoms found for ligand '{self.ligand_name}' in the PDB file")
        
        if not self.protein_atoms:
            raise ValueError("No protein atoms found in the PDB file")
        
        return self.write_output_files()


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Split a protein-ligand complex PDB file into separate components",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python split_complex.py complex.pdb JZ4 output
  python split_complex.py protein_ligand.pdb ATP results
        """
    )
    
    parser.add_argument("pdb_file", help="Input PDB file containing protein-ligand complex")
    parser.add_argument("ligand_name", help="3-letter residue name of the ligand (e.g., JZ4, ATP)")
    parser.add_argument("output_prefix", help="Prefix for output files (will create _protein.pdb and _ligand.pdb)")
    
    args = parser.parse_args()
    
    try:
        splitter = PDBComplexSplitter(args.pdb_file, args.ligand_name, args.output_prefix)
        protein_file, ligand_file = splitter.run()
        
        print("\n" + "="*50)
        print("SPLITTING COMPLETED SUCCESSFULLY")
        print("="*50)
        print(f"Input file: {args.pdb_file}")
        print(f"Ligand name: {args.ligand_name}")
        print(f"Output files:")
        print(f"  Protein: {protein_file}")
        print(f"  Ligand: {ligand_file}")
        print("="*50)
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main() 