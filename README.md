# Running MMGBSA

I used OpenMM and Ambertools to run MMGBSA on a protein-ligand complex. This repository is the product of my testing with these tools. 

Use ambertools to generate `.prmtop` and `.inpcrd` files using the tleap file. You may need to visit the H++ webserver to assign charges to the PDB. I used `obabel` to get mol and mol2 files from the ligand PDB. Use `prepareFiles.py` to simulate the simulation files, and then run the gbscript file with `MMPBSA.py` from ambertools to perform MMGBSA. It will output the calculated binding affinity of the ligand to the receptor protein, as well as a decomposition of contributions from each residue in the receptor. The calculation is based on the Generalized Born Surface Area, which assumes an implicit solvent.   

