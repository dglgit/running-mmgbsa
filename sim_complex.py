"""
Run a MD simulation for a complex, optionally adding a solvent box
"""

import sys, time, argparse

from openff.toolkit import Molecule
from openmmforcefields.generators import SystemGenerator
import openmm
from openmm import app, unit, LangevinIntegrator, Vec3
from openmm.app import PDBFile, Simulation, Modeller, PDBReporter, StateDataReporter, DCDReporter, ForceField
from pdbfixer import PDBFixer
from openmm import Platform
t0 = time.time()


parser = argparse.ArgumentParser(description="simulateComplex")

parser.add_argument("-p", "--protein", required=True, help="Protein PDB file")
parser.add_argument("-l", "--ligand", required=False, help="Ligand molfile")
parser.add_argument("-o", "--output", default='output', help="Base name for output files")
parser.add_argument("-s", "--steps", type=int, default=5e5, help="Number of steps")
parser.add_argument("-z", "--step-size", type=float, default=0.001, help="Step size (ps)")
parser.add_argument("-f", "--friction-coeff", type=float, default=1, help="Friction coefficient (ps)")
parser.add_argument("-i", "--interval", type=int, default=1000, help="Reporting interval")
parser.add_argument("-t", "--temperature", type=int, default=300, help="Temperature (K)")
parser.add_argument("--solvate", action='store_true', help="Add solvent box")
parser.add_argument("--padding", type=float, required=False, help="Padding for solvent box (A)")
parser.add_argument("--water-model", default="tip3p",
                    choices=["tip3p", "spce", "tip4pew", "tip5p", "swm4ndp"],
                    help="Water model for solvation")
parser.add_argument("--positive-ion", default="Na+", help="Positive ion for solvation")
parser.add_argument("--negative-ion", default="Cl-", help="Negative ion for solvation")
parser.add_argument("--ionic-strength", type=float, default="0", help="Ionic strength for solvation")
parser.add_argument("--no-neutralize", action='store_true', help="Don't add ions to neutralize")
parser.add_argument("-e", "--equilibration-steps", type=int, default=2e5, help="Number of equilibration steps")
parser.add_argument("--protein-force-field", default='amber/ff14SB.xml', help="Protein force field")
parser.add_argument("--ligand-force-field", default='gaff-2.11', help="Ligand force field")
parser.add_argument("--water-force-field", default='amber/tip3p_standard.xml', help="Ligand force field")
parser.add_argument("--membrane",action='store_true')
parser.add_argument("--num-solvent",type=int,required=False)

args = parser.parse_args()
print("simulateComplex: ", args)

pdb_in = args.protein
mol_in = args.ligand
output_base = args.output
output_complex = output_base + '_complex.pdb'
output_traj_dcd = output_base + '_traj.dcd'
output_min = output_base + '_minimised.pdb'
num_steps = args.steps
reporting_interval = args.interval
temperature = args.temperature * unit.kelvin
equilibration_steps = args.equilibration_steps
print('Processing', pdb_in, 'and', mol_in, 'with', num_steps, 'steps generating outputs',
      output_complex, output_min, output_traj_dcd)
if mol_in is None:
    print('no ligand used')
# get the chosen or fastest platform
platform = Platform.getPlatformByName('OpenCL')#utils.get_platform()
print(f'platform: {platform}')
#exit()

if mol_in is not None:
    print('Reading ligand')
    ligand_mol = Molecule.from_file(mol_in)

print('Preparing system')
# Initialize a SystemGenerator using the GAFF for the ligand and tip3p for the water.
forcefield_kwargs = {'constraints': app.HBonds, 'rigidWater': True, 'removeCMMotion': True, 'hydrogenMass': 4*unit.amu }
print('forcefield kwargs is:',forcefield_kwargs)
if mol_in is not None:
    system_generator = SystemGenerator(
        forcefields=[args.protein_force_field, args.water_force_field],
        small_molecule_forcefield=args.ligand_force_field,
        molecules=[ligand_mol],
        forcefield_kwargs=forcefield_kwargs)
else:
    system_generator = SystemGenerator(
        forcefields=[args.protein_force_field, args.water_force_field],
        small_molecule_forcefield=args.ligand_force_field,
        forcefield_kwargs=forcefield_kwargs)

# Use Modeller to combine the protein and ligand into a complex
print('Reading protein')
protein_pdb = PDBFixer(filename=pdb_in)
protein_pdb.findMissingResidues()
protein_pdb.findMissingAtoms()
protein_pdb.findNonstandardResidues()
protein_pdb.addMissingAtoms()
protein_pdb.addMissingHydrogens(7.4)
protein_pdb.removeHeterogens(False)

print('Preparing complex')
# The topology is described in the openforcefield API

modeller = Modeller(protein_pdb.topology, protein_pdb.positions)
modeller.addHydrogens(ForceField(args.protein_force_field))
print('System has %d atoms' % modeller.topology.getNumAtoms())

# The topology is described in the openforcefield API
if mol_in is not None:
    print('Adding ligand...')
    numLigands=1
    lig_top = ligand_mol.to_topology()
    for i in range(numLigands):
        modeller.add(lig_top.to_openmm(), lig_top.get_positions().to_openmm())
    print('added',numLigands,"of the ligand")
    print('System has %d atoms' % modeller.topology.getNumAtoms())
boxVectors=[Vec3(15,0,0),Vec3(0,15,0),Vec3(0,0,15)]
# Solvate
if args.solvate:
    print('Adding solvent...')
    # we use the 'padding' option to define the periodic box.
    # we just create a box that has a 10A (default) padding around the complex.
    if args.membrane:
        print('adding membrane')
        modeller.addMembrane(system_generator.forcefield, lipidType='POPC', minimumPadding=1*unit.nanometer,ignoreExternalBonds=True)
    else:
        print('adding water')
        if args.num_solvent is None:
            modeller.addSolvent(system_generator.forcefield, model=args.water_model, padding=args.padding * unit.angstroms,
                            positiveIon=args.positive_ion, negativeIon=args.negative_ion,
                            ionicStrength=args.ionic_strength * unit.molar, neutralize=not args.no_neutralize)
        else:
            print('adding',args.num_solvent,'waters')
            modeller.addSolvent(system_generator.forcefield, model=args.water_model, numAdded=args.num_solvent,
                            positiveIon=args.positive_ion, negativeIon=args.negative_ion,
                            ionicStrength=args.ionic_strength * unit.molar, neutralize=not args.no_neutralize)
    print('System has %d atoms' % modeller.topology.getNumAtoms())

with open(output_complex, 'w') as outfile:
    PDBFile.writeFile(modeller.topology, modeller.positions, outfile)

# Create the system using the SystemGenerator
if mol_in is not None:
    system = system_generator.create_system(modeller.topology, molecules=ligand_mol)
else:
    system = system_generator.create_system(modeller.topology)
friction_coeff = args.friction_coeff / unit.picosecond
step_size = args.step_size * unit.picoseconds
duration = (step_size * num_steps).value_in_unit(unit.nanoseconds)
print('Simulating for {} ns'.format(duration))

integrator = LangevinIntegrator(temperature, friction_coeff, step_size)


if system.usesPeriodicBoundaryConditions():
    #system.setDefaultPeriodicBoxVectors(*boxVectors) 
    print('Default Periodic box: {}'.format(system.getDefaultPeriodicBoxVectors()))
else:
    print('No Periodic Box')
if args.solvate:
    system.addForce(openmm.MonteCarloBarostat(1 * unit.atmospheres, temperature, 25))

addRestraint=False
'''
if addRestraint: 
    restraint = openmm.CustomCentroidBondForce(1, "0.5*k*((x1-x0)^2 + (y1-y0)^2 + (z1-z0)^2)")
    restraint.addPerBondParameter("x0")
    restraint.addPerBondParameter("y0")
    restraint.addPerBondParameter("z0")
    restraint.addGlobalParameter("k", 10000.0)  # kJ/mol/nm^2

    # Define the group (all protein atoms)
    protein_atoms = [atom.index for atom in modeller.topology.atoms() if atom.residue.name != 'HOH']
    restraint.addGroup(protein_atoms)

    # Add the bond to a fixed point (e.g., origin)
    restraint.addBond([0], [0, 0, 0])

    # Add the force to the system
    system.addForce(restraint)
'''
if addRestraint:
    restraint = openmm.CustomExternalForce("k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    restraint.addGlobalParameter("k", 10.0 * openmm.unit.kilojoules_per_mole/openmm.unit.nanometer**2)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    for atom in modeller.topology.atoms():
        if atom.name == "CA":  # Restrain C-alpha atoms
            restraint.addParticle(atom.index, modeller.positions[atom.index])
    system.addForce(restraint)
#for i in range(210, 230):
#    system.setParticleMass(i, 9999999)
system.addForce(openmm.CMMotionRemover(10))
simulation = Simulation(modeller.topology, system, integrator, platform=platform)
context = simulation.context
context.setPositions(modeller.positions)

print('Minimising ...')
simulation.minimizeEnergy()

# Write out the minimised PDB.
with open(output_min, 'w') as outfile:
    PDBFile.writeFile(modeller.topology, context.getState(getPositions=True, enforcePeriodicBox=True).getPositions(), file=outfile, keepIds=True)

# equilibrate
simulation.context.setVelocitiesToTemperature(temperature)
print('Equilibrating ...')
simulation.step(equilibration_steps)

# Run the simulation.
stateFile=open(output_base+"_stateData.txt",'w')
simulation.reporters.append(DCDReporter(output_traj_dcd, reporting_interval, enforcePeriodicBox=True))
simulation.reporters.append(StateDataReporter(sys.stdout, reporting_interval * 5, step=True, potentialEnergy=True, temperature=True, totalEnergy=True, kineticEnergy=True))
simulation.reporters.append(StateDataReporter(stateFile, reporting_interval * 5, step=True, potentialEnergy=True, temperature=True, totalEnergy=True, kineticEnergy=True))
print('Starting simulation with', num_steps, 'steps ...')
t1 = time.time()
simulation.step(num_steps)
t2 = time.time()
print('Simulation complete in {} mins at {}. Total wall clock time was {} mins'.format(
    round((t2 - t1) / 60, 3), temperature, round((t2 - t0) / 60, 3)))
print('Simulation time was', round(duration, 3), 'ns')
stateFile.close()