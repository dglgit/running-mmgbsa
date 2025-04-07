from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout,argv

inpcrd = AmberInpcrdFile(argv[1])
prmtop = AmberPrmtopFile(argv[2], periodicBoxVectors=inpcrd.boxVectors)
method=PME if inpcrd.boxVectors is not None else NoCutoff
system = prmtop.createSystem(implicitSolvent=OBC2,nonbondedMethod=method, constraints=HBonds, implicitSolventSaltConc=0.15*moles/liter)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(inpcrd.positions)
with open(argv[3]+"_initial.pdb", 'w') as outfile:
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), file=outfile, keepIds=True)

print('starting minimization')
simulation.minimizeEnergy()
#simulation.reporters.append(PDBReporter(argv[1].split('.')[0]+'.csv', 1000))
#simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
#        potentialEnergy=True, temperature=True))
simulation.context.setVelocitiesToTemperature(300)
print('starting equilibration')
simulation.step(10000)
with open(argv[3]+'_eqed.pdb', 'w') as outfile:
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), file=outfile, keepIds=True)
output_base=argv[3].split('.')[0]
simulation.reporters.append(DCDReporter(output_base+"_traj.dcd", 1000))
print('starting production run')
simulation.step(2e6)