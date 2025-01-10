
from openmm.app import *
from openmm import *
from openmm.unit import *
import parmed as pmd

pdb_input = "../meso/meso_vmd.pdb"

forcefield = ForceField("charmm36_jul20.xml")

pdb = PDBFile(pdb_input)
modeller = Modeller(pdb.topology, pdb.positions)

# Rename water residues to TIP3 (optional for CHARMM compatibility)
for residue in modeller.topology.residues():
    if residue.name == "MGL":
        residue.name = "MGLYOL"

modeller.addSolvent(forcefield, model="tip3p", boxSize=Vec3(3.0, 3.0, 3.0)*nanometers)
## ... Call some modelling functions here ...
#system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME)
#simulation = Simulation(modeller.topology, system, integrator)
#simulation.context.setPositions(modeller.positions)

## Manually map atom types to CHARMM types
#for atom in structure.atoms:
#    if atom.residue.name == "HOH":
#        if atom.name == "H1" or atom.name == "H2":
#            atom.type = "HT"  # CHARMM hydrogen type for TIP3
#        elif atom.name == "O":
#            atom.type = "OT"  # CHARMM oxygen type for TIP3

PDBFile.writeFile(modeller.topology, modeller.positions, open('output.pdb', 'w'))
