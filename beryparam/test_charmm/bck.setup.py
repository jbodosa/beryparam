from openmm.app import *
from openmm import unit, Vec3
import numpy as np

# Load your PDB file
pdb_filename = '../meso/meso.pdb'
pdb = PDBFile(pdb_filename)

# Load the system
#forcefield = ForceField('amber14-all.xml')  # or other applicable force field
#system = forcefield.createSystem(pdb.topology, nonbondedMethod='PME', nonbondedCutoff=1.0*unit.nanometers)

# Get the coordinates of an atom, e.g., residue 10
target_atom_index = 6  # The index of the reference atom
target_position = pdb.positions[target_atom_index]

# Define the desired distance (e.g., 5 angstroms)
distance = 20.0 * unit.angstroms

# Calculate a new position for the ion along a vector
direction = np.random.randn(3)  # Randomly choose a vector direction
direction /= np.linalg.norm(direction)  # Normalize the vector
ion_position = target_position + distance * Vec3(*direction)

print(f'Placing ion at {ion_position}')

new_chain = pdb.topology.addChain()
new_residue = pdb.topology.addResidue('POT', new_chain)
# Manually create a new atom representing Na+
ion_atom = pdb.topology.addAtom( 'POT', Element.getBySymbol('K'), residue=new_residue)
#pdb.topology.addAtom(ion_atom)
pdb.positions.append(ion_position)

### vacuum system
output_filename = 'vacuum.pdb'
PDBFile.writeFile(pdb.topology, pdb.positions, open(output_filename, 'w'))

a = 6# nm
b = 6# nm
c = 6# nm
# Load CHARMM files
#psf = app.CharmmPsfFile("sod_cla_wat.psf")  # Path to your CHARMM .psf file
#pdb = app.PDBFile("sod_cla_wat_eq_last.pdb")        # Path to your CHARMM .pdb file
#pdb.topology.setPeriodicBoxVectors(Vec3(a, b, c))
ff = CharmmParameterSet("../forcefield/toppar_c36_jul20/top_all36_cgenff.rtf",
                        "../forcefield/toppar_c36_jul20/par_all36_cgenff.prm",
                        "../forcefield/toppar_c36_jul20/toppar_water_ions_charmm.str")

modeller = Modeller(pdb.topology, pdb.positions)
modeller.addSolvent(ff, model='tip3p', boxSize=Vec3(a, b, c)*unit.nanometers)
#system = ff.createSystem(modeller.topology, nonbondedMethod=PME)
### water system
output_filename = 'water.pdb'
PDBFile.writeFile(modeller.topology, pdb.positions, open(output_filename, 'w'))
