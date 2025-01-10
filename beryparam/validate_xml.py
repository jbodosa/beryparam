from openmm.app import *
from openmm import *
from openmm.unit import *

# Load the OpenMM XML file
forcefield = ForceField("charmm36.xml")

# Example: Create a water box
modeller = Modeller(Topology(), [])
modeller.addSolvent(forcefield, boxSize=Vec3(3.0, 3.0, 3.0) * nanometer)

# Create system and simulation
system = forcefield.createSystem(modeller.topology)
print("System successfully created with the converted CHARMM force field!")

