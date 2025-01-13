
import numpy as np
from openmm.app import *
from openmm import *
from openmm.unit import *
from openmm import XmlSerializer


from sys import stdout


# Input and output file names
xml_file = "water_with_restraint.xml"
pdb_file = "output_after.pdb"

# Load the system from the XML file
with open(xml_file, 'r') as f:
    system = XmlSerializer.deserialize(f.read())

# Load the topology (necessary for the simulation setup)
pdb = PDBFile(pdb_file)

forcefield = ForceField("charmm36_jul20.xml")
integrator = LangevinIntegrator(300, 1, 0.001)  # Temp (K), Friction coeff (ps^-1), Step size (ps)

# Create the simulation
simulation = Simulation(pdb.topology, system, integrator)

# Set initial positions
simulation.context.setPositions(pdb.positions)

# Minimize energy (optional)
print("Minimizing energy...")
simulation.minimizeEnergy()

# Iterate through residues and print their names
print("Residue Names in the System:")
for residue in pdb.topology.residues():
    print(f"Residue Name: {residue.name}, Residue ID: {residue.id}")


# Write out the structure to PDB
out_file = "xml_test.pdb"
print(f"Writing structure to {out_file}...")
positions = simulation.context.getState(getPositions=True).getPositions()
with open(out_file, 'w') as f:
    PDBFile.writeFile(pdb.topology, positions, f)

print("Done!")

