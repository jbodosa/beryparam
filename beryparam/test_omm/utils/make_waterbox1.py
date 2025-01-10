from openmm.app import *
from openmm import *
from openmm.unit import *
import parmed as pmd

# Define simulation box size (in nanometers)
box_size = 3.0

# Load the force field (TIP3P water model is common for CHARMM)
forcefield = ForceField("../charmm36_jul20.xml")

# Create an empty topology and system
modeller = Modeller(Topology(), [])

# Add water to the empty box
modeller.addSolvent(forcefield, model="tip3p", boxSize=Vec3(box_size, box_size, box_size) * nanometer)

# Get the system from the modeller
system = forcefield.createSystem(modeller.topology)

# Print details about the water box
print(f"Number of atoms: {len(modeller.positions)}")
print(f"Box dimensions: {box_size}")

# Optionally save the water box to a PDB file
with open("waterbox.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
print("Water box saved to 'waterbox.pdb'")

system = forcefield.createSystem(modeller.topology)

# Convert the system to a ParmEd Structure
structure = pmd.openmm.load_topology(modeller.topology, system, xyz=modeller.positions)

# Manually map atom types to CHARMM types
for atom in structure.atoms:
    if atom.residue.name == "HOH":
        if atom.name == "H1" or atom.name == "H2":
            atom.type = "HT"  # CHARMM hydrogen type for TIP3
        elif atom.name == "O":
            atom.type = "OT"  # CHARMM oxygen type for TIP3

# Write CHARMM-compatible PSF and CRD files
structure.save("waterbox.psf", format="psf", overwrite=True)
structure.save("waterbox.crd", overwrite=True)

print("PSF file saved to 'waterbox.psf'")
print("CRD file saved to 'waterbox.crd'")
