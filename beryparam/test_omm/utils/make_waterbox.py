from openmm.app import *
from openmm import *
from openmm.unit import *
import parmed as pmd

# Parameters
box_size = 3.0 * nanometers  # Dimensions of the water box
water_model = "tip3p"  # Water model compatible with CHARMM36
output_psf = "water_box.psf"
output_crd = "water_box.crd"
output_pdb = "water_box.pdb"

# Load CHARMM36 force field
forcefield = ForceField("../charmm36_jul20.xml")

# Create a water box
modeller = Modeller(Topology(), [])
modeller.addSolvent(forcefield, model=water_model, boxSize=Vec3(box_size, box_size, box_size))

# Create a system
system = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff)

# Save to PSF and CRD using ParmEd
pmd_structure = pmd.openmm.load_topology(modeller.topology, system, modeller.positions)
pmd_structure.save(output_psf, overwrite=True)
pmd_structure.save(output_crd, overwrite=True)
pmd_structure.save(output_pdb, overwrite=True)

print(f"PSF file saved to {output_psf}")
print(f"CRD file saved to {output_crd}")
print(f"PDB file saved to {output_pdb}")

