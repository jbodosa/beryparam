
import MDAnalysis as mda

## Load the PDB file
#u = mda.Universe("recenter.pdb")
#
## Write the coordinates to a CRD file
#with mda.Writer("test.crd", n_atoms=u.atoms.n_atoms) as writer:
#    writer.write(u.atoms)
#

# Load the PDB file
u = mda.Universe("recenter.pdb")

# Guess elements for the atoms
u.add_TopologyAttr('elements')
print(u.atoms.elements)  # Check the guessed elements

## Write the coordinates to a CRD file
#with mda.Writer("test.crd", n_atoms=u.atoms.n_atoms) as writer:
#    writer.write(u.atoms)

# Write to a CHARMM-style CRD file
with mda.coordinates.CRD.CRDWriter("test.crd", n_atoms=u.atoms.n_atoms) as writer:
    writer.write(u.atoms)

