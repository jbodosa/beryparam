import numpy as np
from openmm.app import *
from openmm import *
from openmm.unit import *
import parmed as pmd

#from openmm import CustomBondForce, System, XmlSerializer

pdb_input = "../meso/meso_vmd.pdb"
psf_input = "meso_charmm.psf"

forcefield = ForceField("charmm36_jul20.xml")

pdb = PDBFile(pdb_input)
psf = CharmmPsfFile(psf_input)

modeller = Modeller(psf.topology, pdb.positions)

reference_atom_index = 6  # Adjust to the index of your reference atom
positions = pdb.getPositions(asNumpy=False)
ref_position = positions[reference_atom_index]

#print("ref pos ", ref_position )
# Define desired distance in Å
distance = 2.0 #* nanometers

# Define the normalized direction vector
direction = np.array([1.0, 1.0, 1.0])
direction /= np.linalg.norm(direction)
direction = Vec3(direction[0], direction[1], direction[2])

# Calculate ion's new position
ion_position = ref_position.__add__(direction * distance *nanometers )
#print("ion pos ", ion_position )

# Create a new PDB topology for the ion
from openmm.app import Element

ion_element = Element.getBySymbol('K')

## Manually add the ion to the topology and positions
## Manually update the topology to include the ion
new_topology = Topology()
#
#for chain in psf.topology.chains():
#    new_chain = new_topology.addChain()
#    for residue in chain.residues():
#        new_residue = new_topology.addResidue(residue.name, new_chain)
#        for atom in residue.atoms():
#            new_topology.addAtom(atom.name, atom.element, new_residue)


pot_chain = new_topology.addChain('POT')
pot_residue = new_topology.addResidue('POT', pot_chain)
pot_atom = new_topology.addAtom('POT', ion_element, pot_residue)

# Append the ion's position to the list of existing positions
modeller.add(new_topology, [ion_position ])
#print(modeller.topology)
#print(new_topology)
#print(new_positions)
## Save everything back to a new PDB file
#for pos in modeller.positions:
#    print(pos)

#for atom in modeller.topology.atoms():
#    print("element : " , atom.element,"symbol : " , atom.element.symbol)

modeller.addSolvent(forcefield, model="tip3p", boxSize=Vec3(3.0, 3.0, 3.0)*nanometers, neutralize=False)

PDBFile.writeFile(modeller.topology, modeller.positions, open('output_before.pdb', 'w'))
# Manually map atom types to CHARMM types
for atom in modeller.topology.atoms():
    if atom.residue.name == "HOH":
        atom.residue.name = "TIP3"
        if atom.name == "H1" or atom.name == "H2":
            atom.type = "HT"  # CHARMM hydrogen type for TIP3
        elif atom.name == "O":
            atom.type = "OT"  # CHARMM oxygen type for TIP3

# Manually map atom types to CHARMM types
for atom in modeller.topology.atoms():
    if atom.residue.name == "MGL":
        atom.residue.name = "MGLYOL"


    #with open('meso_with_ion.pdb', 'w') as f:
    #    PDBFile.writeFile(modeller.topology, modeller.positions, f)

PDBFile.writeFile(modeller.topology, modeller.positions, open('output_after.pdb', 'w'))
#
#system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME)
#pmd_structure = pmd.openmm.load_topology(modeller.topology, system, modeller.positions)
#pmd_structure.save("output_psf.psf", format="psf", overwrite=True)
#pmd_structure.save("output_crd.crd", overwrite=True)
#pmd_structure.save("output_pdb.pdb", overwrite=True)


# Define the custom inverted flat-bottom restraint parameters
atom1 = 5  # Index of the first atom in the PSF file
atom2 = 14  # Index of the second atom in the PSF file
threshold = 2.0 * nanometers  # Distance threshold
spring_constant = 5.0 * kilocalories_per_mole / angstrom**2

# Define the inverted flat-bottom potential (apply if r < threshold)
inverted_flat_bottom_restraint = CustomBondForce(
    "step(threshold - r) * 0.5 * k * (threshold - r)^2"
)
inverted_flat_bottom_restraint.addPerBondParameter("threshold")
inverted_flat_bottom_restraint.addPerBondParameter("k")

# Add the restraint between atom1 and atom2
inverted_flat_bottom_restraint.addBond(atom1, atom2, [threshold, spring_constant])
# Create OpenMM system
system_res = forcefield.createSystem(
    topology= modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.2 * nanometers,
    switchDistance=1.0 * nanometers,
    constraints=HBonds
)

# Add the restraint force to the system
system_res.addForce(inverted_flat_bottom_restraint)

# Save the system with the restraint to an XML file
with open("water_with_restraint.xml", "w") as f:
    f.write(XmlSerializer.serialize(system_res))

