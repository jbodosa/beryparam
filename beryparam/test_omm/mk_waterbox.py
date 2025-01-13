import numpy as np
from openmm.app import *
from openmm import *
from openmm.unit import *
import parmed as pmd
from openmm.app import Element

# Read in the pdb
#pdb_input = "../dmp/dmp_c.pdb"
pdb_input = "../meso/meso_charmm.pdb"
# Read in the psf
#psf_input = "../dmp/dmp_convert.psf"
psf_input = "../meso/meso_charmm.psf"

# Read in the convreted charmm forcefield
forcefield = ForceField("charmm36_jul20.xml")

# Load the pdb and psf files
pdb = PDBFile(pdb_input)
psf = CharmmPsfFile(psf_input)
psf.setBox(3* nanometers, 3* nanometers, 3* nanometers, 90.0, 90.0, 90.0 )

# Create the model with correct resnames from psf
modeller = Modeller(psf.topology, pdb.positions)

# Create OpenMM system
system = forcefield.createSystem(
    topology= modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.2 * nanometers,
    switchDistance=1.0 * nanometers,
    constraints=HBonds
)

nonbonded = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]

charges = []
for i in range(system.getNumParticles()):
    charge, sigma, epsilon = nonbonded.getParticleParameters(i)
    charges.append(charge)

sys_charge = np.round(sum(charges).value_in_unit(elementary_charge), decimals=1)
## TEST
#sys_charge = 1.0

def add_inverted_flat_bottom_restraint(system=system, atom1=5, atom2=14):
    # Define the custom inverted flat-bottom restraint parameters
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

    # Add the restraint force to the system
    system.addForce(inverted_flat_bottom_restraint)

    return(system)

def add_ion(modeller=modeller, ref_index=6, atomic_num=19):
    # Add an ion at a defined distance from ref atom
    reference_atom_index = ref_index
    positions = modeller.getPositions()
    ref_position = positions[reference_atom_index]

    distance = 2.0 # nanometers

    # Define the normalized direction vector
    direction = np.array([1.0, 1.0, 1.0])
    direction /= np.linalg.norm(direction)
    direction = Vec3(direction[0], direction[1], direction[2])

    # Calculate ion's new position
    ion_position = ref_position.__add__(direction * distance *nanometers )
    #print("ion pos ", ion_position )

    # Get the element by atomic number
    ion_element = Element.getByAtomicNumber(atomic_num) #getBySymbol('K')
    element_name = ion_element.name

    #    if element_name == "potassium":
    #        ion_segid = "POT"
    #
    #    if element_name == "chlorine":
    #        ion_segid = "CLL"

    new_topology = Topology()
    ion_chain = new_topology.addChain("ION")
    ion_residue = new_topology.addResidue("ION", ion_chain)
    ion_atom = new_topology.addAtom("ION", ion_element, ion_residue)

    # Append the ion's position to the list of existing positions
    modeller.add(new_topology, [ion_position ])
    return(modeller)

if sys_charge ==0:
    pass
elif sys_charge <0:
# Add potassium ion
    num_particles = system.getNumParticles()
    print(num_particles)
    modeller = add_ion(modeller, ref_index=5, atomic_num=19)
    system = add_inverted_flat_bottom_restraint(system, atom1=6, atom2=14 )
    system.getNumConstraints()

elif sys_charge >0:
# Add chloride ions
    num_particles = system.getNumParticles()
    print(num_particles)
    modeller = add_ion(modeller, ref_index=5, atomic_num=17)
    system = add_inverted_flat_bottom_restraint(system, atom1=6, atom2=14 )
    system.getNumConstraints()

modeller.addSolvent(forcefield, model="tip3p", boxSize=Vec3(3.0, 3.0, 3.0)*nanometers, neutralize=False)
PDBFile.writeFile(modeller.topology, modeller.positions, open('mol_waterbox.pdb', 'w'))

# Save the system with the restraint to an XML file
with open("mol_waterbox.xml", "w") as f:
    f.write(XmlSerializer.serialize(system))

#
### Manually map atom types to CHARMM types
##for atom in modeller.topology.atoms():
##    if atom.residue.name == "HOH":
##        atom.residue.name = "TIP3"
##        if atom.name == "H1" or atom.name == "H2":
##            atom.type = "HT"  # CHARMM hydrogen type for TIP3
##        elif atom.name == "O":
##            atom.type = "OT"  # CHARMM oxygen type for TIP3
##
## Manually map atom types to CHARMM types
##for atom in modeller.topology.atoms():
#    #    if atom.residue.name == "MGL":
#    #        atom.residue.name = "MGLYOL"

#
#nonbonded = [f for f in system_res.getForces() if isinstance(f, NonbondedForce)][0]
#
##charge, sigma, epsilon = nonbonded.getParticleParameters(0)
##
##print(charge)
##print(epsilon)
##print(sigma)
#
##bonded = [f for f in system_res.getForces() if isinstance(f, HarmonicBondForce)][0]
##for i in range(bonded.getNumBonds()):
##    particle1, particle2, length, k = bonded.getBondParameters(i)
##    if particle1 == 0 or particle2 == 0:
##        print(f'Particles ({particle1}, {particle2}), length = {length}, k = {k}')
