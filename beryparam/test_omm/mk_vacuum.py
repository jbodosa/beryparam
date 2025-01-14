
from openmm import CustomBondForce, System, XmlSerializer
from openmm.app import Simulation
from openmm import app
from yank.experiment import ExperimentBuilder
from openmmtools import states
import openmm.unit as unit

a=6#nm
b=6#nm
c=6#nm
# Load CHARMM files
psf = app.CharmmPsfFile("sod_cla_vac.psf")  # Path to your CHARMM .psf file
pdb = app.PDBFile("sod_cla_vac.pdb")        # Path to your CHARMM .pdb file
psf.setBox(a*unit.nanometers, b*unit.nanometers, c*unit.nanometers)
params = app.CharmmParameterSet("../forcefield/toppar_c36_jul20/top_all36_cgenff.rtf",
                                "../forcefield/toppar_c36_jul20/par_all36_cgenff.prm",
                                "../forcefield/toppar_c36_jul20/toppar_water_ions_charmm.str")
# CHARMM parameter files

# Define the custom inverted flat-bottom restraint parameters
atom1 = 0  # Index of the first atom in the PSF file
atom2 = 1  # Index of the second atom in the PSF file
threshold = 2.0 * unit.nanometers  # Distance threshold
spring_constant = 5.0 * unit.kilocalories_per_mole / unit.angstrom**2

# Define the inverted flat-bottom potential (apply if r < threshold)
inverted_flat_bottom_restraint = CustomBondForce(
    "step(threshold - r) * 0.5 * k * (threshold - r)^2"
)
inverted_flat_bottom_restraint.addPerBondParameter("threshold")
inverted_flat_bottom_restraint.addPerBondParameter("k")

# Add the restraint between atom1 and atom2
inverted_flat_bottom_restraint.addBond(atom1, atom2, [threshold, spring_constant])

# Create OpenMM system
system = psf.createSystem(
    params,
    nonbondedMethod=app.NoCutoff, #app.PME,
    nonbondedCutoff=1.2 * unit.nanometers,
    switchDistance=1.0 * unit.nanometers,
    constraints=app.HBonds,
)

# Add the restraint force to the system
system.addForce(inverted_flat_bottom_restraint)

# Save the system with the restraint to an XML file
with open("vacuum_with_restraint.xml", "w") as f:
    f.write(XmlSerializer.serialize(system))

