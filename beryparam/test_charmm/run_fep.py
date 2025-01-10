
from openmm import CustomBondForce, System
from openmm.app import Simulation
from yank.experiment import ExperimentBuilder
from openmmtools import states
import openmm.unit as unit

# Step 1: Define the inverted flat-bottom restraint force template
def create_inverted_flat_bottom(particle1, particle2, r_in=2, k=500.0):
    force = CustomBondForce(
        """
        step(r_in - r) * 0.5 * k * (r_in - r)^2;
        r = distance(p1, p2);
        """
    )
    force.addGlobalParameter("r_in", r_in)  # Inner radius (nm)
    force.addGlobalParameter("k", k)  # Spring constant (kJ/mol/nm^2)
    force.addBond(particle1, particle2, [])  # Add specific particle indices
    return force

# Step 2: Create systems with specific restraints for water and vacuum
#
temperature= 298*unit.kelvin
pressure= 1*unit.atmosphere

res_cla1_wat = create_inverted_flat_bottom(0, 20509) #Chloride
system_water = System()
system_water.addForce(res_cla1_wat)  # Restrain 1-20510

#thermodynamic_state_water = states.ThermodynamicState(system=system_water, temperature=temperature, pressure=pressure)
#corr_wat_res1 = res_cla1_wat.get_standard_state_correction(thermodynamic_state_water)
#corr_wat_res2 = res_cla2_wat.get_standard_state_correction(thermodynamic_state_water)

#print(f"Standard State Correction (Water) Res1: {corr_wat_res1} kJ/mol") 
#print(f"Standard State Correction (Water) Res2: {corr_wat_res2} kJ/mol") 
#
res_cla1_vac = create_inverted_flat_bottom(0, 20509)
system_vacuum = System()
system_vacuum.addForce(res_cla1_vac)  # Restrain 1-20510

#thermodynamic_state_vac = states.ThermodynamicState(system=system_vacuum, temperature=temperature, pressure=pressure)
#corr_vac_res1 = res_cla1_vac.get_standard_state_correction(thermodynamic_state_vac)
#corr_vac_res2 = res_cla2_vac.get_standard_state_correction(thermodynamic_state_vac)
#
#print(f"Standard State Correction (vacuum) Res1: {corr_vac_res1} kJ/mol") 
#print(f"Standard State Correction (vacuum) Res2: {corr_vac_res2} kJ/mol") 

# Step 3: Load the provided YAML content
yaml_content = """
options:
  minimize: yes
  verbose: yes
  output_dir: explicit
  default_number_of_iterations: 50
  temperature: 298*kelvin
  pressure: 1*atmosphere
  checkpoint_interval: 10
  #alchemical_pme_treatment: exact

solvents:
  water:
    nonbonded_method: PME
    nonbonded_cutoff: 25*angstroms
    clearance: 16*angstroms
    negative_ion: CLA
  vacuum:
    nonbonded_method: NoCutoff #PME

systems:
  hydration-system:
    phase1_path: [sod_cla_wat.pdb, sod_cla_wat.psf]
    phase2_path: [sod_cla_vac.pdb, sod_cla_vac.psf]
    solvent1: water
    solvent2: vacuum
    solvent_dsl: not resname SOD
    charmm_parameter_files: [../toppar_c36_jul20/top_all36_cgenff.rtf,
                             ../toppar_c36_jul20/par_all36_cgenff.prm,
                             ../toppar_c36_jul20/toppar_water_ions_charmm.str]

protocols:
  hydration-protocol:
    solvent1:
      alchemical_path:
        lambda_electrostatics: [1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00]
        lambda_sterics:        [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]
        #lambda_electrostatics: [1.00, 0.9999, 0.9998, 0.9997, 0.9996, 0.9995, 0.9994, 0.9993, 0.9992, 0.9991, 0.9990, 0.9989, 0.9988, 0.9987, 0.9986, 0.9985, 0.9984, 0.9983, 0.9982, 0.9981, 0.9980, 0.997, 0.996, 0.995, 0.994, 0.993, 0.992, 0.991, 0.98, 0.97, 0.96, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00]
        #lambda_sterics:        [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]
    solvent2:
      alchemical_path:
        lambda_electrostatics: [1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00]
        lambda_sterics:        [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]
        #lambda_electrostatics: [1.00, 0.9999, 0.9998, 0.9997, 0.9996, 0.9995, 0.9994, 0.9993, 0.9992, 0.9991, 0.9990, 0.9989, 0.9988, 0.9987, 0.9986, 0.9985, 0.9984, 0.9983, 0.9982, 0.9981, 0.9980, 0.997, 0.996, 0.995, 0.994, 0.993, 0.992, 0.991, 0.98, 0.97, 0.96, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00]
        #lambda_sterics:        [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]
        #lambda_electrostatics: [1.00, 0.9999, 0.9998, 0.9997, 0.9996, 0.9995, 0.9994, 0.9993, 0.9992, 0.9991, 0.9990, 0.998, 0.997, 0.996, 0.995, 0.994, 0.993, 0.992, 0.991, 0.98, 0.97, 0.96, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00]
        #lambda_sterics:        [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]
        #lambda_electrostatics: [1.00, 0.999, 0.998, 0.997, 0.996, 0.995, 0.994, 0.993, 0.992, 0.991, 0.98, 0.97, 0.96, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00]
        #lambda_sterics:        [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]
        #lambda_electrostatics: [1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00]
        #lambda_sterics:        [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]


experiments:
  system: hydration-system
  protocol: hydration-protocol
"""
# Save the YAML content to a file
with open("explicit.yaml", "w") as f:
    f.write(yaml_content)

# Step 4: Initialize and run the experiment using ExperimentBuilder
experiment = ExperimentBuilder('explicit.yaml')
experiment.run_experiments()  # Run the Yank experiments

## Step 4: Build the Yank experiment
#experiment = ExperimentBuilder(yaml_content)
#experiment.setup()  # Initialize Yank experiment
#
## Step 5: Attach the appropriate system to the simulation based on the phase
#simulation_water = Simulation(system_water, None, None)  # For water phase
#simulation_vacuum = Simulation(system_vacuum, None, None)  # For vacuum phase
#
