import parmed as pmd


# Define your CHARMM files
rtf_file = "forcefield/toppar_c36_jul20/top_all36_prot.rtf"  # Replace with your topology file
prm_file = "forcefield/toppar_c36_jul20/par_all36_prot.prm"  # Replace with your parameter file
str_file = "forcefield/toppar_c36_jul20/toppar_water_ions_charmm.str"  # Replace with your parameter file

# Load CHARMM files
#rtf_list =' '
#prm_list =' '
#str_list =' '
#
#rtf_list+=' "../forcefield/toppar_c36_jul20/top_all36_prot.rtf", '
#prm_list+=' "../forcefield/toppar_c36_jul20/par_all36m_prot.prm", '
#rtf_list+=' "../forcefield/toppar_c36_jul20/top_all36_lipid.rtf",  '
#prm_list+=' "../forcefield/toppar_c36_jul20/par_all36_lipid.prm",  '
#rtf_list+=' "../forcefield/toppar_c36_jul20/top_all36_na.rtf",  '
#prm_list+=' "../forcefield/toppar_c36_jul20/par_all36_na.prm",  '
#rtf_list+=' "../forcefield/toppar_c36_jul20/top_all36_carb.rtf",  '
#prm_list+=' "../forcefield/toppar_c36_jul20/par_all36_carb.prm",  '
#rtf_list+=' "../forcefield/toppar_c36_jul20/top_all36_cgenff.rtf", '
#prm_list+=' "../forcefield/toppar_c36_jul20/par_all36_cgenff.prm", '
#str_list+=' "../forcefield/toppar_c36_jul20/toppar_water_ions_charmm.str",  '
#str_list+=' "../forcefield/toppar_c36_jul20/stream/lipid/toppar_all36_lipid_model_yalun.str" '


charmm_ff = pmd.charmm.CharmmParameterSet(
 "forcefield/toppar_c36_jul20/top_all36_prot.rtf",
 "forcefield/toppar_c36_jul20/par_all36m_prot.prm",
 "forcefield/toppar_c36_jul20/top_all36_lipid.rtf",
 "forcefield/toppar_c36_jul20/par_all36_lipid.prm",
 "forcefield/toppar_c36_jul20/top_all36_na.rtf",
 "forcefield/toppar_c36_jul20/par_all36_na.prm",
 "forcefield/toppar_c36_jul20/top_all36_carb.rtf",
 "forcefield/toppar_c36_jul20/par_all36_carb.prm",
 "forcefield/toppar_c36_jul20/top_all36_cgenff.rtf",
 "forcefield/toppar_c36_jul20/par_all36_cgenff.prm",
 "forcefield/toppar_c36_jul20/toppar_water_ions_charmm.str",
 "forcefield/toppar_c36_jul20/stream/lipid/toppar_all36_lipid_model_yalun.str")

print(charmm_ff)

# Write to OpenMM XML format
output_file = "charmm36.xml"

try:
    # ParmEd provides a `write` function specific to OpenMM XML
    openmm_ff = pmd.openmm.OpenMMParameterSet.from_parameterset(charmm_ff)
    openmm_ff.write(output_file)
    print(f"CHARMM force field successfully converted to OpenMM XML: {output_file}")
except Exception as e:
    print(f"Error converting to OpenMM XML: {e}")

#charmm_ff.write(output_file, format="openmm")
#print(f"CHARMM force field converted to OpenMM XML: {output_file}")

