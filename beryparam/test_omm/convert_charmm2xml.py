import parmed as pmd


# Define your CHARMM files


charmm_ff = pmd.charmm.CharmmParameterSet(
 "../forcefield/toppar_c36_jul20/top_all36_prot.rtf",
 "../forcefield/toppar_c36_jul20/par_all36m_prot.prm",
 "../forcefield/toppar_c36_jul20/top_all36_lipid.rtf",
 "../forcefield/toppar_c36_jul20/par_all36_lipid.prm",
 "../forcefield/toppar_c36_jul20/top_all36_na.rtf",
 "../forcefield/toppar_c36_jul20/par_all36_na.prm",
 "../forcefield/toppar_c36_jul20/top_all36_carb.rtf",
 "../forcefield/toppar_c36_jul20/par_all36_carb.prm",
 "../forcefield/toppar_c36_jul20/top_all36_cgenff.rtf",
 "../forcefield/toppar_c36_jul20/par_all36_cgenff.prm",
 "../forcefield/toppar_c36_jul20/toppar_water_ions_charmm.str",
 "../forcefield/toppar_c36_jul20/stream/lipid/toppar_all36_lipid_model_yalun.str")

print(charmm_ff)
charmm_ff.write(top="top.rtf", par = "par.prm") #, str="stream.str")

# Write to OpenMM XML format
output_file = "charmm36_jul20.xml"

try:
    # ParmEd provides a `write` function specific to OpenMM XML
    openmm_ff = pmd.openmm.OpenMMParameterSet.from_parameterset(charmm_ff)
    openmm_ff.write(output_file)
    print(f"CHARMM force field successfully converted to OpenMM XML: {output_file}")
except Exception as e:
    print(f"Error converting to OpenMM XML: {e}")


