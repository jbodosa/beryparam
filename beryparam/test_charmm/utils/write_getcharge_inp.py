import os
from datetime import datetime, timezone
import logging


logger = logging.getLogger(__name__)

ION_LIST = ["SOD", "POT", "CAL","CLA"]
#####################################
# CHARMM inp file to write psf #
#####################################

### Writing a CHARMM input file
### Write the mol psf script
def write_getcharge_inp(crd_infile, psf_infile):
    output_file = "ions.inp"

    with open(output_file, 'w') as inp:
        # Get username
        user = os.getenv('USERNAME') or os.getenv('USER')
        # Get current UTC date and time
        current_utc = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M")

        inp.write(f"* GENERATED BY SYSGEN\n")
        inp.write(f"* THIS IS TO WRITE THE PSF FOR SYS \n")
        inp.write(f"* Made by {user} on {current_utc}\n")
        inp.write(f"* \n")
        inp.write(f"\n")
        inp.write(f"DIMENS CHSIZE 2000000\n")
        #        inp.write(f"!DIMENS CHSIZE 5000000 MAXRES 3000000\n")
        inp.write(f"\n")
        inp.write(f"bomlev -3 \n")
        inp.write(f"\n")

        inp.write(f"! Read charmm ff files \n")
        inp.write(f"stream toppar.str \n")
        inp.write(f"\n")

        #        for mol_segid, mol_n in  zip(mol_segid_list, mol_n_list):
        #            if mol_segid not in ION_LIST:
        #                # If not ion
        #                inp.write(f"! read the sequence\n")
        #                inp.write(f"read sequence {mol_segid} {mol_n}\n")
        #                inp.write(f"generate {mol_segid}\n")
        #                inp.write(f"\n")
        #            else:
        #                # If it is ion
        #                inp.write(f"! read the sequence\n")
        #                inp.write(f"read sequence {mol_segid} {mol_n}\n")
        #                inp.write(f"generate {mol_segid} noangle nodihedral\n")
        #                inp.write(f"\n")

        inp.write("! Read the connex from the PSF file\n")
        inp.write(f"open read unit 30 card name {psf_infile}\n")
        inp.write("read psf card unit 30 \n")
        inp.write("close unit 30\n")
        inp.write("\n")

        inp.write("! Read the coordinates from the CRD file\n")
        inp.write(f"open read unit 31 card name {crd_infile}\n")
        inp.write("read coor card unit 31 \n")
        inp.write("close unit 31\n")
        inp.write("\n")
        inp.write("coor stat sele all end \n")
        inp.write("calc cgtot = int ( ?cgtot ) \n")
        inp.write("\n")

        inp.write("open write unit 90 card name psf_param.str \n")
        inp.write("write title unit 90 \n")
        inp.write("* int ncharge = @cgtot \n")
        inp.write("* float xmax = ?xmax \n")
        inp.write("* float ymax = ?ymax \n")
        inp.write("* float zmax = ?zmax \n")
        inp.write("* float xmin = ?xmin \n")
        inp.write("* float ymin = ?ymin \n")
        inp.write("* float zmin = ?zmin \n")
        inp.write("* \n")
        inp.write("\n")
        inp.write("stop")
        inp.close()
        logger.debug("Wrote ions inp")
    return(f"Wrote inp for ions")

