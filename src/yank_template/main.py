#!/usr/bin/env python3

import logging
from system import *
from utils.waterbox import *

# Clean up first
file_path = ['recenter.pdb', 'recenter.psf', 'recenter.crd', 'recenter_off.pdb', 'box_ion.pdb', 'box_ion.psf', 'box_crd', 'POT.pdb', 'CLA.pdb', "setup.log", "app.log"]

for f in file_path:
    if os.path.exists(f):
        os.remove(f)
        print(f"File '{f}' deleted successfully.")
        #logger.info(f"File '{f}' deleted successfully.")

def main():
    ####################
    ####################
    # Logger #
    ####################
    ####################

    # create logger
    logging.basicConfig(filename='setup.log', level=logging.DEBUG, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%m-%d-%Y %H:%M:%S >')
    #    logging.getLogger(__name__)
    #    logging.setLevel(logging.DEBUG)
    #
    #    # create console handler and set level to debug
    #    ch = logging.FileHandler('setup.log')
    #    ch.setLevel(logging.DEBUG)
    #
    #    # create formatter
    #    formatter = logging.Formatter(fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%m-%d-%Y %H:%M:%S >')
    #
    #    # add formatter to ch
    #    ch.setFormatter(formatter)
    #
    #    # add ch to logger
    #    #
    #    logging.addHandler(ch)
    # Example usage:

    logging.info("Start logging ...")
    system = System()
    system.read_pdb("../dmp/dmp.pdb")
    #system.convert_pdb2crd("dmp_convert.crd", use_CHARMM=False, resname_fix=True, resname_old_list=["MGL"], resname_new_list=["MGLYOL"])
    system.convert_pdb2crd("convert.crd", use_CHARMM=False, resname_fix=False, resname_old_list=[], resname_new_list=[])
    #system.read_crd("meso_convert.crd")
    system.recenter_crd(crd_infile="convert.crd")
    system.write_psf(psf_outfile="recenter.psf")
    system.add_ions()
    write_waterbox_inp(30,30,30)

if __name__ == '__main__':
    main()

