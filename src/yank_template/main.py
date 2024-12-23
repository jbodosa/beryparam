#!/usr/bin/env python3

import logging
from system import *

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
    #system.read_pdb("../meso/meso_charmm.pdb")
    system.read_pdb("../dmp/dmp.pdb")
    system.convert_pdb2crd("dmp_convert.crd", use_CHARMM=False)
    #   system.read_crd("../meso/meso_wrong.crd")
    #   system.read_psf("../meso/meso_wrong.psf")
    #system.write_pdb("dmp_c.pdb")
    #   system.write_psf("meso_wrong.psf")
    #   system.write_crd("meso_wrong.crd")

if __name__ == '__main__':
    main()

