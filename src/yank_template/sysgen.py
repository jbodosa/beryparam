#!/usr/bin/env python3
#

# Imports
import os
import subprocess
import shutil
import logging

from openmm.app import *
from utils import *

# Remote zaratan
PACKMOL="/home/jbodosa/scratch/data/exec/packmol/packmol"
CHARMM="/home/jbodosa/scratch/data/exec/gnu/charmm"

## Local M1
#PACKMOL="/Users/jbodosa/Documents/UMD/Rotation/Lab/Work/data/packages/packmol-20.14.4-docs1/packmol"
#CHARMM="/Users/jbodosa/Documents/UMD/Rotation/Lab/Work/data/packages/gnu/charmm"

####################
# Logger #
####################

# create logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# create console handler and set level to debug
ch = logging.FileHandler('setup.log')
ch.setLevel(logging.DEBUG)

# create formatter
formatter = logging.Formatter(fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%m-%d-%Y %H:%M:%S >')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)


class sysgen:
    def __init__(self):
        """
            Read charge input file and output the charge.
        """
        logger.info("sysgen started ...")

    ###############
    ###############
    #   READERS     #
    ###############
    ###############

    # Read the PDB file generated by any software
    def read_pdb(self, pdb_file):
        """
          Read the pdb file into openmm.
          :param pdb_file: path to the pdb file
          :param box_dim: (a,b,c) tuple of box edges
          :param box_dim_unit: Angstrom or nanometer unit
        """
        #self.pdb_file = pdb_file
        #self.box_dim = box_dim
        #self.box_dim_unit = box_dim_unit

        try:
            pdb = PDBFile(pdb_file)
            self.pdb = pdb
            return(pdb)
        except FileNotFoundError:
            return("{file} does not exist.".format(file=pdb_file))
        except Exception as e:
            return("Got error : {str(e)}".format(e=e))

    # Read the CRD file generated by any software
    def read_crd(self, crd_file):
        """
          Read the crd file into openmm.
          :param crd_file: path to the crd file
        """

        try:
            crd = charmmcrdfiles.CharmmCrdFile(crd_file)
            self.crd = crd
            return f"CRD file read from '{crd_file}'"
        except FileNotFoundError:
            return("{file} does not exist.".format(file=crd_file))
        except Exception as e:
            return("Got error : {str(e)}".format(e=e))

    # Read the PSD file generated by any software
    def read_psf(self, psf_file):
        """
          Read the psf file into openmm.
          :param psf_file: path to the psf file
        """

        try:
            psf = charmmcrdfiles.CharmmPsfFile(psf_file)
            self.psf = psf
            return f"PSF file read from '{psf_file}'"
        except FileNotFoundError:
            return("{file} does not exist.".format(file=psf_file))
        except Exception as e:
            return("Got error : {str(e)}".format(e=e))

    # Create a CHARMM crd file from any pdb
    def write_crd(self, output_path, use_CHARMM=False):
        """
            Read the pdb file and write out the crd file.
        """
        try:
            if not use_CHARMM:
                pdb = self.pdb
                if isinstance(pdb, PDBFile):
                    with open(output_path, 'w') as crd_file:
                        # Write the CRD file header
                        crd_file.write(f"* CHARMM Coordinate File\n")
                        crd_file.write(f"* Generated by sysgen\n")
                        crd_file.write(f"*\n")
                        # "NUMATOMS_EXT": "{0:10d} EXT\n",
                        crd_file.write(f"{len(list(pdb.topology.atoms())):10d}  EXT\n")

                        # Write the atomic data
                        for i, atom in enumerate(pdb.topology.atoms()):
                            position = pdb.positions[i]
                            #segid = getattr(atom.residue, 'id', 'MGLYOL')

                            # Picked the format from MDAnalysis
                            # https://docs.mdanalysis.org/1.1.0/_modules/MDAnalysis/coordinates/CRD.html#CRDWriter
                            #
                            # "ATOM_EXT": ("{serial:10d}{totRes:10d}  {resname:<8.8s}  {name:<8.8s}"
                            # "{pos[0]:20.10f}{pos[1]:20.10f}{pos[2]:20.10f}  "
                            # "{chainID:<8.8s}  {resSeq:<8d}{tempfactor:20.10f}\n"),
                            #
                            # The resname and segname are wrong
                            atom.residue.name="MGLYOL"
                            atom.residue.id="MGLYOL"
                            crd_file.write(f"{i + 1:10d}{atom.residue.index + 1:10d}  {atom.residue.name:<8.8s}  {atom.name:<8.8s}{position.x:>20.10f}{position.y:20.10f}{position.z:20.10f}  {atom.residue.id:<8.8s}  {atom.residue.index + 1:<8d}{0:20.10f}\n")
                crd_file.close()
                return f"CRD file written to '{output_path}'"
        except FileNotFoundError:
            return("{file} does not exist.".format(file=pdb_file))
        except Exception as e:
            return("Got error : {str(e)}".format(e=e))

        # Create a CHARMM psf file from a crd
    def write_psf(self, output_path):
        """
            Write out the psf file for the given crd file.
        """

        try:
            crd = self.crd
            if isinstance(crd, charmmcrdfiles.CharmmCrdFile):
                try:
                    # Call write recenter.inp
                    recenter_out = write_recenter_inp()
                    #print(recenter_out)
                    logger.debug(recenter_out)
                    # Call CHARMM to recenter the sys and write psf
                    command = CHARMM+" < recenter.inp > recenter.out"
                    result = subprocess.run(command, shell=True, executable="/bin/bash",  capture_output=True, text=True, check=True)
                    return(result.stdout)
                except subprocess.CalledProcessError as e:
                    #print(f"Command failed with return code {e.returncode}")
                    logger.debug(f"Command failed with return code {e.returncode}")
                    #print(f"Error output: {e.stderr}")
                    logger.debug(f"Error output: {e.stderr}")
        except FileNotFoundError:
            return("{file} does not exist.".format(file=crd_file))
        except Exception as e:
            return("Got error : {str(e)}".format(e=e))


    ###########################
    # Add params to system #
    ###########################

    def set_param(self, sys_param):
        param_dict = parse_param(sys_param)
        for key, value in param_dict.items():
#            print(key, value)
            setattr(self, key, value)
        return(param_dict)

    ###########################
    # For neutral No counter-ions using packmol #
    ###########################
    ## Probably change this to CHARMM later
    # CLEAN

    def no_ions(self):
        pass

        ## Make a copy of ion.pdb and rename it correctly
        #replace_resname(self.counter_ion, "ion/ion.pdb", f"{self.counter_ion}.pdb")
        ## Write the packmol input file
        #pack_status = write_pack_inp(abs(int(self.ncharge)), self.counter_ion, self.ion_dist ) # How many counter ions and the dist

        #command = PACKMOL+" < box_ion.inp > box_ion.out"
        #result = subprocess.run(command, shell=True, executable="/bin/bash",  capture_output=True, text=True, check=True)
        #print(result.stdout)
        #return(f"Added {abs(int(self.ncharge))} number of {self.counter_ion} to the system.")


    ###########################
    # Add counter-ions using packmol #
    ###########################
    ## Probably change this to CHARMM later

    def place_ions(self):

        # Make a copy of ion.pdb and rename it correctly
        replace_resname(self.counter_ion, "ion/ion.pdb", f"{self.counter_ion}.pdb")
        # Write the packmol input file
        pack_status = write_pack_inp(abs(int(self.ncharge)), self.counter_ion, self.ion_dist ) # How many counter ions and the dist

        command = PACKMOL+" < box_ion.inp > box_ion.out"
        result = subprocess.run(command, shell=True, executable="/bin/bash",  capture_output=True, text=True, check=True)
        #print(result.stdout)
        logger.debug(result.stdout)
        # FIX
        # Need to recenter first
        self.read_pdb(pdb_file = 'box_ion.pdb') #, box_dim=60, box_dim_unit="Ang")
        self.write_crd('box_ion.crd', use_CHARMM=False)
        self.read_crd(crd_file = 'box_ion.crd')
        self.write_psf('box_ion.psf')
        return(f"Added {abs(int(self.ncharge))} number of {self.counter_ion} to the system.")


    def add_ions(self, ion_dist): # dist: distance to place away from sys/mol

        # OVER_HERE
        self.ion_dist = ion_dist
        #ncharge = int(self.ncharge)
        # TEST
        # Test line below
        self.ncharge = 1 # Test different charges
        self.ncharge = int(self.ncharge) # Need it to be int here

        if self.ncharge == 0:
            #print("Neutral system")
            logger.info("Neutral system")
            # TODO
            # Handle neutral system no counter-ion
            shutil.copy('recenter.pdb', 'box_ion.pdb')
            shutil.copy('recenter.crd', 'box_ion.crd')
            shutil.copy('recenter.psf', 'box_ion.psf')
        elif self.ncharge < 0 : # Negative charge
            #print(f"System has {self.ncharge} charge")
            logger.info(f"System has {self.ncharge} charge")
            # Add positive counter-ions POT/SOD
            self.counter_ion = "POT"
            ion_status = self.place_ions()
        elif self.ncharge > 0: # Positive charge
            #print(f"System has {self.ncharge} charge")
            logger.info(f"System has {self.ncharge} charge")
            # Add negative counter-ions POT/SOD
            self.counter_ion = "CLA"
            ion_status = self.place_ions()


        # Should return the box_ion.pdb
        return("Added ions")


# Example usage:


# Clean up first
file_path = ['recenter.pdb', 'recenter.psf', 'recenter.crd', 'recenter_off.pdb', 'box_ion.pdb', 'box_ion.psf', 'box_crd', 'POT.pdb', 'CLA.pdb']

for f in file_path:
    if os.path.exists(f):
        os.remove(f)
        #print(f"File '{f}' deleted successfully.")
        logger.info(f"File '{f}' deleted successfully.")

system = sysgen()
pdb = system.read_pdb(pdb_file = '../meso/meso.pdb') #, box_dim=60, box_dim_unit="Ang")
crd_status = system.write_crd('input.crd', use_CHARMM=False)
crd = system.read_crd(crd_file = 'input.crd')
psf = system.write_psf('input.psf')

param_dict = system.set_param("sys_param.str")
#print(system.__dict__)
system.add_ions(ion_dist =2)
#print("Done")
logger.info("Done")

