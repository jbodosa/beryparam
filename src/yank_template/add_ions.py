#!/usr/bin/env python3
#

## Is the molecule/system neutral ?
# If not then add ions
# Else exit with simple no
## Run packmol script to add ions
#First check the charge of molecule/system
# Run in subprocess $PACKMOL <ion_box.inp >ion_box.out


#```
#    Check the syscharge file in the structre dir. If neutral molecule then don't add charge else check if positive
#    or negative and add correct counter ion (CLA/POT) or ask for counter-ion name.
#```

# Imports
import os
import subprocess
from openmm.app import *
# Read the input file with charge
charge_file = "../meso/syscharge.dat"

class SysGen:
    def __init__(self, pdb_file, box_dim, box_dim_unit):
        """
            Read charge input file and output the charge.
            :param pdb_file: path to the pdb file
            :param box_dim: (a,b,c) tuple of box edges
            :param box_dim_unit: Angstrom or nanometer unit
        """
        self.pdb_file = pdb_file
        self.box_dim = box_dim
        self.box_dim_unit = box_dim_unit

    ### Module to read and output the charge from charge_path
    ##def get_charge(self):
    ##    try:
    ##        with open(self.file_path, 'r') as f:
    ##            charge=f.read()
    ##        return(charge)
    ##    except FileNotFoundError:
    ##        return("{file} does not exit.".format(file=self.file_path))
    ##    except Exception as e:
    ##        return("Got error : {str(e)}".format(e=e))

    # Read the PDB file generated by any software
    def read_pdb(self):
        """
          Read the pdb file into openmm.
        """
        try:
            pdb = PDBFile(self.pdb_file)
            return(pdb)
        except FileNotFoundError:
            return("{file} does not exist.".format(file=self.pdb_file))
        except Exception as e:
            return("Got error : {str(e)}".format(e=e))

    def write_crd(self, output_path):
        """
            Read the pdb file and write out the crd file.
        """
        try:
            pdb = self.read_pdb()
            if isinstance(pdb, PDBFile):
                with open(output_path, 'w') as crd_file:
                    # Write the CRD file header
                    crd_file.write(f"* CHARMM Coordinate File\n")
                    crd_file.write(f"* Generated by SysGen\n")
                    crd_file.write(f"*\n")
                    # "NUMATOMS_EXT": "{0:10d} EXT\n",
                    crd_file.write(f"{len(list(pdb.topology.atoms())):10d}  EXT\n")

                    # Write the atomic data
                    for i, atom in enumerate(pdb.topology.atoms()):
                        position = pdb.positions[i]
                        segid = getattr(atom.residue, 'id', 'HETATOM')

                        # Picked the format from MDAnalysis
                        # https://docs.mdanalysis.org/1.1.0/_modules/MDAnalysis/coordinates/CRD.html#CRDWriter
                        #
                        # "ATOM_EXT": ("{serial:10d}{totRes:10d}  {resname:<8.8s}  {name:<8.8s}"
                        # "{pos[0]:20.10f}{pos[1]:20.10f}{pos[2]:20.10f}  "
                        # "{chainID:<8.8s}  {resSeq:<8d}{tempfactor:20.10f}\n"),
                        #
                        # The resname and segname are wrong
                        crd_file.write(f"{i + 1:10d}{atom.residue.index + 1:10d}  {atom.residue.name:<8.8s}  {atom.name:<8.8s}{position.x:>20.10f}{position.y:20.10f}{position.z:20.10f}  {segid:<8.8s}  {atom.residue.index + 1:<8d}{0:20.10f}\n")
            return f"CRD file written to '{output_path}'"
        except Exception as e:
            raise e

# Example usage:
file_reader = SysGen(pdb_file = '../meso/meso.pdb', box_dim=60, box_dim_unit="Ang")
pdb = file_reader.read_pdb()
crd_status = file_reader.write_crd('output.crd')
print(crd_status)

