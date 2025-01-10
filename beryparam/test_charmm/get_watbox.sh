#!/bin/bash
#

PACKMOL="/home/jbodosa/scratch/data/exec/packmol/packmol"
CHARMM="/home/jbodosa/scratch/data/exec/gnu/charmm"

$PACKMOL <ion_box.inp >ion_box.out

$CHARMM <ion_wat.inp >ion_wat.out
$CHARMM <ion_vac.inp >ion_vac.out
