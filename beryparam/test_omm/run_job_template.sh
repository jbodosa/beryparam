#!/bin/bash
#SBATCH -A ACCOUNT_NAME
#SBATCH -t TIME
#SBATCH --partition=gpu
#SBATCH -n 1
#SBATCH --gpus=a100:1
#SBATCH --job-name=JOBNAME
#SBATCH --output="run.out"
#
##SBATCH --gpus=a100_1g.5gb:1

source ~/.bashrc.conda
conda activate Dock-MD-FEP

## Running yank
echo "Running FEP on yank"

## MPI mode
build_mpirun_configfile "yank script --yaml=config.yaml"
mpiexec.hydra -f hostfile -configfile configfile

# Serial mode
#yank script --yaml=config.yaml

# Analyze the data
echo "Analyzing data..."
yank analyze --store=run/experiments

num=$(cat expno.txt)
num=$(($num + 1))
mv run run_$num
mv run.out run_$num.out
cp config.yaml config_$num.yaml

#python get_fkhist.py $num

rm expno.txt
echo $num >expno.txt
