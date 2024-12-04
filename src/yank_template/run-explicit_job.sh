#!/bin/bash
#SBATCH -A energybio-eng
#SBATCH -t 10:00:00
#SBATCH --partition=gpu
#SBATCH -n 1
#SBATCH --gpus=a100:2
#SBATCH --job-name=JNAME
#SBATCH --output="explicit.out"
#
##SBATCH --gpus=a100_1g.5gb:1

source ~/.bashrc.conda
conda activate Dock-MD-FEP

# Run the simulation
echo "Running simulation..."
python mk_vacuum.py 
python mk_water.py 

## MPI mode
build_mpirun_configfile "yank script --yaml=explicit.yaml"

mpiexec.hydra -f hostfile -configfile configfile
# Serial mode
#yank script --yaml=explicit.yaml

date
# Analyze the data
echo "Analyzing data..."
yank analyze --store=explicit/experiments

num=$(cat expno.txt)
num=$(($num + 1))
mv explicit explicit_$num
mv explicit.out explicit_$num.out
cp explicit.yaml explicit_$num.yaml

python get_fkhist.py $num

rm expno.txt
echo $num >expno.txt
