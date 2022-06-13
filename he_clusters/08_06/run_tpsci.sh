#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -t 02:00:00
#SBATCH -p normal_q
#SBATCH -A nmayhall_group
##SBATCH --exclusive
#SBATCH --job-name 08_06
FILE=tpsci

# to see a list of available modules.
module load Julia/1.5.1-linux-x86_64
module load Anaconda3/2020.11
source activate tpsci

# Change to the directory from which the job was submitted
cd $SLURM_SUBMIT_DIR

date
julia --project=../../../FermiCG/.  -t 16  $FILE.jl  > $FILE.out
date

exit;
