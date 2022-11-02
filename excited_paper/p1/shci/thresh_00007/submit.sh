#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH -t 60:00:00
#SBATCH -p largemem_q
#SBATCH -A nmayhall_group
#SBATCH --job-name p1_shci_7
#SBATCH --exclusive

module reset
module load gcc/8.2.0
module load openmpi/gcc/64/4.0.4
export OMP_NUM_THREADS=1


cd $SLURM_SUBMIT_DIR

# run job
mpirun -n 24 /home/vibin1/tinkercliff_shci/shci_excited/shci/shci

exit;
