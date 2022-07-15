#!/bin/bash

#SBATCH --nodes=1
#SBATCH -p normal_q
#SBATCH --cpus-per-task=16
#SBATCH -t 10:00:00
#SBATCH --account=nmayhall_group
##SBATCH --account=nmayhall_group-paid
#SBATCH --job-name N24
#SBATCH --exclusive

#export NTHREAD=16
#export JULIAENV=/home/nbraunsc/FermiCG/


# allow time to load modules, this takes time sometimes for some reason
sleep 10

hostname

module reset
module load site/tinkercliffs-rome/easybuild/setup  #only for infer
module load site/tinkercliffs/easybuild/setup  #only for infer
#module load Python/3.8.6-GCCcore-10.2.0
#module load Julia/1.7.2-linux-x86_64
module load Anaconda3/2020.07
module load gcc/8.2.0

source activate mim_env


echo "Usage: sbatch submit.sh {input file} {data file} {data file}"

# set these to 1 if you are using julia threads heavily to avoid oversubscription
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

export INFILE=$1
export OUTFILE="${INFILE}.out"
export WORKDIR=$(pwd)

echo $INFILE
echo $OUTFILE
echo $WORKDIR
echo $TMPDIR


#julia --project=$JULIAENV -t $NTHREAD $INFILE >& $OUTFILE
python nat_orb_active_space.py > nat_orb_active_space.out

exit

