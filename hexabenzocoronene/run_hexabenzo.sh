#!/bin/bash

#SBATCH --nodes=1
#SBATCH -p preemptable_q
##SBATCH -p normal_q
#SBATCH --cpus-per-task=16
#SBATCH -t 01:00:00
#SBATCH --account=nmayhall_group
#SBATCH --job-name hexabenzo
##SBATCH --exclusive

# . /etc/bashrc


if [ -z ${HOME+x} ];
then
    export HOME=$(echo ~)
    source /etc/profile
    source /etc/bashrc
    source $HOME/.bashrc
fi

hostname
sleep 10

module reset
module load site/tinkercliffs-rome/easybuild/setup  #only for infer
module load site/tinkercliffs/easybuild/setup  #only for infer
module load Anaconda3/2020.07
module load gcc/8.2.0

source activate fermi_python

export INFILE=$1
export OUTFILE="${INFILE}.out"
export WORKDIR=$(pwd)

echo $INFILE
echo $OUTFILE
echo $WORKDIR
echo $TMPDIR

cp $INFILE $TMPDIR/
if [ "$2" ]
then
        cp $2 $TMPDIR/
fi
if [ "$3" ]
then
        cp $3 $TMPDIR/
fi
cd $TMPDIR


#Start an rsync command which runs in the background and keeps a local version of the output file up to date
touch $OUTFILE
while true; do rsync -av $OUTFILE $WORKDIR/"${INFILE}.${SLURM_JOB_ID}.out"; sleep 60; done &

python $INFILE >& $OUTFILE

cp $OUTFILE $WORKDIR/"${INFILE}.out"
rm $WORKDIR/"${INFILE}.${SLURM_JOB_ID}.out"

mkdir $WORKDIR/"${INFILE}.${SLURM_JOB_ID}.scr"

cp -r * $WORKDIR/"${INFILE}.${SLURM_JOB_ID}.scr"
rm -r *

#moving standard output slurm file to specific job directory
mv $WORKDIR/"slurm-${SLURM_JOB_ID}.out" $WORKDIR/"${INFILE}.${SLURM_JOB_ID}.scr"

exit

