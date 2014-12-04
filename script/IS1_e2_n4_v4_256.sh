#!/bin/bash
#PBS -l nodes=32:ppn=16
#PBS -l walltime=72:00:00
#PBS -m abe
#PBS -M changgoo@princeton.edu
#
# Calculate the number of processors
#PBS -o IS1_256_SN.log
#PBS -e IS1_256_SN.err
#

NPROCS=`wc -l < $PBS_NODEFILE`
JOBID=N256_SN2
RUNDIR=/scratch/network/changgoo/IS1_e2_n4_v4_rhosd01

cd $RUNDIR
cp /home/changgoo/athena_StarParticle/bin/athena_$JOBID $RUNDIR/athena_$JOBID

echo Time is `date`
echo Directory is `pwd`
echo Number of Processors is $NPROCS

module load fftw
module load openmpi

#mpirun -np $NPROCS /home/changgoo/athena_StarParticle/bin/athena -i /home/changgoo/athena_StarParticle/bin/athinput.starpar_ti_tiger job/problem_id=$JOBID domain1/Nx1=256 domain1/Nx2=256 domain1/Nx3=256 domain1/AutoWithNProc=$NPROCS 
mpirun -np $NPROCS $RUNDIR/athena_$JOBID -r rst/N256_feedback.0002.rst job/problem_id=$JOBID
