#!/bin/bash
#PBS -l nodes=2:ppn=16
#PBS -l walltime=24:00:00
#PBS -m abe
#PBS -M changgoo@princeton.edu
#
# Calculate the number of processors
#PBS -o IS1_128.log
#PBS -e IS1_128.err
#

NPROCS=`wc -l < $PBS_NODEFILE`
JOBID=N128

cd /scratch/network/changgoo/IS1_e2_n4_v4_rhosd01

echo Time is `date`
echo Directory is `pwd`
echo Number of Processors is $NPROCS

module load fftw
module load openmpi

#mpirun -np $NPROCS /home/changgoo/athena_StarParticle/bin/athena -i /home/changgoo/athena_StarParticle/bin/athinput.starpar_ti_tiger job/problem_id=$JOBID domain1/Nx1=128 domain1/Nx2=128 domain1/Nx3=64 domain1/AutoWithNProc=$NPROCS 
mpirun -np $NPROCS /home/changgoo/athena_StarParticle/bin/athena -r rst/N128.0004.rst
