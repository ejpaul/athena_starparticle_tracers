#!/bin/bash
#PBS -l nodes=16:ppn=16
#PBS -l walltime=24:00:00
#PBS -m abe
#PBS -M changgoo@princeton.edu
#
# Calculate the number of processors
#PBS -o IS1_L256.log
#PBS -e IS1_L256.err
#

NPROCS=`wc -l < $PBS_NODEFILE`
JOBID=L256_rst

cd /scratch/network/changgoo/IS1_e2_n4_v4_rhosd01

echo Time is `date`
echo Directory is `pwd`
echo Number of Processors is $NPROCS

module load fftw
module load openmpi

#mpirun -np $NPROCS /home/changgoo/athena_StarParticle/bin/athena -i /home/changgoo/athena_StarParticle/bin/athinput.starpar_ti_tiger job/problem_id=$JOBID domain1/Nx1=256 domain1/Nx2=256 domain1/Nx3=64 domain1/AutoWithNProc=$NPROCS domain1/x1min=-128 domain1/x1max=128 domain1/x2min=-128 domain1/x2max=128
mpirun -np $NPROCS /home/changgoo/athena_StarParticle/bin/athena -r rst/L256.0002.rst job/problem_id=$JOBID output3/dt=10 output3/time=110 log/iflush=1
