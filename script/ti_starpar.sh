#!/bin/bash
#PBS -l nodes=4:ppn=16
#PBS -l walltime=24:00:00
#PBS -m abe
#PBS -M changgoo@princeton.edu
#
# Calculate the number of processors
#PBS -o TI_disk_256.log
#PBS -e TI_disk_256.err
#

NPROCS=`wc -l < $PBS_NODEFILE`
JOBID=N256_kp4_n8_K0

cd /scratch/network/changgoo/TI_disk

echo Time is `date`
echo Directory is `pwd`
echo Number of Processors is $NPROCS

module load fftw
module load openmpi

mpirun -np $NPROCS /home/changgoo/athena_StarParticle/bin/athena -i /home/changgoo/athena_StarParticle/bin/athinput.starpar_ti job/problem_id=$JOBID domain1/Nx1=256 domain1/Nx2=256 domain1/Nx3=64 domain1/AutoWithNProc=$NPROCS
