#!/bin/bash
#PBS -l nodes=16:ppn=16
#PBS -l walltime=24:00:00
#PBS -m abe
#PBS -M changgoo@princeton.edu
#
# Calculate the number of processors
#PBS -o TI_disk_512.log
#PBS -e TI_disk_512.err
#

NPROCS=`wc -l < $PBS_NODEFILE`
JOBID=N512_kp4_n8_K0

cd /scratch/network/changgoo/TI_disk

echo Time is `date`
echo Directory is `pwd`
echo Number of Processors is $NPROCS

module load fftw
module load openmpi

mpirun -np $NPROCS /home/changgoo/athena_StarParticle/bin/athena -i /home/changgoo/athena_StarParticle/bin/athinput.starpar_ti job/problem_id=$JOBID domain1/Nx1=512 domain1/Nx2=512 domain1/Nx3=32 domain1/AutoWithNProc=$NPROCS domain1/x1min=-256 domain1/x1max=256 domain1/x2min=-256 domain1/x2max=256
