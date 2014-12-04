#!/bin/bash
#PBS -l nodes=8:ppn=16
#PBS -l walltime=24:00:00
#PBS -m abe
#PBS -M changgoo@princeton.edu
#
# Calculate the number of processors
#PBS -o TI_disk_L128_N256_rst.log
#PBS -e TI_disk_L128_N256_rst.err
#

NPROCS=`wc -l < $PBS_NODEFILE`
JOBID=N256_dx05_kp4_n8_K0

cd /scratch/network/changgoo/TI_disk_rst

echo Time is `date`
echo Directory is `pwd`
echo Number of Processors is $NPROCS

module load fftw
module load openmpi

mpirun -np $NPROCS /home/changgoo/athena_StarParticle/bin/athena -r N256_dx05_kp4_n8_K0.0002.rst output1/num=0
