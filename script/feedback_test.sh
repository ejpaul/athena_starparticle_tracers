#!/bin/bash
#PBS -l nodes=4:ppn=16
#PBS -l walltime=1:00:00
#PBS -m abe
#PBS -M changgoo@princeton.edu
#
# Calculate the number of processors
#PBS -o n1_SN_2pc_cool_ST.log
#PBS -e n1_SN_2pc_cool_ST.err
#

NPROCS=`wc -l < $PBS_NODEFILE`
JOBID=n1_SN_2pc_cool_ST_v00
RUNDIR=/scratch/network/changgoo/feedback_test_new
SOURCEDIR=/home/changgoo/athena_starparticle/
BIN=$SOURCEDIR/bin/athena
INPUT=$SOURCEDIR/athinput/athinput.starpar_feedback_test

if ! mkdir -p $RUNDIR; then
	echo "Cannot create $RUNDIR"
	exit 1;
fi;

cd $RUNDIR

echo Time is `date`
echo Directory is `pwd`
echo Number of Processors is $NPROCS

module load fftw
module load openmpi

mpirun -np $NPROCS $BIN -i $INPUT job/problem_id=$JOBID domain1/AutoWithNProc=$NPROCS  problem/amp=0.0 problem/n0=1 problem/n1=1 problem/rin=1000. problem/rout=1000. problem/tHII=-1 problem/tSN=0.0 problem/ton=0.0 log/file_open=0 domain1/x1min=-64 domain1/x1max=64 domain1/x2min=-64 domain1/x2max=64 domain1/x3min=-64 domain1/x3max=64 domain1/Nx1=64   domain1/Nx2=64   domain1/Nx3=64   time/tlim=2.0
