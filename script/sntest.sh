#!/bin/bash

xmin=0
xmax=128
Nx=4
for rin in 32
do
	id=n1_SN_Lx${xmax}_Nx${Nx}_rin${rin}_TE_h_oct
	./athena_TEh -i athinput.sntest_SNoct -d sntest_TI job/problem_id=$id problem/rin=$rin problem/rout=$rin log/file_open=1 domain1/Nx1=${Nx} domain1/Nx2=${Nx} domain1/Nx3=${Nx} domain1/AutoWithNProc=1 problem/heat_ratio=0.5 problem/fm=0.0 job/maxout=1 output2/dt=0.01 output1/dt=1.e-4 domain1/x1min=$xmin domain1/x1max=$xmax domain1/x2min=$xmin domain1/x2max=$xmax domain1/x3min=$xmin domain1/x3max=$xmax
done
