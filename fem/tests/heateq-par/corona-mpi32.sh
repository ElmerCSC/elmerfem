#!/bin/csh
#$ -S /bin/csh
#$ -N mpi-test
#$ -o test.out
#$ -cwd
#$ -l h_rt=00:10:00
#$ -l h_vmem=200M
#$ -pe cre 2

setenv LD_LIBRARY_PATH ".:/opt/SUNWhpc/lib:/wrk/vierinen/src/buildtools/build.corona1.2005-06-09-10-49-12/fem/src/"

printf "JOB_ID: %s\n" $JOB_ID
printf "PE_HOSTFILE: %s\n" $PE_HOSTFILE

printf "Running it\n"
mprun -np 2 ../../src/ElmerSolver_mpi

printf "Done...\n" 
