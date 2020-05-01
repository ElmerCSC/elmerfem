#!/bin/bash

# make mesh
ElmerGrid 1 2 mesh2D -metis 2
# compile required USFs
make

# Parameters
lambda='0.0e00 1.0e03 1.0e04 1.0e05 5.0e05 1.0e06 5.0e06 1.0e07 5.0e07 1.0e08 1.0e09 1.0e10 1.0e11'
DATAFILE="..\/DATA\/MacAyeal_VELOCITIES_NOISE.txt"

rm -rf LCurve.dat

# looop over reg. coefs.
c=0
for i in $lambda
do
  c=$((c+1))

  echo $i
  # get .sif file
  NAME=OPT_"$c"
  sed  "s/<Lambda>/"$i"/g;s/<NAME>/$NAME/g;s/<OBS_FILE>/$DATAFILE/g" SIF/OPTIM.sif > OPTIM_$c.sif

  # run 
  echo OPTIM_$c.sif > ELMERSOLVER_STARTINFO
  mpirun -np 2 ElmerSolver_mpi

  # psot process
  python ../SCRIPTS/MakeReport.py $NAME
  echo $(tail -n 1 Cost_"$NAME".dat | awk '{print $3}') $(tail -n 1 CostReg_"$NAME".dat | awk '{print $2}') $i >> LCurve.dat 
done
