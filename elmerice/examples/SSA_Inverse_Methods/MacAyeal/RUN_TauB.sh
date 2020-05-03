#!/bin/bash

## you can choose to run serial or parallel:
## here is the number of partitions
NP=1

# make mesh
ElmerGrid 1 2 mesh2D 
##
if [ $NP -gt 1 ]
then
  ElmerGrid 2 2 mesh2D -metis $NP
fi
# compile required USFs
make

# Parameters
lambda='0.0e00 1.0e03 1.0e04 1.0e05 5.0e05 1.0e06 5.0e06 1.0e07 5.0e07 1.0e08 1.0e09 1.0e10 1.0e11'
DATAFILE="..\/DATA\/MacAyeal_VELOCITIES_NOISE.txt"

LCURVE=LCurve_TB.dat
rm -rf $LCURVE

c=0
for i in $lambda
do
  c=$((c+1))

  echo $i
  NAME=OPT_TB_"$c"
  sed  "s/<Lambda>/"$i"/g;s/<NAME>/$NAME/g;s/<OBS_FILE>/$DATAFILE/g" SIF/OPTIM_TAUB.sif > OPTIM_TB_$c.sif

  echo OPTIM_TB_$c.sif > ELMERSOLVER_STARTINFO
  if [ $NP -gt 1 ]
  then
    mpirun -np 2 ElmerSolver_mpi
  else
    ElmerSolver
  fi

  python ../SCRIPTS/MakeReport.py $NAME
  echo $(tail -n 1 Cost_"$NAME".dat | awk '{print $3}') $(tail -n 1 CostReg_"$NAME".dat | awk '{print $2}') $i $c >> $LCURVE
done
