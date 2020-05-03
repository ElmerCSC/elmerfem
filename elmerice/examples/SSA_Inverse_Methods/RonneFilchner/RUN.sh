#!/bin/bash
## compile required USFs
make USFs_RonneFilchner

## Data file
DATAFILE="..\/DATA\/RonneFilchner.nc"

## regularisation parameters
lambda='0.0e00 1.0e03 1.0e04 2.0e04 5.0e04 1.0e05 2.0e05 5.0e05 1.0e06 2.0e06 5.0e06 1.0e07 1.0e08 1.0e09 1.0e10 1.0e11'
rm -rf LCurve.dat

c=0
for i in $lambda
do
  c=$((c+1))

  echo $i
  NAME=OPT_"$c"
  sed  "s/<Lambda>/"$i"/g;s/<NAME>/$NAME/g;s/<OBS_FILE>/$DATAFILE/g" SIF/OPTIM_ETA.sif > OPTIM_$c.sif

  echo OPTIM_$c.sif > ELMERSOLVER_STARTINFO
  # Has to be parallel on 2 partition to restart initial file
  mpirun -np 2 ElmerSolver_mpi

  python ../SCRIPTS/MakeReport.py $NAME
  echo $(tail -n 1 Cost_"$NAME".dat | awk '{print $3}') $(tail -n 1 CostReg_"$NAME".dat | awk '{print $2}') $i $c >> LCurve.dat
done
~
