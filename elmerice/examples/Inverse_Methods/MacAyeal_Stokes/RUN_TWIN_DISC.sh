#!/bin/bash
# make mesh
if [ ! -d "mesh2D" ]; then
  ElmerGrid 1 2 mesh2D
fi
# compile required USFs
make

#parameters
lambda='0.0e00'
DATAFILE="..\/DATA\/MacAyeal_VELOCITIES.txt"
NAME="TWINC"

# get .sif file
sed  "s/<Lambda>/"$lambda"/g;s/<NAME>/$NAME/g;s/<OBS_FILE>/$DATAFILE/g" SIF/OPTIM_DISC.sif > OPTIM_$NAME.sif

# run 
echo OPTIM_$NAME.sif > ELMERSOLVER_STARTINFO
ElmerSolver

# post process
python ../SCRIPTS/MakeReport.py $NAME
