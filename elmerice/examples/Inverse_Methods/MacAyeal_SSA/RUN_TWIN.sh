#!/bin/bash
# make mesh
ElmerGrid 1 2 mesh2D 
# compile required USFs
make

#parameters
lambda='0.0e00'
DATAFILE="..\/DATA\/MacAyeal_VELOCITIES.txt"
NAME="TWIN"

# get .sif file
sed  "s/<Lambda>/"$lambda"/g;s/<NAME>/$NAME/g;s/<OBS_FILE>/$DATAFILE/g" SIF/OPTIM.sif > OPTIM_$NAME.sif

# run 
echo OPTIM_$NAME.sif > ELMERSOLVER_STARTINFO
ElmerSolver

# post process
python ../SCRIPTS/MakeReport.py $NAME
