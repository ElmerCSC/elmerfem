#!/bin/bash

# make mesh
if [ ! -d "mesh2D" ]; then
  ElmerGrid 1 2 mesh2D 
fi
# compile required USFs
make

#parameters
lambda='0.0e00'
NAME="TWIND"

# get .sif file
cp SIF/DIRECT.sif .
sed  "s/<Lambda>/"$lambda"/g;s/<NAME>/$NAME/g" SIF/OPTIM_CONT.sif > OPTIM_$NAME.sif

# run 
#echo DIRECT.sif > ELMERSOLVER_STARTINFO
#ElmerSolver

echo OPTIM_$NAME.sif > ELMERSOLVER_STARTINFO
ElmerSolver

# post process
python ../SCRIPTS/MakeReport.py $NAME
