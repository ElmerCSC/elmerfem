#!/bin/bash

# mesh resolution
res=2000.0

export SRC_DIR=../src

## Make mesh
. $SRC_DIR/MakeMesh.sh $res

## compile
make -f $SRC_DIR/Makefile

## get data
python3 $SRC_DIR/MakeObs.py

## run init
ElmerSolver INIT_OPT.sif

## run assimilation for different regularisation parameters
lambda=(0.001 0.1 0.5 1.0 2.5 5.0 10.0 50.0)

c=0
for i in "${lambda[@]}"
do
	c=$((c+1))
	ElmerSolver Optimisation.sif -ipar 1 $c -rpar 1 $i
done


. ./MakeLCurve.sh
