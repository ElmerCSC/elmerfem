#!/bin/bash
 
## define the src dir
export SRC_DIR=../src

## mesh resolution
res=2000.0

## Make mesh
. $SRC_DIR/MakeMesh.sh $res

## compile required USFs
make -f $SRC_DIR/Makefile

## generate synthetic obsevations
python3 $SRC_DIR/MakeObs.py

## run validation
for i in {1..10}
do
	ElmerSolver Validation.sif -ipar 1 $i
done


