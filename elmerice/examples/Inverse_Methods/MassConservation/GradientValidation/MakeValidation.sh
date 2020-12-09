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

## run validation
for i in {1..10}
do
	ElmerSolver Validation.sif -ipar 1 $i
done


