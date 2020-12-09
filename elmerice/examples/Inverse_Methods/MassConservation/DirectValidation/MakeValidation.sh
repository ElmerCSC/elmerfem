#!/bin/bash
 
export SRC_DIR=../src

## compile
make -f $SRC_DIR/Makefile RAMP

res=(10000.0 5000.0 2000.0 1000.0 500.0)

for i in "${res[@]}"
do
	name=$(echo $i/1 | bc);
	. $SRC_DIR/MakeMesh.sh $i rectangle_"$name".msh
	ElmerSolver RAMP.sif -ipar 1 $name

	ElmerGrid 2 2 rectangle_"$name" -increase -out rectangle_"$name"_2nd
	ElmerSolver RAMP_2nd.sif -ipar 1 $name
done

cat RAMP_*.dat | sort -k1 > Convergence.dat
cat RAMP2nd_*.dat | sort -k1 > Convergence2nd.dat

