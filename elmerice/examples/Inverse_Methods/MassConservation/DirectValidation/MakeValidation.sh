#!/bin/bash
 
## define the src dir
export SRC_DIR=../src

## compile required USFs
make -f $SRC_DIR/Makefile RAMP

## loop over several resolutions
res=(10000.0 5000.0 2000.0 1000.0 500.0)

for i in "${res[@]}"
do
	name=$(echo $i/1 | bc);

	## mesh the mesh (linear elements) and run the test
	. $SRC_DIR/MakeMesh.sh $i rectangle_"$name".msh
	ElmerSolver RAMP.sif -ipar 1 $name

	## increase element order and run the test
	ElmerGrid 2 2 rectangle_"$name" -increase -out rectangle_"$name"_2nd
	ElmerSolver RAMP_2nd.sif -ipar 1 $name
done

## sort the Save Scalars results as a function of the number of dofs
cat RAMP_*.dat | sort -k1 > Convergence.dat
cat RAMP2nd_*.dat | sort -k1 > Convergence2nd.dat

