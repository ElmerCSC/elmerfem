#!/bin/bash
# check that there is one input argument
if [ -z "$1" ]
then
    echo "No argument supplied"
    echo "provide required mesh resolution as argument"
    exit
fi

# argument is the required mesh resolution
res=$1
# convert to int for the output name
printf -v int %.0f "$res"

if [ -z "$2" ]
then
  output=rectangle.msh
else
  output=$2
fi

# make gmsh mesh
gmsh -format msh2 -2 -setnumber lc $res $SRC_DIR/rectangle.geo -o $output

# convert to Elmer
ElmerGrid 14 2 $output -autoclean
