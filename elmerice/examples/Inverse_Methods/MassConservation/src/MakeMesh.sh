#!/bin/bash
#########################################
# Bash script to generate the Elmer mesh
#  from rectangle .geo
#
# INPUTS:
#	- SRC_DIR: environment variable with the path to this directory
#	- Arguments:
#		- $1 : element size
#		- $2 : output mesh name (.msh) (OPTIONAL)
#########################################
# check that there is one input argument
if [ -z "$1" ]
then
    echo "No argument supplied"
    echo "provide required mesh resolution as argument"
    exit
fi

# argument is the required mesh resolution
res=$1

# if a second argument is present it is the gmsh mesh file
#  default is rectangle.msh
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
