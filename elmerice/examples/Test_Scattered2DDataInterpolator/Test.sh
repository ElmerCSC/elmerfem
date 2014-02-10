#!/bin/bash

###
make clean

### make the mesh
ElmerGrid 1 2 mesh2D

### run the test
ElmerSolver case.sif
