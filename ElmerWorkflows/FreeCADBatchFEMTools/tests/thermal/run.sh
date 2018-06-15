#!/bin/bash

rm -r thermal thermal.sif thermal.unv scalars.dat* ELMERSOLVER_STARTINFO 2> /dev/null

FreeCAD -c $PWD/thermal_generation.py

echo "thermal.sif" > ELMERSOLVER_STARTINFO
ElmerSolver

