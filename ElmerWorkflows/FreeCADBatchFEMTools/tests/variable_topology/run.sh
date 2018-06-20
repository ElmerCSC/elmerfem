#!/bin/bash

rm -r variable_topology variable_topology.sif variable_topology.unv scalars.dat* ELMERSOLVER_STARTINFO 2> /dev/null

FreeCAD -c $PWD/variable_topology_generation.py

sleep 1

echo "variable_topology.sif" > ELMERSOLVER_STARTINFO
ElmerSolver

