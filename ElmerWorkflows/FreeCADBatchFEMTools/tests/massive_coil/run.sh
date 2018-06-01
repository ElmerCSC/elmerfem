#!/bin/bash

rm circuits_harmonic_massive/1962.unv circuits_harmonic_massive/TEST.PASSED circuits_harmonic_massive/ELMERSOLVER_STARTINFO 2> /dev/null
rm -r circuits_harmonic_massive/1962 2> /dev/null

FreeCAD -c $PWD/massive_generation.py

cd circuits_harmonic_massive
echo "sif/1962.sif" > ELMERSOLVER_STARTINFO
ElmerSolver
cd ..
