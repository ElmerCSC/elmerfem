#!/bin/bash

# A test case for FreeCAD automatic scripting with thermo-structural (Elmer test is a modified fem/tests/ThermalBiMetal2)
# original date: November 2019
# Author: Eelis Takala
# email: eelis.takala@gmail.com
FreeCAD -c ${PWD}/cylinders.py
ElmerGrid 8 2 cylinders.unv -autoclean
ElmerSolver
#rm cylinders TEST.PASSED -r
