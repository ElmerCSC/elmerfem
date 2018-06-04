# FreeCADBatchFEMTools - A library for using FreeCAD for FEM preprocessing in batch mode
FreeCADBatchFEMTools.py is a set of tools for handling batch type preprocessing for 
Elmer Solver. The library was written by Trafotek Oy in collaboration with VTT Technical 
Reseach Centre of Finland Ltd supported by Business Finland Oy.

# Features
- Inteded for batch mode usage (automatic preprocessing)
- Creating mesh objects for FreeCAD
- Manipulating the mesh sizes of solids
- Creating meshes using Gmsh
- Running ElmerGrid and exporting mesh for ElmerSolver
- Selecting and naming bodies and boundaries for ElmerSolver

# Requirements
- FreeCAD 0.18 (tested with Libs: 0.18R13824 (Git))
- Gmsh 3
- ElmerGrid (Elmer 8.3)

# Tests

Tests are located at ./tests folder.

# Notes
- To run scipts in batch mode that use FreeCADBatchFEMTools, use:
$ FreeCAD -c $PWD/script_name.py

# Authors
- Eelis Takala, Trafotek Oy
- Janne Keranen, VTT Technical Reseach Centre of Finland Ltd 

