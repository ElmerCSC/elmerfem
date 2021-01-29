# src

Source codes and input files

## Content

* **Configuration files**:
	- *PARAMETERS.sif* : text file with MATC parameters to run the tests

* **Elmer Codes**:
	- *RAMP.F90* : A collection of user functions to compute the analytical solutions for the ramp ice shelf
	- *Random.F90* : A user function to initialise random variables
	- *Makefile* : Makefile for the compilation

* **Meshing utilities**:
	- *rectangle.geo*: template gmsh .geo file to mesh a rectangle
	- *MakeMesh.sh*: bash script to create the mesh and convert it to Elmer

* **Others**:
	- *MakeObs.py* : Python3 script to generate synthetic thickness observations


