Test Case for 3D Calving
====

This is a simple test case for the 3D calving model.


What to do if this test fails
====

It's very likely that this test is failing due to gmsh dependency. Build gmsh without MPI support and then put the following in the simulation section of your sif:

Gmsh Path = String "/your/new/gmsh/installation/gmsh"

