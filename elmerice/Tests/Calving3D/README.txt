Test Case for 3D Calving
====

More detailed documentation on Elmer's 2D Calving can be found at:

http://elmerice.elmerfem.org/wiki/doku.php?id=problems:calving


What to do if this test fails
====

It's very likely that this test is failing due to gmsh dependency. Build gmsh without MPI support and then put the following in the simulation section of your sif:

Gmsh Path = String "/your/new/gmsh/installation/gmsh"

