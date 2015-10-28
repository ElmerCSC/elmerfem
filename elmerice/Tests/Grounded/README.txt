To Execute the test:
--------------------
Compile
elmerf90 ./PROG/bedrock.f90 ./PROG/fbed.f90 -o bedrock

Make the mesh: 
ElmerGrid 1 2 Cube.grd

Execute
ElmerSolver grounded.sif

Nodes in contact at the bed should evolve from 0 to all the nodes at the base

Results established:
------------------
19.03.2015
Laure Tavard,LGGE
Froggy cluster (CIMENT: Grenoble University HPC centre)
Revision 58f71b4
