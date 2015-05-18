# Test for DG Solvers

To run the test:
----------------
elmerf90 ./PROG/DGtoNodalVariable.f90 -o DGtoNodalVariable
elmerf90 ./PROG/InitializeDGVariable.f90 -o InitializeDGVariable
ElmerGrid 1 2 cube.grd
ElmerSolver density.sif

Results established:
------------------
19.03.2015
Laure Tavard,LGGE
Froggy cluster (CIMENT: Grenoble University HPC centre)
Revision 58f71b4

