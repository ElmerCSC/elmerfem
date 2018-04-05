
To execute this test:

1/ Make the mesh on the frontal interactively:
ElmerGrid 1 2 mesh_B.grd

(should create the directory mesh_B/ and mesh.* files inside)

2/ Run the FE simulations :             
! For the Glen's law in the distrib
ElmerSolver test_glenDistrib.sif
! For the MATC definition of the Glen's law, see elmerice/examples
ElmerSolver test_glenMATC.sif 

Both should give the same results. Note the activation of the Newton linearization 
which allows a very fast convergence of the non-linear iteration down to 1.0e-8

Results established:
------------------
19.03.2015
Laure Tavard,LGGE
Froggy cluster (CIMENT: Grenoble University HPC centre)
Revision 58f71b4
------------------
05.04.2018
Updated F. Gillet-Chaulet. 
- Change ref. norm and tol. for stokes
- test result from flowdepth


