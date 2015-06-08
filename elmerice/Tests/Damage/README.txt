How to run the test:
--------------------
Require two calls of the FreeSurfaceSolver. 
It is then needed to make a copy of the object file FreeSurfaceSolver.so (or FreeSurfaceSolver.dylib for mac)

cp $ELMER_HOME/share/elmersolver/lib/FreeSurfaceSolver.so FreeSurfaceSolver1
ElmerGrid 1 2 mesh.grd
ElmerSolver damage.sif

Results established:
---------------------
19.03.2015
Laure Tavard,LGGE
Froggy cluster (CIMENT: Grenoble University HPC centre)
Revision 58f71b4

