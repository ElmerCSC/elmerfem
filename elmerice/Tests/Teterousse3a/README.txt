# Test with 2 FreeSurface Solver (Bottom and Top)
#
cp $ELMER_HOME/share/elmersolver/lib/FreeSurfaceSolver.so FreeSurface1
elmerf90 ./PROG/USF_TR.f90 -o USF_TR
ElmerGrid 14 2 teterousse.msh -autoclean -order 1.0 0.1 0.01
ElmerSolver teterousse3a.sif

Results established:
------------------
19.03.2015
Laure Tavard,LGGE
Froggy cluster (CIMENT: Grenoble University HPC centre)
Revision 58f71b4    
