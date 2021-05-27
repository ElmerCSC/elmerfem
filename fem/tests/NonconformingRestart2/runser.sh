elmerf90 -o InterpSolver.so InterpSolver.F90

ElmerGrid 1 2 mesh_coarse.grd 
ElmerSolver coarse.sif

ElmerGrid 1 2 mesh_fine.grd 
ElmerSolver_mpi fine.sif
