How to run the test:
--------------------
ElmerGrid 1 2 cube.grd -partition 4 1 1 -periodic 1 1 0 
mpirun -np 4 ElmerSolver_mpi test.sif

The solution in 'Stress' must be indetical to 'StressAna' on the bottom boundary of the mesh.  

