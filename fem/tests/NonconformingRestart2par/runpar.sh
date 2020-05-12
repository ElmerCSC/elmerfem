ElmerGrid 1 2 mesh_coarse.grd -partition 2 1 1
echo 'coarse.sif' > ELMERSOLVER_STARTINFO
mpirun -np 2 ElmerSolver_mpi 

ElmerGrid 1 2 mesh_fine.grd -partition 3 1 1
echo 'fine.sif' > ELMERSOLVER_STARTINFO
mpirun -np 3 ElmerSolver_mpi 
