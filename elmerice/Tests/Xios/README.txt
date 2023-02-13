# Test XIOS
#  => Generate a netcdf output file using XIOS then check that we can read the file and that the connectivity is correct 
#  => XIOS configiurations files: iodef.xml and context_elmerice.xml
# 
To run the test:
1. Compile Check.F90 (requires NETCDF)
2. Generate the mesh for N partitions:
	ElmerGrid 1 2 rectangle -metis N
3. Run the test:
	- XIOS attached mode:
		mpirun -np N ElmerSolver_mpi Case.sif
	- XIOS detached mode
		mpirun -np N ElmerSolver_mpi Case.sif : -np N2 xios_server

Results established:
------------------
04.05.2022
Fabien Gillet-Chaulet, IGE
