How to run the test:
--------------------
ElmerGrid 2 2 Syn1Mesh -partition 4 1 1
mpirun -np 4 ElmerSolver_mpi SaveGridDataNetCDFTest.sif 


How to check the results:
-------------------------
ncdiff -O  NetCDFTest_ref.nc Syn1Mesh/NetCDFTest.nc -o diff.nc && ncap2 -A  -s "maxzs=zs.max();print(maxzs)" diff.nc

Notes:

1. ncdiff and ncap2 are NCO tool

2. The max difference between the fields (ref VS created) , should not exceed 10e-6
