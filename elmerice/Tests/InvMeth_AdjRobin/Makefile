# Test with Adjoint_Robin Solver
# Robin inverse method; optimisation of the slip coef.
run:	
	$(FORT) ./PROG/GradientValidation.f90 -o ./GradientValidation
	$(FORT) ./PROG/USF_Init.f90 -o ./USF_Init
	$(ELMER_GRID) 1 2 mesh2D -metis 4 4
	mpirun -n 4 ElmerSolver_mpi
	chmod u+x script.sh
	./script.sh

clean:
	/bin/rm gradientnormadjoint_adj_robin.dat Cost_Adj_Robin.dat
	/bin/rm *.dat *.dat.names
	/bin/rm GradientValidation USF_Init
	/bin/rm InvMeth_AdjRobin.txt
	/bin/rm  *.log *.out
	/bin/rm -r mesh2D
