# Test with Adjoint_Beta Solver
# Control inverse method; optimisation of the slip coef.
run:	
	$(FORT) ./PROG/GradientValidation.f90 -o GradientValidation
	$(FORT) ./PROG/USF_Init.f90 -o USF_Init
	$(ELMER_GRID) 1 2 mesh2D -metis 4 4
	mpirun -n 4 ElmerSolver_mpi
	chmod u+x script.sh
	./script.sh

clean:
	/bin/rm gradientnormadjoint_adj_beta.dat Cost_Adj_Beta.dat
	/bin/rm *.dat *.dat.names
	/bin/rm GradientValidation USF_Init
	/bin/rm InvMeth_AdjBeta.txt
	/bin/rm *.log *.out
	/bin/rm -r mesh2D
