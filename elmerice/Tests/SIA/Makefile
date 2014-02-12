# Test for GroundedSolver and GetHydrostaticLoads Solvers
#
run:	
	$(FORTnosh) ./PROG/MshGlacierSynthetic.f90 ./PROG/fbed.f90 ./PROG/fsurf.f90 -o ./MshGlacierSynthetic
	$(ELMER_GRID) 1 2 mesh_A.grd
	./MshGlacierSynthetic
	$(ELMER_SOLVER) sia.sif

clean:
	/bin/rm *.dat *.dat.names *.log 
	/bin/rm  MshGlacierSynthetic
	/bin/rm SIA.*
	/bin/rm -r mesh_A
