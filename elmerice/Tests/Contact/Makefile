# Test for GroundedSolver and GetHydrostaticLoads Solvers
#
run:	
	$(ELMER_GRID) 1 2 ice.grd
	$(ELMER_SOLVER) cavity.sif

clean:
	/bin/rm *.dat *.log *.dat.names
	/bin/rm Contact.*
	/bin/rm -r ice
