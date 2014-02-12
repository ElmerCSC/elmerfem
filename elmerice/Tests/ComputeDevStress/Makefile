# Test with ComputeDevStress Solver
#
run:
	$(ELMER_GRID) 1 2 rectangle.grd
	$(ELMER_SOLVER) computeDevStress.sif

clean:
	/bin/rm *.dat *.log *.dat.names
	/bin/rm ComputeDevStress.*
	/bin/rm -r rectangle
