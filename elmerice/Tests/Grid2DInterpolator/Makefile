# Test with Grid2DInterpolator Solver
#
run:	
	$(ELMER_GRID) 14 2 teterousse1a.msh -autoclean -order 1.0 0.1 0.01
	$(ELMER_SOLVER) teterousse1a.sif

clean:
	/bin/rm *.dat *.dat.names *.log 
	/bin/rm Grid2DInterpolator.*
	/bin/rm -r teterousse1a
