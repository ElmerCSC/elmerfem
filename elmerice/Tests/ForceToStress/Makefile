# Test for ForceToStress Solvers
#
run:	
	$(ELMER_GRID) 1 2 cube.grd
	$(ELMER_SOLVER) test.sif

clean:
	/bin/rm *.dat *.dat.names *.log 
	/bin/rm ForceToStress.*
	/bin/rm -r cube
