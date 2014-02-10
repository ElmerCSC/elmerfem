# Test for IntegratedVelocity Solvers
#
run:	
	$(ELMER_GRID) 1 2 rectangle.grd
	$(ELMER_SOLVER) test.sif

clean:
	/bin/rm *.dat *.dat.names *.log 
	/bin/rm IntegratedVelocity.*
	/bin/rm -r rectangle
