# Test for Thickness Solver(3D)
#
run:	
	$(ELMER_GRID) 1 2 cube.grd
	$(ELMER_SOLVER) NS_3D.sif

clean:
	/bin/rm *.dat *.dat.names *.log 
	/bin/rm ThicknessSolver.*
	/bin/rm -r cube
