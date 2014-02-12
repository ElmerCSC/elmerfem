# Test for Grounded Solvers
#
run:	
	$(FORT) ./PROG/bedrock.f90 ./PROG/fbed.f90 -o bedrock
	$(ELMER_GRID) 1 2 Cube.grd
	$(ELMER_SOLVER) grounded.sif

clean:
	/bin/rm *.dat *.dat.names *.log 
	/bin/rm bedrock
	/bin/rm Grounded.*
	/bin/rm -r Cube
