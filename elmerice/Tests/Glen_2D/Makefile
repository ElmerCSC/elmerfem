# Test for ElmerIceFunction Glen2D
#
run:	
	$(ELMER_GRID) 1 2 mesh_B.grd
	$(ELMER_SOLVER) test_glenDistrib.sif

clean:
	/bin/rm *.dat *.dat.names *.log 
	/bin/rm Glen_2D.*
	/bin/rm -r mesh_B
