# Test for ElmerIceFunction Buoyancy
#
run:	
	$(ELMER_GRID) 1 2 mesh.grd
	$(ELMER_SOLVER) buoyancy.sif

clean:
	/bin/rm *.dat *.log *.dat.names
	/bin/rm Buoyancy.*
	/bin/rm -r mesh
