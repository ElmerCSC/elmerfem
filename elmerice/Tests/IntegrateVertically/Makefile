# Test for IntegrateVertically Solver
#
run:	
	$(ELMER_GRID) 1 2 mesh.grd
	$(ELMER_SOLVER) integratevertically.sif

clean:
	/bin/rm *.log  *.dat *.dat.names
	/bin/rm IntegrateVertically
	/bin/rm IntegrateVertically.*
	/bin/rm -r mesh
