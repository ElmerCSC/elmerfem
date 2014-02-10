# Test with DEnsity Solver
#
run:
	$(ELMER_GRID) 1 2 mesh.grd
	$(ELMER_SOLVER) density.sif

clean:
	/bin/rm *.dat *.dat.names *.log 
	/bin/rm Density.*
	/bin/rm -r mesh
