# Test for DG Solvers
#
run:	
	$(FORT) ./PROG/DGtoNodalVariable.f90 -o DGtoNodalVariable
	$(FORT) ./PROG/InitializeDGVariable.f90 -o InitializeDGVariable
	$(ELMER_GRID) 1 2 cube.grd
	$(ELMER_SOLVER) density.sif

clean:
	/bin/rm *.dat *.dat.names *.log 
	/bin/rm InitializeDGVariable DGtoNodalVariable
	/bin/rm DGsolver.*
	/bin/rm -r cube
