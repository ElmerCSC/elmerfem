# Test with FlowHeight Solver
#
run:	
	$(FORTnosh) ./PROG/MshGlacierSynthetic.f90 ./PROG/fbed.f90 ./PROG/fsurf.f90 -o MshGlacierSynthetic
	$(ELMER_GRID) 1 2 mesh.grd
	./MshGlacierSynthetic
	$(ELMER_SOLVER) flowheight.sif

clean:
	/bin/rm MshGlacierSynthetic
	/bin/rm *.log *.dat *.dat.names
	/bin/rm FlowHeight.* 
	/bin/rm -r mesh
