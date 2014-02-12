# Test with FlowDepth Solver
#
run:	
	$(FORTnosh) ./PROG/MshGlacierSynthetic.f90 ./PROG/fbed.f90 ./PROG/fsurf.f90 -o MshGlacierSynthetic
	$(ELMER_GRID) 1 2 mesh.grd
	./MshGlacierSynthetic
	$(ELMER_SOLVER) flowdepth.sif

clean:
	/bin/rm *.dat *.dat.names *.log 
	/bin/rm MshGlacierSynthetic
	/bin/rm FlowDepth.*
	/bin/rm -r mesh
