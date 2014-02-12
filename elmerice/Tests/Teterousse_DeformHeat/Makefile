# Test for teterousse modified with ComputeStrinRate Solver
#
run:	
	$(FORT) ./PROG/USF_TR.f90 -o USF_TR
	$(ELMER_GRID) 14 2 teterousse.msh -autoclean -order 1.0 0.1 0.01
	$(ELMER_SOLVER) test.sif

clean:
	/bin/rm *.dat *.dat.names *.log 
	/bin/rm  USF_TR
	/bin/rm Teterousse_DeformHeat.*
	/bin/rm -r teterousse


