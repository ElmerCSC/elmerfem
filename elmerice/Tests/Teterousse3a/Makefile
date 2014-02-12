# Test with 2 FreeSurface Solver (Bottom and Top)
#
run:	
	$(FORT) ./PROG/USF_TR.f90 -o USF_TR 
	$(ELMER_GRID) 14 2 teterousse.msh -autoclean -order 1.0 0.1 0.01
	$(ELMER_SOLVER) teterousse3a.sif

clean:
	/bin/rm *.dat *.dat.names *.log 
	/bin/rm  USF_TR
	/bin/rm Teterousse3a.*
	/bin/rm -r teterousse
