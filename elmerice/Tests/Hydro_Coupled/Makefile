# Test for Hydrology Solvers
#

run:	
	$(ELMER_GRID) 1 2 cube.grd -boundbound 1 5 7 -boundbound 2 5 7 -boundbound 3 5 7 -boundbound 4 5 8
	$(ELMER_SOLVER) Coupled.sif

clean:
	/bin/rm *.dat *.log *.dat.names
	/bin/rm Hydro_Coupled.*
	/bin/rm -r cube
