# Solves Rayleigh-Taylor instability with levelset 
#
run:
	$(ELMER_GRID) 1 2 beam.grd 
	$(ELMER_SOLVER)

clean:
	/bin/rm -r beam
	/bin/rm test.log temp.log mon.out 
