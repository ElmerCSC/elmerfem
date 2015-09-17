run:
	$(ELMER_GRID) 14 2 blunt.msh -autoclean 
	$(ELMER_GRID) 1 2 base.grd
	$(ELMER_GRID) 2 2 blunt -in base -unite -out mesh
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log temp.log mon.out so_locations
	/bin/rm -r blunt base mesh
