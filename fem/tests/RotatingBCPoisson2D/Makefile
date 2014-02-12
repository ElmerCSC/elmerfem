# A test for simple 2D rotating BCs 
run:
	$(ELMER_GRID) 14 2 mortar.msh -autoclean
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log temp.log mon.out
	/bin/rm -r mortar
