# A test for simple 2D rotating BCs 
run:
	$(ELMER_GRID) 1 2 square.grd
	$(ELMER_GRID) 1 2 cube.grd
	$(ELMER_GRID) 2 2 cube -in square -unite -out mesh
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log temp.log mon.out
	/bin/rm -r square cube mesh
