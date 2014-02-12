# Test for repeated reading of existing solution data 
#
run:
	$(ELMER_GRID) 1 2 Step
	/bin/cp flow.dat Step/
	/bin/cp flow.dat.pos Step/
	$(ELMER_SOLVER)

clean:
	/bin/rm test.log temp.log mon.out
	/bin/rm -r Step
