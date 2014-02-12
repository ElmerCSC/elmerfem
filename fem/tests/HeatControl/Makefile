# Test heater control for a simple square with given target temperature
#
run:
	$(ELMER_GRID) 1 2 square
	$(ELMER_SOLVER)

clean:
	/bin/rm test.log temp.log mon.out
	/bin/rm -r square
