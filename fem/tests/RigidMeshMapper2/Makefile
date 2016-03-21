# Test rigid mesh mapper with a rotating rectangle
#
run:
	$(ELMER_GRID) 1 2 rectangles
	$(ELMER_SOLVER)

clean:
	/bin/rm test.log temp.log mon.out
	/bin/rm -r rectangles
