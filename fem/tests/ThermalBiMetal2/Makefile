# Test for bimetal bending under temperature lift
#
run:
	$(ELMER_GRID) 1 2 beams
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log temp.log mon.out
	/bin/rm -r beams
