# Test for natural convection. The heating is set on vertical sides 
# to speed up the test case.
#
run:
	$(ELMER_GRID) 1 2 square
	$(ELMER_SOLVER)

clean:
	/bin/rm test.log temp.log mon.out
	/bin/rm -r square
