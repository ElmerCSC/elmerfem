# Backward facing step / mesh grading when refining by
# splitting elements
#
run:
	$(ELMER_GRID) 1 2 Step;
	$(ELMER_SOLVER)

clean:
	/bin/rm test.log temp.log mon.out
	/bin/rm -r Step
