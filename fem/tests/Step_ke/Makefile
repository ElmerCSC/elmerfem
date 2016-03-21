# Backward facing step / k-epsilon model test
#
run:
	$(ELMER_GRID) 1 2 Step_ke;
	$(ELMER_SOLVER)

clean:
	/bin/rm test.log temp.log mon.out
	/bin/rm -r Step_ke
