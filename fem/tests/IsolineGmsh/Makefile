# Elastic displacement limited from the center
#
run:
	$(ELMER_GRID) 1 2 angle
	$(ELMER_SOLVER)

clean:
	/bin/rm test.log temp.log mon.out
	/bin/rm -r isoline angle
	/bin/rm mesh.geo
