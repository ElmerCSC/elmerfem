# Linear elasticity solver test in 2D / adaptivity
#
run:
	$(ELMER_GRID) 1 2 mesh -triangles
	$(ELMER_SOLVER)

clean:
	/bin/rm test.log temp.log mon.out
	/bin/rm -r mesh RefinedMesh*
