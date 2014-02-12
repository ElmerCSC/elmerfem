# A test case demonstrating the use of permutated extrusion in ElmerGrid
# This time in two different phases. The 1st phase could be replaced
# with any meshing tool. 
#
run:
	$(ELMER_GRID) 1 2 squares.grd
	$(ELMER_GRID) squares2cubes.eg
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log temp.log mon.out
	/bin/rm -r squares
	/bin/rm -r cubes
