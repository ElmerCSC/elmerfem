run:
	ElmerGrid 14 2 blunt.msh -autoclean 
	ElmerGrid 1 2 base.grd
	ElmerGrid 2 2 blunt -in base -unite -out mesh
	ElmerGrid extrude.eg
	ElmerSolver


clean:
	/bin/rm test.log temp.log mon.out so_locations
	/bin/rm -r blunt base mesh
