# A test for simple 3D rotating BCs 
# Mesh created by extrusion
run:
	$(ELMER_GRID) 14 2 mortar_quarter.msh -autoclean
	$(ELMER_GRID) extrude_quarter.eg
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log temp.log mon.out
	/bin/rm -r mortar_quarter
	/bin/rm -r mortar3d
