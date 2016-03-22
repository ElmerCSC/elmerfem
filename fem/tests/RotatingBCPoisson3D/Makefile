# A test for simple 3D rotating BCs 
# Mesh created by extrusion
run:
	$(ELMER_GRID) 14 2 mortar.msh -autoclean
	$(ELMER_GRID) extrude.eg
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log temp.log mon.out
	/bin/rm -r mortar
	/bin/rm -r mortar3d
