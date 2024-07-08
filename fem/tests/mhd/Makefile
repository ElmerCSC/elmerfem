# magnetohydrodynamics with two bodies
#
run:
	$(F90) -c  exmag.f90
	$(LD) -o exmag$(SHL_EXT) exmag$(OBJ_EXT) $(LIBS)
	$(ELMER_GRID) 1 2 mesh.grd
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log temp.log mon.out exmag$(OBJ_EXT) exmag$(SHL_EXT)
	/bin/rm -r mesh
