# Flow-structure interaction (FSI) test in 2D, bending beam
#
run:
	$(F90) -c FsiStuff.f90
	$(LD) -o FsiStuff$(SHL_EXT) FsiStuff$(OBJ_EXT) $(LIBS)
	$(ELMER_GRID) 1 2 fsi
	$(ELMER_SOLVER)

clean:
	/bin/rm FsiStuff$(OBJ_EXT) FsiStuff$(SHL_EXT) test.log temp.log mon.out so_locations
	/bin/rm -r fsi
