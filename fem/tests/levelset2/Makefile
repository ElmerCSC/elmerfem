# Filling of a volume by a denser fluid using level set approach
#
run:
	$(F90) -c MyStuff.f90
	$(LD) -o MyStuff$(SHL_EXT) MyStuff$(OBJ_EXT) $(LIBS)
	$(ELMER_GRID) 1 2 fill
	$(ELMER_SOLVER)

clean:
	/bin/rm MyStuff$(OBJ_EXT) MyStuff$(SHL_EXT) test.log temp.log mon.out so_locations
	/bin/rm -r fill
