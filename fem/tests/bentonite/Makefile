# A model for heat and moisture behaviour of bentonite clay
#
run:
	$(F90) -c BentoniteSolver.F90
	$(LD) -o BentoniteSolver$(SHL_EXT) BentoniteSolver$(OBJ_EXT) $(LIBS)
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log temp.log mon.out BentoniteSolver$(SHL_EXT) BentoniteSolver$(OBJ_EXT) so_locations
