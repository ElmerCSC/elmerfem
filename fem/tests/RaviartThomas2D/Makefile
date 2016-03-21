# A test for the mixed approximation of a 2-D Poisson problem with the lowest-order Raviart-Thomas elements
#
run:
	$(F90) -c RTSolve.f90
	$(LD) -o RTSolve$(SHL_EXT) RTSolve$(OBJ_EXT) $(LIBS)
	$(ELMER_SOLVER)

clean:
	/bin/rm test.log temp.log mon.out RTSolve$(SHL_EXT) RTSolve$(OBJ_EXT) so_locations
