# An interpolation test using H(curl)-conforming finite elements
#
run:
	$(F90) -c EdgeElementSolve.f90
	$(LD) -o EdgeElementSolve$(SHL_EXT) EdgeElementSolve$(OBJ_EXT) $(LIBS)
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log temp.log mon.out EdgeElementSolve$(SHL_EXT) EdgeElementSolve$(OBJ_EXT) so_locations

