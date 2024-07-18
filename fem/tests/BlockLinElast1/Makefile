# A poisson equation solver using the coupled system routines
#
run:
	$(F90) -c BlockLinElast.f90 
	$(LD) -o BlockLinElast$(SHL_EXT) BlockLinElast$(OBJ_EXT) $(LIBS)
	$(ELMER_GRID) 1 2 square
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log temp.log mon.out BlockLinElast$(SHL_EXT) BlockLinElast$(OBJ_EXT) so_locations
	/bin/rm -r square
