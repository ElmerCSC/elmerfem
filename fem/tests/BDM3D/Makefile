# A test for the mixed approximation of a 3-D Poisson problem with the BDM elements of degree k=1
#
run:
	$(F90) -c BDMSolve.f90
	$(LD) -o BDMSolve$(SHL_EXT) BDMSolve$(OBJ_EXT) $(LIBS)
	$(ELMER_SOLVER)

clean:
	/bin/rm test.log temp.log mon.out BDMSolve$(SHL_EXT) BDMSolve$(OBJ_EXT) so_locations
