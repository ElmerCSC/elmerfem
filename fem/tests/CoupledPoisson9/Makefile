# A poisson equation solver using the coupled system routines
#
run:
	$(F90) -c PoissonAssembly.f90 
	$(LD) -o PoissonAssembly$(SHL_EXT) PoissonAssembly$(OBJ_EXT) $(LIBS)
	$(ELMER_GRID) 1 2 squares
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log temp.log mon.out PoissonAssembly$(SHL_EXT) PoissonAssembly$(OBJ_EXT) so_locations
	/bin/rm -r squares
