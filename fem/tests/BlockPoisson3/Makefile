# A poisson equation solver using the coupled system routines
#
run:
	$(F90) -c BlockPoissonAssembly.f90 
	$(LD) -o BlockPoissonAssembly$(SHL_EXT) BlockPoissonAssembly$(OBJ_EXT) $(LIBS)
	$(ELMER_GRID) 1 2 square
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log temp.log mon.out BlockPoissonAssembly$(SHL_EXT) BlockPoissonAssembly$(OBJ_EXT) so_locations
	/bin/rm -r square
