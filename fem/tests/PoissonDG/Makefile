# poisson equation by discontinuous galerkin
#
run:
#	$(ELMER_GRID) 1 2 square
	$(F90) -c PoissonDG.F90
	$(LD) -o PoissonDG$(SHL_EXT) PoissonDG$(OBJ_EXT) $(LIBS)
	$(ELMER_SOLVER)

clean:
	/bin/rm test.log temp.log mon.out dg.ep PoissonDG$(SHL_EXT) PoissonDG$(OBJ_EXT)
#	/bin/rm -r square
