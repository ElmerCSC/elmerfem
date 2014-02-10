# poisson equation by discontinuous galerkin
#
run:
	$(ELMER_GRID) 1 2 square
	$(F90) -c Zalesak.f90
	$(LD) -o Zalesak$(SHL_EXT) Zalesak$(OBJ_EXT) $(LIBS)
	$(ELMER_SOLVER)

clean:
	/bin/rm test.log temp.log mon.out dg.ep Zalesak$(SHL_EXT) Zalesak$(OBJ_EXT)
	/bin/rm -r square
