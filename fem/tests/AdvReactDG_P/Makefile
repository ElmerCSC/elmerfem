# advection-reaction equation by discontinuous galerkin & p element basis
#
run:
	$(F90) -c AdvReactDG_P.F90
	$(LD) -o AdvReactDG_P$(SHL_EXT) AdvReactDG_P$(OBJ_EXT) $(LIBS)
	$(ELMER_GRID) 1 2 square
	$(ELMER_SOLVER)

clean:
	/bin/rm test.log temp.log AdvReactDG_P$(OBJ_EXT) AdvReactDG_P$(SHL_EXT)
	/bin/rm -r square

