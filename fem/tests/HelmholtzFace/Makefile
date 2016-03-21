# Solves vector Helmholtz equation using face elements
#
run:
	$(F90) -c Helmholtz.f90
	$(LD) -o Helmholtz$(SHL_EXT) Helmholtz$(OBJ_EXT) $(LIBS)
	$(ELMER_GRID) 1 2 square.grd -triangles
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log temp.log mon.out Helmholtz$(SHL_EXT) Helmholtz$(OBJ_EXT) so_locations  ratkaisu.ep
	/bin/rm -r square
