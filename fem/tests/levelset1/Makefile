# Solves vector Helmholtz equation using face elements
#
run:
	$(F90) -c InitSquare.f90
	$(LD) -o InitSquare$(SHL_EXT) InitSquare$(OBJ_EXT) $(LIBS)
	$(ELMER_GRID) 1 2 square.grd 
	$(ELMER_SOLVER)

clean:
	/bin/rm -r square
	/bin/rm test.log temp.log mon.out f.dat f.dat.names InitSquare$(SHL_EXT) InitSquare$(OBJ_EXT) so_locations 
