run:
	$(F90) -c TangentForce.f90
	$(LD) -o TangentForce$(SHL_EXT) TangentForce$(OBJ_EXT) $(LIBS)
	$(ELMER_GRID) 1 2 square.grd
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log temp.log mon.out TangentForce$(SHL_EXT) TangentForce$(OBJ_EXT) so_locations
	/bin/rm -r square
