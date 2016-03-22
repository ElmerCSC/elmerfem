# Computes Tresca-stress for a simple 2d geometry
#
run:
	$(F90) -c Tresca.f90
	$(LD) -o Tresca$(SHL_EXT) Tresca$(OBJ_EXT) $(LIBS)
	$(ELMER_SOLVER)

clean:
	/bin/rm test.log temp.log mon.out Tresca$(SHL_EXT) Tresca$(OBJ_EXT) so_locations
