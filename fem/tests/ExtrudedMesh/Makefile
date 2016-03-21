# Extruded mesh test
#
run:
	$(F90) -c Poisson.f90
	$(LD) -o Poisson$(SHL_EXT) Poisson$(OBJ_EXT) $(LIBS)
	$(ELMER_GRID) 1 2 square
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log temp.log mon.out Poisson$(SHL_EXT) Poisson$(OBJ_EXT) so_locations
	/bin/rm -r square
