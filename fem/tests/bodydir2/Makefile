# A periodic boundary condition test
#
run:
	$(F90) -c Blowers.f90
	$(LD) -o Blowers$(SHL_EXT) Blowers$(OBJ_EXT) $(LIBS)
	$(ELMER_GRID) 1 2 cross
	$(ELMER_SOLVER)


clean:
	/bin/rm -r cross Blowers$(SHL_EXT) Blowers$(OBJ_EXT) so_locations
	/bin/rm test.log temp.log mon.out
