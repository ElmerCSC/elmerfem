# A poisson equation solver using p elements
#
run:
	$(F90) -c WaveEq.f90
	$(LD) -o WaveEq$(SHL_EXT) WaveEq$(OBJ_EXT) $(LIBS)
	$(ELMER_GRID) 1 2 square
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log temp.log mon.out WaveEq$(SHL_EXT) WaveEq$(OBJ_EXT) so_locations
	/bin/rm -r square
