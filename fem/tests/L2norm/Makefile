run:
	$(F90) -c Difference.f90
	$(LD) -o Difference$(SHL_EXT) Difference$(OBJ_EXT) $(LIBS)
	$(ELMER_SOLVER)

clean:
	/bin/rm test.log temp.log mon.out Difference$(SHL_EXT) Difference$(OBJ_EXT) so_locations
