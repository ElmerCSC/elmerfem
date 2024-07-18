# test 1st order time integration
#
run:
	$(F90) -c TimeIntegrate.f90 Source.f90
	$(LD) -o TimeIntegrate$(SHL_EXT) TimeIntegrate$(OBJ_EXT) $(LIBS)
	$(LD) -o Source$(SHL_EXT) Source$(OBJ_EXT) $(LIBS)
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log temp.log mon.out TimeIntegrate$(SHL_EXT) TimeIntegrate$(OBJ_EXT)
	/bin/rm Source$(SHL_EXT) Source$(OBJ_EXT) so_locations 
