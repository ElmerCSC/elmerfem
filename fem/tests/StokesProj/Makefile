# Navier-Velocity solver with pressure projection
#
run:
	$(F90) -c Velocity.f90 Pressure.f90
	$(LD) -o Velocity$(SHL_EXT) Velocity$(OBJ_EXT) $(LIBS)
	$(LD) -o Pressure$(SHL_EXT) Pressure$(OBJ_EXT) $(LIBS)
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log temp.log mon.out Velocity$(SHL_EXT) Velocity$(OBJ_EXT) Pressure$(SHL_EXT) Pressure$(OBJ_EXT) so_locations
