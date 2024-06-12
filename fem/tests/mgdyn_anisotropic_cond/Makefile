# A poisson equation solver using p elements
#
run:
	$(F90) -c currents.f90
	$(LD) -o currents$(SHL_EXT) currents$(OBJ_EXT) $(LIBS)
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log mon.out currents$(SHL_EXT) currents$(OBJ_EXT) so_locations toroidcurrents.mod
