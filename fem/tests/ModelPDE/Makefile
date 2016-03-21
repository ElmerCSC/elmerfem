# A model PDE solver for the convection-diffusion equation
#
run:
	$(F90) -c ModelPDE.f90
	$(LD) -o ModelPDE$(SHL_EXT) ModelPDE$(OBJ_EXT) $(LIBS)
#	$(ELMER_GRID) 1 2 square
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log temp.log mon.out ModelPDE$(SHL_EXT) ModelPDE$(OBJ_EXT) so_locations
#	/bin/rm -r square
