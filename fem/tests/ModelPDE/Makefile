# A model PDE solver for the convection-diffusion equation
#
run:
	$(F90) -c AdvDiff.f90
	$(LD) -o AdvDiff$(SHL_EXT) AdvDiff$(OBJ_EXT) $(LIBS)
#	$(ELMER_GRID) 1 2 square
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log temp.log mon.out AdvDiff$(SHL_EXT) AdvDiff$(OBJ_EXT) so_locations
#	/bin/rm -r square
