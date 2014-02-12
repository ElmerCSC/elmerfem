# Advection-Diffusion equation with flux corrected transport (FCT)
#
run:
	$(F90) -c InitDist.f90 
	$(LD) -o InitDist$(SHL_EXT) InitDist$(OBJ_EXT) $(LIBS)
	$(ELMER_GRID) 1 2 square
	$(ELMER_SOLVER)

clean:
	/bin/rm temp.log test.log mon.out
	/bin/rm -r square
	/bin/rm InitDist$(SHL_EXT) InitDist$(OBJ_EXT)
