# A piece of artery with inlet pump and outlet 1D model
#
run:
	$(F90) -c Motion.f90 
	$(LD) -o Motion$(SHL_EXT) Motion$(OBJ_EXT) $(LIBS)
	$(ELMER_GRID) 1 2 contra.grd
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log temp.log mon.out Motion$(SHL_EXT) Motion$(OBJ_EXT) so_locations
	/bin/rm -r contra
