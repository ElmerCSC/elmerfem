# A contact problem where hard sphere meets elastic brick
#
run:
	$(F90) -c ContactBC.f90
	$(LD) -o ContactBC$(SHL_EXT) ContactBC$(OBJ_EXT) $(LIBS)
	$(ELMER_GRID) 1 2 cube.grd
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log temp.log mon.out ContactBC$(SHL_EXT) ContactBC$(OBJ_EXT) so_locations
	/bin/rm -r cube
