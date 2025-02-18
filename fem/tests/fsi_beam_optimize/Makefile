# A FSI optimization case using restart
#
run:
	$(F90) -c WallProfiles.f90 
	$(LD) -o WallProfiles$(SHL_EXT) WallProfiles$(OBJ_EXT) $(LIBS)
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log temp.log mon.out WallProfiles$(SHL_EXT) WallProfiles$(OBJ_EXT) so_locations
	/bin/rm case_b.dat case_b.dat.names mesh/case_b.ep
