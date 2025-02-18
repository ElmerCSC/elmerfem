# A generic turing system
#
run:
	$(F90) -c Turing.f90
	$(LD) -o Turing$(SHL_EXT) Turing.o $(LIBS)
	$(ELMER_GRID) 1 2 square
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log temp.log mon.out Turing$(SHL_EXT) Turing.o so_locations
	/bin/rm -r square
