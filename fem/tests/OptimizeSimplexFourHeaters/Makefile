# Optimized heater 
#
run:
	$(F90) -c CostFunction.f90
	$(LD) -o CostFunction$(SHL_EXT) CostFunction$(OBJ_EXT) $(LIBS)
	$(ELMER_GRID) 1 2 spots
	$(ELMER_SOLVER)

clean:
	/bin/rm CostFunction$(OBJ_EXT) CostFunction$(SHL_EXT) test.log temp.log mon.out so_locations optimize.dat optimize-best.dat
	/bin/rm -r spots
