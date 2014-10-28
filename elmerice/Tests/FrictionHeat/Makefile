# Test with Friction Heat
#
run:
	$(ELMER_GRID) 1 2 cube.grd
	$(FORT) -o DummySolver.so DummySolver.f90	
	$(ELMER_SOLVER) test.sif

clean:
	/bin/rm *.dat *.dat.names *.log FrictionHeat.txt
	/bin/rm -r cube
	/bin/rm DummySolver.so