# A simple heat equation test
#
run:
	$(ELMER_GRID) 1 2 global
	$(ELMER_GRID) 1 2 local
	$(ELMER_GRID) 1 2 reference
	$(ELMER_SOLVER)


clean:
	/bin/rm test.log temp.log mon.out 
	/bin/rm -r global local reference
