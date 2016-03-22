# Water pressure in land dyke computed from steady-state Richards equation
# with van Genuchten material laws

run:
	$(ELMER_GRID) 1 2 land_gauss
	$(ELMER_SOLVER)

clean:
	/bin/rm test.log temp.log mon.out
	/bin/rm -r land_gauss
