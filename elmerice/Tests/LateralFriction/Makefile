# Test for ElmerIceFunction LateralFriction
#
run:	
	$(ELMER_GRID) 1 2 rectangle.grd
	$(ELMER_SOLVER) lateralFriction.sif

clean:
	/bin/rm *.dat *.dat.names *.log 
	/bin/rm LateralFriction.*
	/bin/rm -r rectangle
