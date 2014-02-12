# Test for ElmerIceFunction ShapeFactor
#
run:	
	$(ELMER_GRID) 1 2 rectangle.grd
	$(ELMER_SOLVER) shapefactor2d.sif

clean:
	/bin/rm *.dat *.dat.names *.log 
	/bin/rm ShapeFactor.*
	/bin/rm -r rectangle
