# Test the Coulomb friction law with SSA
#
run:	
	$(ELMER_GRID) 1 2 rectangle.grd
	$(ELMER_SOLVER) ismip_SSA_2D_Coulomb.sif

clean:
	/bin/rm *.dat *.dat.names *.log 
	/bin/rm SSA_Coulomb.*
	/bin/rm -r rectangle
