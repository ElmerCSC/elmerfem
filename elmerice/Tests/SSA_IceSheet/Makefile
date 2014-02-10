# Grow an ice cap
#
run:	
	$(ELMER_GRID) 1 2 rectangle.grd
	$(ELMER_SOLVER) ssa_IceSheet.sif

clean:
	/bin/rm *.dat *.dat.names *.log 
	/bin/rm SSA_IceSheet.*
	/bin/rm -r rectangle

