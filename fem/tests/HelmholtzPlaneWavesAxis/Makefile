# Simple test for plane waves in 2D and axisymmetric case
#

run:
	$(ELMER_GRID) 1 2 rect.grd
	$(ELMER_SOLVER)

clean:
	/bin/rm mon.out
	/bin/rm -rf temp.log test.log 
	/bin/rm -r rect