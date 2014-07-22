# Test for Porous Solvers
#
run:	
	$(GMSH) ./PROG/mesh_504.geo -3 -o mesh_504.msh
	$(ELMER_GRID) 14 2 mesh_504.msh -autoclean
	$(ELMER_SOLVER) cube504_g.sif

clean:
	/bin/rm mesh_504.msh
	/bin/rm Porous.*
	/bin/rm -r mesh_504
	/bin/rm *.dat *.dat.names *.log 

