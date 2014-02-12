# Electrostatics and adaptivity with remeshing
#

run:
	/bin/mkdir -p kulma
	$(ELMER_MESH2D) kulma.mif kulma
	/bin/cp kulma.mif kulma
	$(ELMER_SOLVER)

clean:
	/bin/rm test.log temp.log bgmesh fort.12 mon.out
	/bin/rm -rf kulma
	/bin/rm -rf RefinedMesh*
