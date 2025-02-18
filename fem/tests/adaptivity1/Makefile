# Adaptive solver test with remeshing
#

run:
	/bin/mkdir -p Ldomain
	$(ELMER_MESH2D) Ldomain.mif .
	/bin/mv mesh.* Ldomain
	/bin/cp Ldomain.mif Ldomain
	$(ELMER_SOLVER)

clean:
	/bin/rm test.log temp.log bgmesh fort.12 mon.out
	/bin/rm -rf Ldomain
	/bin/rm -rf RefinedMesh*
