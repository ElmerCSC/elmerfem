# Transient diffusion test
#
run:
	$(F90) -c AdvDiffHandles.f90 
	$(LD) -o AdvDiffHandles$(SHL_EXT) AdvDiffHandles$(OBJ_EXT) $(LIBS)
	$(ELMER_GRID) 1 2 Step
	$(ELMER_SOLVER)

clean:
	/bin/rm temp.log test.log mon.out
	/bin/rm -r Step
	/bin/rm AdvDiffHandles$(SHL_EXT) AdvDiffHandles$(OBJ_EXT)