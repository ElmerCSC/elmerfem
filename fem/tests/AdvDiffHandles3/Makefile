# Transient diffusion test
#
run:
	$(F90) -c AdvDiffHandles3.f90 
	$(LD) -o AdvDiffHandles3$(SHL_EXT) AdvDiffHandles3$(OBJ_EXT) $(LIBS)
	$(ELMER_GRID) 1 2 Step
	$(ELMER_SOLVER)

clean:
	/bin/rm temp.log test.log mon.out
	/bin/rm -r Step
	/bin/rm AdvDiffHandles3$(SHL_EXT) AdvDiffHandles3$(OBJ_EXT)