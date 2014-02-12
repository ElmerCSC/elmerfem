# A (navier-)stokes solver with p-bubbles and/or p2/p1 elements
#
# feasible elements (these are just examples):
# -------------------------------------------
# 303b1    - tri with one bubble
# 303e1b1  - p2/p1 triangle
#
# 404b4    - quad with four bubbles
# 404e1b1  - q2/q1 quad
#
#  etc.
#
run:
	$(F90) -c Stokes.f90
	$(LD) -o Stokes1$(SHL_EXT) Stokes$(OBJ_EXT) $(LIBS)
	$(ELMER_SOLVER)

clean:
	/bin/rm test.log temp.log mon.out Stokes1$(SHL_EXT) Stokes$(OBJ_EXT) so_locations
