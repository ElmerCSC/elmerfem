INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( case.sif case.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO angle.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

# The single-task run is kept so that the parallel result is compared against a
# serial one computed with the very same solver settings, not just against the
# stored reference norm.
ADD_ELMER_TEST(CoupledPoissonAitken NPROCS 1 2)
