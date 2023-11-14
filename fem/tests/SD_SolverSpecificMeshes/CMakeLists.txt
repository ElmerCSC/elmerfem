INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( ptest.sif ptest.sif COPYONLY)
ADD_ELMERTEST_MODULE(SD_SolverSpecificMeshes Poisson Poisson.f90)

file(COPY Poisson.f90 square.grd ELMERSOLVER_STARTINFO square square2 DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(SD_SolverSpecificMeshes)
ADD_ELMER_LABEL(SD_SolverSpecificMeshes quick serendipity)
