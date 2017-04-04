INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( case.sif case.sif COPYONLY)
ADD_ELMERTEST_MODULE(CoupledPoisson9 PoissonAssembly PoissonAssembly.f90)

file(COPY ELMERSOLVER_STARTINFO squares.grd PoissonAssembly.f90 DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(CoupledPoisson9)
ADD_ELMER_LABEL(CoupledPoisson9 quick)
