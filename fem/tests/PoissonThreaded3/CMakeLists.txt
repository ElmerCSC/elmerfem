INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

IF(NOT(WITH_OpenMP)) # This test is currently broken when OMP is enabled
CONFIGURE_FILE( ttest.sif ttest.sif COPYONLY)
ADD_ELMERTEST_MODULE(PoissonThreaded3 Poisson Poisson.F90)

file(COPY Poisson.F90 cube.grd ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(PoissonThreaded3 LABELS threaded)
ELSE()
  MESSAGE(WARNING "PoissonThreaded3 test is currently disabled if OpenMP is enabled")
ENDIF(NOT(WITH_OpenMP))
