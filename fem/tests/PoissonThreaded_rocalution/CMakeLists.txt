IF(HAVE_ROCALUTION)

  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE( ttest.sif ttest.sif COPYONLY)
  ADD_ELMERTEST_MODULE(PoissonThreaded_rocalution Poisson Poisson.F90)

  file(COPY Poisson.F90 cube.grd ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_ELMER_TEST(PoissonThreaded_rocalution LABELS threaded rocalution)

ENDIF(HAVE_ROCALUTION)
