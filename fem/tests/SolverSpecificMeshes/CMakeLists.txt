IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE( ptest.sif ptest.sif COPYONLY)
  ADD_ELMERTEST_MODULE(SolverSpecificMeshes Poisson Poisson.f90)

  file(COPY Poisson.f90 square.grd ELMERSOLVER_STARTINFO square square2 DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_ELMER_TEST(SolverSpecificMeshes LABELS quick)
ENDIF()
