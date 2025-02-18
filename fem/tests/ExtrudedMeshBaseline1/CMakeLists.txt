IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE( case_baseline.sif case_baseline.sif COPYONLY)
  ADD_ELMERTEST_MODULE(ExtrudedMeshBaseline1 Poisson Poisson.f90)

  file(COPY Poisson.f90 square.grd ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_ELMER_TEST(ExtrudedMeshBaseline1 LABELS quick useextrude)
ENDIF()
