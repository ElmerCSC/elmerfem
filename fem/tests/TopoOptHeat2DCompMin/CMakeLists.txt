IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(case.sif case.sif COPYONLY)

  file(COPY tree.grd ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_ELMER_TEST(TopoOptHeat2DCompMin LABELS slow)
  # This test's assembly is not thread safe; pin it to a single OpenMP thread
  # so it runs deterministically (runtest.cmake honors a pre-set OMP_NUM_THREADS).
  SET_TESTS_PROPERTIES(TopoOptHeat2DCompMin PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
ENDIF()
