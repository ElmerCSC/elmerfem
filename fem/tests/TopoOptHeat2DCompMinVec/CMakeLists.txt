IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(case.sif case.sif COPYONLY)

  file(COPY tree.grd ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  # Same as TopoOptHeat2DCompMin but with vector assembly, which takes the
  # local matrix copy path of the vectorized HeatSolve.
  ADD_ELMER_TEST(TopoOptHeat2DCompMinVec LABELS slow)
  SET_TESTS_PROPERTIES(TopoOptHeat2DCompMinVec PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
ENDIF()
