IF(HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE( case.sif case.sif COPYONLY)

  # mesh (from the static test mgdyn_bh) and B-H curve
  file(COPY ELMERSOLVER_STARTINFO HB mesh_tet DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_ELMER_TEST(mgdyn_harmonic_bh_newton LABELS whitney harmonic slow)
ENDIF()
