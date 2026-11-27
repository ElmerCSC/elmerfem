IF(HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE( case.sif case.sif COPYONLY)

  # mesh and B-H curve of the static test mgdyn_bh
  file(COPY ELMERSOLVER_STARTINFO ${CMAKE_CURRENT_SOURCE_DIR}/../mgdyn_bh/mesh_tet ${CMAKE_CURRENT_SOURCE_DIR}/../mgdyn_bh/HB
       DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_ELMER_TEST(mgdyn_harmonic_bh_newton LABELS whitney harmonic slow)
ENDIF()
