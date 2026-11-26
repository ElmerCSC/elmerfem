IF(HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE( case.sif case.sif COPYONLY)

  # B-H curve of the static test mgdyn_bh
  file(COPY ELMERSOLVER_STARTINFO strip ${CMAKE_CURRENT_SOURCE_DIR}/../mgdyn_bh/HB
       DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_ELMER_TEST(mgdyn2D_harmonic_bh_newton LABELS mgdyn harmonic quick)
ENDIF()
