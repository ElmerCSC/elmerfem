IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(sif/case.sif sif/case.sif COPYONLY)

  file(COPY ELMERSOLVER_STARTINFO coil.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  file(COPY sif/circuits-a.definitions sif/circuits-b.definitions DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")

  ADD_ELMER_TEST(circuits2D_two_meshes LABELS 2D circuits harmonic mgdyn)
ENDIF()
