IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(sif/2077.sif sif/2077.sif COPYONLY) 

  file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  file(COPY sif/2077-circuits.definitions DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  file(COPY 2077/mesh.boundary 2077/mesh.header 2077/mesh.nodes 2077/mesh.elements 2077/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/2077")

  ADD_ELMER_TEST(circuits2D_harmonic_stranded_homogenization LABELS 2D circuits homogenization harmonic mgdyn)
ENDIF()
