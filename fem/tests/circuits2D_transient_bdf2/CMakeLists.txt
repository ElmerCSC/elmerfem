IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(sif/bdf2.sif sif/bdf2.sif COPYONLY)

  file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  file(COPY sif/bdf2-circuit.definitions DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  file(COPY 1381/mesh.boundary 1381/mesh.header 1381/mesh.nodes 1381/mesh.elements 1381/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/1381")

  ADD_ELMER_TEST(circuits2D_transient_bdf2 LABELS 2D circuits transient mgdyn)
ENDIF()
