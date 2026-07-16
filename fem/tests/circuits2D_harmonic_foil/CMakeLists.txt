IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(sif/1383-components-foilwinding.sif sif/1383-components-foilwinding.sif COPYONLY) 

  file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  file(COPY sif/1383-circuits-components.definitions DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  file(COPY 1383/mesh.boundary 1383/mesh.header 1383/mesh.nodes 1383/mesh.elements 1383/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/1383")

  ADD_ELMER_TEST(circuits2D_harmonic_foil LABELS 2D circuits harmonic mgdyn)
ENDIF()
