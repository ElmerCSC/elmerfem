IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(sif/coil.sif sif/coil.sif COPYONLY)

  file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  file(COPY sif/coil-circuits.definitions DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  file(COPY coil/mesh.boundary coil/mesh.header coil/mesh.nodes coil/mesh.elements coil/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/coil")
  file(COPY res/dat/gittrack DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/res/dat")

  ADD_ELMER_TEST(circuits_transient_stranded_full_coil LABELS 3D transient circuits stranded mgdyn whitney)
ENDIF()
