IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(waveguide.sif waveguide.sif COPYONLY)

  file(COPY block.sif ELMERSOLVER_STARTINFO shoebox_tets DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_ELMER_TEST(VectorHelmholtzWaveguideQuadBlock LABELS whitney)
ENDIF()
