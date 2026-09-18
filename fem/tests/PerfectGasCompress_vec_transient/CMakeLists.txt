IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE( PerfectGasCompress_vec_transient.sif PerfectGasCompress_vec_transient.sif COPYONLY)

  file(COPY PerfectGasCompress_vec_transient.grd ELMERSOLVER_STARTINFO  DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_ELMER_TEST(PerfectGasCompress_vec_transient LABELS matc quick)
ENDIF()
