IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE( PerfectGasCompress_vec.sif PerfectGasCompress_vec.sif COPYONLY)

  file(COPY PerfectGasCompress_vec.grd ELMERSOLVER_STARTINFO  DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_ELMER_TEST(PerfectGasCompress_vec LABELS matc quick)
ENDIF()
