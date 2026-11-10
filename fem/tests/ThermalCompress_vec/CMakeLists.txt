IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE( ThermalCompress_vec.sif ThermalCompress_vec.sif COPYONLY)

  file(COPY ThermalCompress_vec.grd ELMERSOLVER_STARTINFO  DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_ELMER_TEST(ThermalCompress_vec LABELS matc quick)
ENDIF()
