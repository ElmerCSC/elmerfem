IF(WITH_LUA)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(case.sif case.sif COPYONLY)

  file(COPY ELMERSOLVER_STARTINFO square.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_ELMER_TEST(TwelveSolvers2D LABELS lua)
ENDIF(WITH_LUA)
