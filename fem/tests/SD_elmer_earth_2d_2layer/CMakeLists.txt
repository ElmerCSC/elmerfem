IF(WITH_LUA)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(2layer.sif 2layer.sif COPYONLY)

  file(COPY 2layer.msh timeload.lua ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_ELMER_TEST(SD_elmer_earth_2d_2layer)
  ADD_ELMER_LABEL(SD_elmer_earth_2d_2layer elasticity viscoelastic transient vtu lua serendipity)
ENDIF(WITH_LUA)
