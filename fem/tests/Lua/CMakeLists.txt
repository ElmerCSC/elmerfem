IF(WITH_LUA)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE( case.sif case.sif COPYONLY)
  CONFIGURE_FILE( defaults.lua defaults.lua COPYONLY)
  ADD_ELMERTEST_MODULE(Lua_evol ModelPDEevol ModelPDEevol.F90)
  ADD_ELMERTEST_MODULE(Lua_thr ModelPDEthreaded ModelPDEthreaded.F90)

  file(COPY square.grd ELMERSOLVER_STARTINFO ModelPDEthreaded.F90 ModelPDEevol.F90 DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_ELMER_TEST(Lua LABELS lua quick)
ENDIF(WITH_LUA)
