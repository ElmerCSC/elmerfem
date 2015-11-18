INCLUDE(${TEST_SOURCE}/../test_macros.cmake)

FILE(COPY ${BINARY_DIR}/fem/src/modules/FreeSurfaceSolver${SHLEXT} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
FILE(RENAME FreeSurfaceSolver${SHLEXT} FreeSurface1${SHLEXT})

EXECUTE_PROCESS(COMMAND ${ELMERGRID_BIN} 14 2 teterousse.msh -autoclean -order 1.0 0.1 0.01)

RUN_ELMERICE_TEST()
