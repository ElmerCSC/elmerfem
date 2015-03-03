INCLUDE(${TEST_SOURCE}/../test_macros.cmake)

FILE(COPY ${TEST_SOURCE}/../../../../install/share/elmersolver/lib/FreeSurfaceSolver.so DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
FILE(RENAME FreeSurfaceSolver.so FreeSurface1.so)

EXECUTE_PROCESS(COMMAND ${ELMERGRID_BIN} 14 2 teterousse.msh -autoclean -order 1.0 0.1 0.01)

RUN_ELMERICE_TEST()
