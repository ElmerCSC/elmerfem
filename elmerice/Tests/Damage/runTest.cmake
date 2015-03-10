INCLUDE(${TEST_SOURCE}/../test_macros.cmake)

FILE(COPY ${TEST_SOURCE}/../../../../install/share/elmersolver/lib/FreeSurfaceSolver.so DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
FILE(RENAME FreeSurfaceSolver.so FreeSurfaceSolver1.so)

EXECUTE_PROCESS(COMMAND ${ELMERGRID_BIN} 1 2 mesh.grd)

RUN_ELMERICE_TEST()

