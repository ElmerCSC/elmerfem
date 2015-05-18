INCLUDE(${TEST_SOURCE}/../test_macros.cmake)

EXECUTE_PROCESS(COMMAND ${ELMERGRID_BIN} 1 2 cube.grd -boundbound 1 5 7 -boundbound 2 5 7 -boundbound 3 5 7 -boundbound 4 5 8)

RUN_ELMERICE_TEST()
