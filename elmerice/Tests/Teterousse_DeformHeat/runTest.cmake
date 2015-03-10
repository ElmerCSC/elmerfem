INCLUDE(${TEST_SOURCE}/../test_macros.cmake)

EXECUTE_PROCESS(COMMAND ${ELMERGRID_BIN} 14 2 teterousse.msh -autoclean -order 1.0 0.1 0.01)

RUN_ELMERICE_TEST()
