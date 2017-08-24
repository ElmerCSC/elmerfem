INCLUDE(${TEST_SOURCE}/../test_macros.cmake)

EXECUTE_PROCESS(COMMAND ${ELMERGRID_BIN} 14 2 square.msh -autoclean)

RUN_ELMERICE_TEST()
