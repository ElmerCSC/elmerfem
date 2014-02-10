
include(${TEST_SOURCE}/../test_macros.cmake)

execute_process(COMMAND ${ELMERGRID_BIN} 1 2 1d -decimals 20)

RUN_ELMER_TEST()