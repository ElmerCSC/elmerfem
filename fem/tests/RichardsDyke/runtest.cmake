include(${TEST_SOURCE}/../test_macros.cmake)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 land_dyke)
RUN_ELMER_TEST_NEW()
