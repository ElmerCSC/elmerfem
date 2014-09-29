include(${TEST_SOURCE}/../test_macros.cmake)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 helix.grd)
RUN_ELMER_TEST_NEW()
