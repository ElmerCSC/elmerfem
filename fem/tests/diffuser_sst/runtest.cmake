include(${TEST_SOURCE}/../test_macros.cmake)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 diffuser -scale 0.127 0.127 0.127)
RUN_ELMER_TEST()
