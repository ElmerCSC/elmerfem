include(${TEST_SOURCE}/../test_macros.cmake)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 square)
EXECUTE_ELMER_SOLVER(init.sif)
RUN_ELMER_TEST()
