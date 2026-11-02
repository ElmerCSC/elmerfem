include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 cubes)
RUN_ELMER_TEST()
