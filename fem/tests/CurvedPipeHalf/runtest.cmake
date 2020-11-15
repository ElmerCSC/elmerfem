include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 curved_pipe_half)
RUN_ELMER_TEST()
