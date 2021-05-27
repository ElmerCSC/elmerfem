include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 angle)
execute_process(COMMAND ${ELMERSOLVER_BIN} part.sif)
RUN_ELMER_TEST()
