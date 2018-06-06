execute_process(COMMAND ${ELMERGRID_BIN} 1 2 angle_coarse)
execute_process(COMMAND ${ELMERSOLVER_BIN} coarse.sif)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 angle_fine)
RUN_ELMER_TEST()
