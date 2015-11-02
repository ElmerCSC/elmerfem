execute_process(COMMAND ${ELMERGRID_BIN} 1 2 angle.grd)
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 angle -out angles -clone 1 1 3 -clonesize 0 0 0.5)
RUN_ELMER_TEST()
