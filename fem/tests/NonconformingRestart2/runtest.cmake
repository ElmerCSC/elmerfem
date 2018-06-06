execute_process(COMMAND ${ELMERGRID_BIN} 1 2 mesh_coarse)
execute_process(COMMAND ${ELMERSOLVER_BIN} coarse.sif)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 mesh_fine)
RUN_ELMER_TEST()
