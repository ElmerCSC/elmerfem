include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 mesh_fine)
RUN_ELMER_TEST()
