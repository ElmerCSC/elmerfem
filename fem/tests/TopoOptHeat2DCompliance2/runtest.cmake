include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 tree2.grd )
RUN_ELMER_TEST()
