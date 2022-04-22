include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 cylinders.grd)
RUN_ELMER_TEST()
