include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 Step_v2f)
RUN_ELMER_TEST()
