include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 shoebox_prisms)
RUN_ELMER_TEST()
