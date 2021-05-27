include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 land_dyke)
RUN_ELMER_TEST()
