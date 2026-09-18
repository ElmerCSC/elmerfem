include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 PerfectGasCompress_vec_transient)
RUN_ELMER_TEST()
