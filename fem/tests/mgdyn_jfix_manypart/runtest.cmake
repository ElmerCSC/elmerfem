include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 cube -metis 8 3 -nooverwrite)
RUN_ELMER_TEST()
