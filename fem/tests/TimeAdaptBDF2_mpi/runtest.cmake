include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 sq)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 sq -partition 2 1 1 -nooverwrite)
RUN_ELMER_TEST()
