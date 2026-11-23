include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 rect)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 rect -partition 1 2 1 -nooverwrite)
RUN_ELMER_TEST()
