include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 angle.grd -partition 4 1 1 -nooverwrite)
RUN_ELMER_TEST()
