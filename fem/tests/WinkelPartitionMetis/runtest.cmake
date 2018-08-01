include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 winkel.grd -metis 8 -nooverwrite)

RUN_ELMER_TEST()
