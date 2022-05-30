include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 quads.grd -nooverwrite)
execute_process(COMMAND ${ELMERSOLVER_BIN} case2d.sif)
RUN_ELMER_TEST()
