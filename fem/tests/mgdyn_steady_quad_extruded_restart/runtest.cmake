include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 quads.grd -nooverwrite)
EXECUTE_ELMER_SOLVER(case2d.sif)
RUN_ELMER_TEST()
