include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 slide_squares.grd)
RUN_ELMER_TEST()
