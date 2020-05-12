include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 squares.grd)
execute_process(COMMAND ${ELMERGRID_BIN} squares2cubes.eg)

RUN_ELMER_TEST()
