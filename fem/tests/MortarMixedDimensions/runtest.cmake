include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 square.grd)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 cube.grd)
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 cube -in square -unite -out mesh)

RUN_ELMER_TEST()
