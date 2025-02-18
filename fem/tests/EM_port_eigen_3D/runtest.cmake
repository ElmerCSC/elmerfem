include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 port.grd)
execute_process(COMMAND ${ELMERGRID_BIN} extrude.eg)
RUN_ELMER_TEST()
