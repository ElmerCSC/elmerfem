include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 sinus0)
execute_process(COMMAND ${ELMERGRID_BIN} extrude.eg)
RUN_ELMER_TEST()
