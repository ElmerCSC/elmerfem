include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 mesh -triangles)
RUN_ELMER_TEST()
