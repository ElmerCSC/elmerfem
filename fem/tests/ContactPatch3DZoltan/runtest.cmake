include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 14 2 cubes.msh)
RUN_ELMER_TEST()
