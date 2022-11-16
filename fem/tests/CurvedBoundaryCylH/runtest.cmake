include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 14 2 cylinder_in_cube.msh -autoclean)
RUN_ELMER_TEST()
