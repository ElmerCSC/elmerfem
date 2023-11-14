include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 14 2 squares_inside_square.msh -autoclean)
RUN_ELMER_TEST()
