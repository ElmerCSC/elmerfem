include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 14 2 toroid2d.msh -autoclean )
RUN_ELMER_TEST()
