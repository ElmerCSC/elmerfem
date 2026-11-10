include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 14 2 2layer.msh -autoclean -scale 1000 1000 1)
RUN_ELMER_TEST()
