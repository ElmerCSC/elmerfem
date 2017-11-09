execute_process(COMMAND ${ELMERGRID_BIN} 14 2 Loop.msh -autoclean)
RUN_ELMER_TEST()
