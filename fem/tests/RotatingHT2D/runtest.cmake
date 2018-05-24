execute_process(COMMAND ${ELMERGRID_BIN} 14 2 teeth.msh -autoclean )
RUN_ELMER_TEST()
