execute_process(COMMAND ${ELMERGRID_BIN} 14 2 BeamFlow.msh -autoclean )
RUN_ELMER_TEST()
