include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 14 2 actuator.msh)
RUN_ELMER_TEST()
