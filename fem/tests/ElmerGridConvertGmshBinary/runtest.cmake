include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 14 2 binary.msh -autoclean -out mesh)
RUN_ELMER_TEST()
