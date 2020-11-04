include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 8 2 box.unv -autoclean -out mesh)
RUN_ELMER_TEST()
