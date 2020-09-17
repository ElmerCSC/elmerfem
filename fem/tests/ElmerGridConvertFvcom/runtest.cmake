include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 16 2 fvcom_test.2dm -out mesh)
RUN_ELMER_TEST()
