include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 9 2 testmesh.mphtxt -out mesh)
RUN_ELMER_TEST()
