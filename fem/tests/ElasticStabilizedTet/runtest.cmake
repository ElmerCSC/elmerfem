include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 14 2 tet.msh -out tetmesh)
RUN_ELMER_TEST()
