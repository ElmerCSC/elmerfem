include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 mesh_tet -metis 2 3 -nooverwrite)
RUN_ELMER_TEST()
