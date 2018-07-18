include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 mesh_fine -partition ${MPIEXEC_NTASKS} 1 1 -nooverwrite)
RUN_ELMER_TEST()
