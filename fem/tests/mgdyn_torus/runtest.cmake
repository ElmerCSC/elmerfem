include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 torus2 -partdual -metis ${MPIEXEC_NTASKS} 3 -nooverwrite)
RUN_ELMER_TEST()
