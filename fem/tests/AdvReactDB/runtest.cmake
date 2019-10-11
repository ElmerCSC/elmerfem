include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 square -metis ${MPIEXEC_NTASKS} 3 -halo -nooverwrite)
RUN_ELMER_TEST()
