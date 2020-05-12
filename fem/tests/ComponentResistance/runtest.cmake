include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 blocks)
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 blocks -partdual -metis ${MPIEXEC_NTASKS} 3 -nooverwrite)
RUN_ELMER_TEST()
