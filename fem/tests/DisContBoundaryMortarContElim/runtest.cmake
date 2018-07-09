include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 square)
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 square -partdual -metis ${MPIEXEC_NTASKS} 3 -halobc -nooverwrite)
RUN_ELMER_TEST()
