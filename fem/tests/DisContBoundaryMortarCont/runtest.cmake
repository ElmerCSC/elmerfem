execute_process(COMMAND ${ELMERGRID_BIN} 1 2 square)
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 square -partdual -metis ${MPIEXEC_NTASKS} 3 -halobc)
RUN_ELMER_TEST()
