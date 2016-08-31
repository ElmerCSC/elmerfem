execute_process(COMMAND ${ELMERGRID_BIN} 1 2 angle -partdual -metis ${MPIEXEC_NTASKS} 3 -nooverwrite)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 square -partdual -metis ${MPIEXEC_NTASKS} 3 -nooverwrite)
RUN_ELMER_TEST()
