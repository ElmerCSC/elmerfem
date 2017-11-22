execute_process(COMMAND ${ELMERGRID_BIN} 2 2 umagnet -partdual -metis ${MPIEXEC_NTASKS} -nooverwrite)
RUN_ELMER_TEST()
