execute_process(COMMAND ${ELMERGRID_BIN} 2 2 cylinder -partdual -metis ${MPIEXEC_NTASKS} 3 -nooverwrite)
RUN_ELMER_TEST()
