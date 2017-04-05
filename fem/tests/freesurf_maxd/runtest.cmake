execute_process(COMMAND ${ELMERGRID_BIN} 1 2 Mesh -partdual -metis ${MPIEXEC_NTASKS} 4 -nooverwrite)
RUN_ELMER_TEST()
