execute_process(COMMAND ${ELMERGRID_BIN} 1 2 winkel -metis ${MPIEXEC_NTASKS} 3 -nooverwrite)

RUN_ELMER_TEST()
