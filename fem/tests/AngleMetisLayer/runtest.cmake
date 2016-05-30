execute_process(COMMAND ${ELMERGRID_BIN} 1 2 angle.grd -partdual -metis ${MPIEXEC_NTASKS} 3 -connect 1 2 -nooverwrite)

RUN_ELMER_TEST()
