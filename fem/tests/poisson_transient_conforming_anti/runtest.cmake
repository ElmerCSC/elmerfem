include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 square.grd -connect 1 2 3 4 -partdual -metiskway ${MPIEXEC_NTASKS} -nooverwrite)
RUN_ELMER_TEST()
