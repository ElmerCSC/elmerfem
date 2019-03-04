include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 square -partdual -metiskway ${MPIEXEC_NTASKS} -nooverwrite -halo)
RUN_ELMER_TEST()
