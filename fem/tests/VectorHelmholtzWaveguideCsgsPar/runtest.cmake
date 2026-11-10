include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 shoebox_prisms -metis ${MPIEXEC_NTASKS} 3 -partdual -nooverwrite)
RUN_ELMER_TEST()
