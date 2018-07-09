include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 circle_in_square -partdual -metis ${MPIEXEC_NTASKS} 3 -connect 5 6 7 8)

RUN_ELMER_TEST()
