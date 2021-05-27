include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 cubes.grd)
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 cubes -partdual -metis ${MPIEXEC_NTASKS} 3)

RUN_ELMER_TEST()
