include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 rect.grd -partdual -metiskway ${MPIEXEC_NTASKS})
RUN_ELMER_TEST()
