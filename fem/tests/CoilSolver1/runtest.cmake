include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 14 2 Curve.msh -autoclean)
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 Curve -partdual -metis ${MPIEXEC_NTASKS} 3 -nooverwrite)
RUN_ELMER_TEST()
