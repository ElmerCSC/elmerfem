include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 sinus3d.grd -autoclean -partdual -metiskway ${MPIEXEC_NTASKS} -nooverwrite)
RUN_ELMER_TEST()
