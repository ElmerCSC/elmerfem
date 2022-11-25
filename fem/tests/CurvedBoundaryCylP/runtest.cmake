include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 14 2 cylinder_in_cylinder.msh -autoclean -partdual -metiskway ${MPIEXEC_NTASKS} -nooverwrite)
RUN_ELMER_TEST()
