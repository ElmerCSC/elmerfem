include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 14 2 cubes_inside_cube.msh -autoclean -partdual -metiskway ${MPIEXEC_NTASKS} -nooverwrite)
RUN_ELMER_TEST()
