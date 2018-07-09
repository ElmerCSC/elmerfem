include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 angle.grd -out angles -clone 1 1 3 -clonesize 0 0 0.5 -nooverwrite)
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 angles -partition 1 1 ${MPIEXEC_NTASKS} -nooverwrite)
RUN_ELMER_TEST()
