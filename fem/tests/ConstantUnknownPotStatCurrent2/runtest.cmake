include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 angle -metiskway ${MPIEXEC_NTASKS} -nooverwrite)

RUN_ELMER_TEST()
