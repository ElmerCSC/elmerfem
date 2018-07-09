include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 winkel -partdual -metisrec ${MPIEXEC_NTASKS} -nooverwrite)

RUN_ELMER_TEST()
