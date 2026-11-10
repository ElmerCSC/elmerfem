include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 shoebox_tets -partdual -metiskway ${MPIEXEC_NTASKS} -nooverwrite)

RUN_ELMER_TEST()
