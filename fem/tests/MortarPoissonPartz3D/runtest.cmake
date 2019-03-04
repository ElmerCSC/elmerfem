include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 14 2 mesh.msh -autoclean)
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 mesh -partdual -metis ${MPIEXEC_NTASKS} 3 -connect 4 10 -partconnect 3 -haloz -nooverwrite)

RUN_ELMER_TEST()
