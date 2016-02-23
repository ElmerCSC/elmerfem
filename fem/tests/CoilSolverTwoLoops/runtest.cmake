execute_process(COMMAND ${ELMERGRID_BIN} 14 2 Loop.msh -autoclean -out TwoLoops -clone 2 1 1 -cloneinds)
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 TwoLoops -partdual -metis ${MPIEXEC_NTASKS} 3 -nooverwrite)
RUN_ELMER_TEST()
