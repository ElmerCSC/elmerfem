include(test_macros)

# Partition the serial mesh in place. ElmerGrid is a no-op for a single task,
# where the solver reads 1381/ directly.
IF(${MPIEXEC_NTASKS} GREATER 1)
  execute_process(COMMAND ${ELMERGRID_BIN} 2 2 1381 -partdual -metis ${MPIEXEC_NTASKS} 3 -nooverwrite)
ENDIF()

RUN_ELMER_TEST()
