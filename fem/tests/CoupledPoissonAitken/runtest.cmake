include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 angle)

# Partition the serial mesh in place. Not needed for a single task, where the
# solver reads angle/ directly.
IF(${MPIEXEC_NTASKS} GREATER 1)
  execute_process(COMMAND ${ELMERGRID_BIN} 2 2 angle -metis ${MPIEXEC_NTASKS} -removeunused)
ENDIF()

RUN_ELMER_TEST()
