INCLUDE(${TEST_SOURCE}/../test_macros.cmake)

EXECUTE_PROCESS(COMMAND ${ELMERGRID_BIN} 1 2 square-tri.grd)

# On Windows with MSMPI, warm up the process manager (smpd) before launching
# ElmerSolver, otherwise the first mpiexec call intermittently produces no output.
IF(WIN32 AND WITH_MPI)
  EXECUTE_PROCESS(COMMAND "${MPIEXEC}" ${MPIEXEC_NUMPROC_FLAG} 1 hostname)
ENDIF()

RUN_ELMERICE_TEST()

