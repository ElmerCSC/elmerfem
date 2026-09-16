INCLUDE(${TEST_SOURCE}/../test_macros.cmake)

IF(WIN32 AND WITH_MPI)
  EXECUTE_PROCESS(COMMAND "${MPIEXEC}" ${MPIEXEC_NUMPROC_FLAG} 1 hostname)
ENDIF()

RUN_ELMERICE_TEST()

FILE(STRINGS "test-stdout.log" _nonconvergence
  REGEX "Coupled system did not converge")
IF(_nonconvergence)
  MESSAGE(FATAL_ERROR "Permafrost_Lunardini had a coupled nonconvergence warning")
ENDIF()
