INCLUDE(${TEST_SOURCE}/../test_macros.cmake)

FILE(COPY ${BINARY_DIR}/fem/src/modules/FreeSurfaceSolver${SHLEXT} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
FILE(RENAME FreeSurfaceSolver${SHLEXT} FreeSurface1${SHLEXT})
EXECUTE_PROCESS(COMMAND ${ELMERGRID_BIN} 14 2 teterousse.msh -autoclean -order 1.0 0.1 0.01)

# The same problem twice: teterousse3a.sif as shipped, on the legacy FlowSolve, then
# stabilized.sif on IncompressibleNSVec with the equal-order pair. Both carry
# the timers so the logs show what each costs.
EXECUTE_ELMER_SOLVER(stabilized.sif)
IF(NOT EXISTS "TEST.PASSED")
  MESSAGE(FATAL_ERROR
    "the equal-order variant produced no TEST.PASSED at all -- it did not run to "
    "completion. See stabilized.sif-stdout.log and -stderr.log in this directory.")
ENDIF()
FILE(READ "TEST.PASSED" _res)
IF(NOT _res EQUAL "1")
  SET(_cmp "")
  IF(EXISTS "stabilized.sif-stdout.log")
    FILE(STRINGS "stabilized.sif-stdout.log" _lines REGEX "CompareToReferenceSolution")
    STRING(REPLACE ";" "\n  " _cmp "${_lines}")
  ENDIF()
  MESSAGE(FATAL_ERROR
    "the equal-order variant failed its reference norm\n  ${_cmp}")
ENDIF()

RUN_ELMERICE_TEST()
