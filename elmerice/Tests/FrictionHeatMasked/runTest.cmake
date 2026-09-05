INCLUDE(${TEST_SOURCE}/../test_macros.cmake)

EXECUTE_PROCESS(COMMAND ${ELMERGRID_BIN} 1 2 cube.grd)

# The same problem twice: test.sif as shipped, on the legacy FlowSolve, then
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
