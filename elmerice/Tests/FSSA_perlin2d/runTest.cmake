INCLUDE(${TEST_SOURCE}/../test_macros.cmake)


# The same problem twice with the method switched in between: Stokes_prognostic_nonlWeert.sif as
# shipped, on the MINI element, then stabilized.sif, the equal-order pair whose
# pressure mode is held by "Pressure Stabilization". Both carry the timers, so
# the logs show what the cheaper element costs and saves.
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
