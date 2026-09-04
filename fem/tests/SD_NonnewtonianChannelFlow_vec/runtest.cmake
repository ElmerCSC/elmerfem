include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 square)

# The same problem twice with the method switched in between: case.sif as
# shipped, on the MINI element, and stabilized.sif, the equal-order pair whose
# pressure mode is held by "Pressure Stabilization" instead of by a bubble. Both
# carry "Solver Timing", so the logs show what the cheaper element costs and
# saves -- assembly, linear system and total, per variant.
#
# Only the closing RUN_ELMER_TEST() inspects TEST.PASSED and every run
# overwrites it, so the first is checked here by hand. That also lets a failure
# name the variant that broke rather than just the test.
EXECUTE_ELMER_SOLVER(stabilized.sif)
IF(NOT EXISTS "TEST.PASSED")
  MESSAGE(FATAL_ERROR
    "the equal-order variant produced no TEST.PASSED at all -- it did not run to completion. "
    "See stabilized.sif-stdout.log and -stderr.log in this directory.")
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

RUN_ELMER_TEST()
