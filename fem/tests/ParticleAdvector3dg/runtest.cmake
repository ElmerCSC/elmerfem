include(test_macros)

# The same problem twice: case.sif on FlowSolve, then stabilized.sif on
# IncompressibleNSVec with the equal-order pair. Both baselines are already
# equal order -- the shipped case sets "Stabilize = True" -- so this compares two
# implementations of one scheme, and solver 2 then advects particles through the
# result, checking the flow is usable downstream and not just close in a norm.
#
# Only the closing RUN_ELMER_TEST() inspects TEST.PASSED and every run overwrites
# it, so the first is checked here by hand.
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
