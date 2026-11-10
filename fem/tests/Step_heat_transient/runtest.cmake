include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 Step)

# The same transient convection-heat problem run twice with HeatSolveVec's
# method switched in between: Step_heat_transient.sif as shipped, on the
# condensed p-bubble, and stabilized.sif, the equal-order pair stabilized by
# SUPG ("Stabilize = True") instead. Only the closing RUN_ELMER_TEST()
# inspects TEST.PASSED and every run overwrites it, so the first is checked
# here by hand -- that also lets a failure name the variant that broke rather
# than just the test.
EXECUTE_ELMER_SOLVER(stabilized.sif)
IF(NOT EXISTS "TEST.PASSED")
  MESSAGE(FATAL_ERROR
    "the SUPG variant produced no TEST.PASSED at all -- it did not run to completion. "
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
    "the SUPG variant failed its reference norm\n  ${_cmp}")
ENDIF()

RUN_ELMER_TEST()
