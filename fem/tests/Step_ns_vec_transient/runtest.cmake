include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 Step)

# The same transient condensed-MINI problem twice, differing only in how many
# nonlinear iterations each timestep is allowed. case_singleiter.sif gets one and
# case.sif gets the floor of two that the solver's _Init asks for, which is the
# minimum that makes the condensed bubble coefficients consistent with the nodal
# solution they were condensed against. The two norms are 0.85 % apart, and that
# gap is the whole point: if the floor stops being applied, case.sif drifts onto
# case_singleiter.sif's answer and fails its own reference norm.
#
# Only the closing RUN_ELMER_TEST() inspects TEST.PASSED and every run
# overwrites it, so the first is checked here by hand -- the same pattern as
# Step_stokes_vec. That also lets a failure name the variant that broke.
EXECUTE_ELMER_SOLVER(case_singleiter.sif)
IF(NOT EXISTS "TEST.PASSED")
  MESSAGE(FATAL_ERROR
    "the single-iteration variant produced no TEST.PASSED at all -- it did not run to completion. "
    "See case_singleiter.sif-stdout.log and -stderr.log in this directory.")
ENDIF()
FILE(READ "TEST.PASSED" _res)
IF(NOT _res EQUAL "1")
  SET(_cmp "")
  IF(EXISTS "case_singleiter.sif-stdout.log")
    FILE(STRINGS "case_singleiter.sif-stdout.log" _lines REGEX "CompareToReferenceSolution")
    STRING(REPLACE ";" "\n  " _cmp "${_lines}")
  ENDIF()
  MESSAGE(FATAL_ERROR
    "the single-iteration variant failed its reference norm\n  ${_cmp}")
ENDIF()

RUN_ELMER_TEST()
