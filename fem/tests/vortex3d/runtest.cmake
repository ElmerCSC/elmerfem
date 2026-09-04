include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 Step3d)

# The same 3D problem twice. Unusually, BOTH runs are equal-order stabilised:
# Step3d.sif takes "Stabilization Method = Stabilized" through the legacy
# FlowSolve, stabilized.sif takes "Pressure Stabilization" through
# IncompressibleNSVec. So this compares two independent implementations of one
# scheme rather than two different discretisations, which is what makes it worth
# keeping as a 3D check.
#
# Only the closing RUN_ELMER_TEST() inspects TEST.PASSED and every run
# overwrites it, so the first is checked here by hand.
EXECUTE_ELMER_SOLVER(stabilized.sif)
FILE(READ "TEST.PASSED" _res)
IF(NOT _res EQUAL "1")
  SET(_cmp "")
  IF(EXISTS "stabilized.sif-stdout.log")
    FILE(STRINGS "stabilized.sif-stdout.log" _lines REGEX "CompareToReferenceSolution")
    STRING(REPLACE ";" "\n  " _cmp "${_lines}")
  ENDIF()
  MESSAGE(FATAL_ERROR
    "the IncompressibleNSVec equal-order variant failed its reference norm\n  ${_cmp}")
ENDIF()

RUN_ELMER_TEST()
