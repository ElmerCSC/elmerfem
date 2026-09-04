INCLUDE(${TEST_SOURCE}/../test_macros.cmake)

EXECUTE_PROCESS(COMMAND ${ELMERGRID_BIN} 14 2 teterousse.msh -autoclean -order 1.0 0.1 0.0)

# The same problem twice: test.sif as shipped, on the legacy FlowSolve, then
# stabilized.sif on IncompressibleNSVec with the equal-order pair. Both carry
# the timers so the logs show what each costs.
#
# ELMER_MODULES_PATH must be set by hand: RUN_ELMERICE_TEST sets it, the bare
# EXECUTE_ELMER_SOLVER in the same macro file does not, and these load
# ElmerIceSolvers.
SET(ENV{ELMER_MODULES_PATH} "${BINARY_DIR}/elmerice/Solvers:${BINARY_DIR}/elmerice/Solvers/GridDataReader:${BINARY_DIR}/elmerice/Solvers/ScatteredDataInterpolator:${BINARY_DIR}/elmerice/Solvers/MeshAdaptation_2D:${BINARY_DIR}/elmerice/UserFunctions")

EXECUTE_ELMER_SOLVER(stabilized.sif)
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
