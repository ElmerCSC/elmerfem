INCLUDE(${TEST_SOURCE}/../test_macros.cmake)


# The same problem twice with the method switched in between: Stokes_prognostic_nonlWeert.sif as
# shipped, on the MINI element, then stabilized.sif, the equal-order pair whose
# pressure mode is held by "Pressure Stabilization". Both carry the timers, so
# the logs show what the cheaper element costs and saves.
#
# ELMER_MODULES_PATH has to be set by hand here. RUN_ELMERICE_TEST sets it, but
# the bare EXECUTE_ELMER_SOLVER in the same file does not, and these cases load
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
