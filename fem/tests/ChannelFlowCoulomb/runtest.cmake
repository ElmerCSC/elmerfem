include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 rect)
EXECUTE_ELMER_SOLVER(taylorhood.sif)
IF(NOT EXISTS "TEST.PASSED")
  MESSAGE(FATAL_ERROR
    "the Taylor-Hood variant produced no TEST.PASSED at all -- it did not run to completion. "
    "See taylorhood.sif-stdout.log and -stderr.log in this directory.")
ENDIF()
FILE(READ "TEST.PASSED" _res)
IF(NOT _res EQUAL "1")
  SET(_cmp "")
  IF(EXISTS "taylorhood.sif-stdout.log")
    FILE(STRINGS "taylorhood.sif-stdout.log" _lines REGEX "CompareToReferenceSolution")
    STRING(REPLACE ";" "\n  " _cmp "${_lines}")
  ENDIF()
  MESSAGE(FATAL_ERROR
    "the Taylor-Hood P2/P1 variant failed its reference norm\n  ${_cmp}")
ENDIF()

RUN_ELMER_TEST()