include(test_macros)

# One element of each 3D family -- tetrahedron, pyramid, prism, hexahedron --
# assembled by IncompressibleNSVec's equal-order pair. A smoke test that every
# family the solver can meet actually assembles; see case.sif for the two
# defects that made it worth having, and for why the MINI variant is not run
# alongside it.
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