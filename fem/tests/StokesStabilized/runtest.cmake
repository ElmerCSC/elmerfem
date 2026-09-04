include(test_macros)

# Equal-order velocity/pressure in IncompressibleNSVec, held stable by the
# "Pressure Stabilization" term instead of by a bubble. Three families, three
# sifs, one ctest entry -- the same idea on quadrilaterals, triangles and
# tetrahedra, so what differs between them is the per-family stabilization
# constant and nothing else.
#
# The quadrilaterals and triangles come from the one Step.grd, the triangles by
# splitting the same quadrilaterals, which is what makes those two comparable.
# They go to different output directories because they would otherwise overwrite
# each other. The tetrahedra are a separate mesh: nothing in the suite generates
# tetrahedra for a flow case, so tet.msh is borrowed from ElasticStabilized.

execute_process(COMMAND ${ELMERGRID_BIN} 1 2 Step.grd -out stepquad)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 Step.grd -triangles -out steptri)
execute_process(COMMAND ${ELMERGRID_BIN} 14 2 tet.msh -out tetmesh)

# Only the closing RUN_ELMER_TEST() inspects TEST.PASSED and every run
# overwrites it, so the earlier families are checked here by hand -- which also
# lets a failure name the family that broke.
MACRO(CHECK_FAMILY SIF WHAT)
  FILE(READ "TEST.PASSED" _res)
  IF(NOT _res EQUAL "1")
    SET(_cmp "")
    IF(EXISTS "${SIF}-stdout.log")
      FILE(STRINGS "${SIF}-stdout.log" _lines REGEX "CompareToReferenceSolution")
      STRING(REPLACE ";" "\n  " _cmp "${_lines}")
    ENDIF()
    MESSAGE(FATAL_ERROR
      "pressure stabilization on ${WHAT} failed its reference norm\n  ${_cmp}")
  ENDIF()
ENDMACRO()

EXECUTE_ELMER_SOLVER(tri.sif)
CHECK_FAMILY(tri.sif "triangles")

EXECUTE_ELMER_SOLVER(tet.sif)
CHECK_FAMILY(tet.sif "tetrahedra")

# The quadrilateral goes through RUN_ELMER_TEST(), which reads the sif named in
# ELMERSOLVER_STARTINFO and reports the timing the harness expects.
RUN_ELMER_TEST()
