include(test_macros)

# Equal-order velocity/pressure in IncompressibleNSVec, held stable by the
# "Pressure Stabilization" term instead of by a bubble. Two families, two sifs,
# one ctest entry -- the same problem on a quadrilateral and a triangle mesh, so
# what differs between them is the per-family stabilization constant and nothing
# else.
#
# Both meshes come from the one Step.grd, the triangles by splitting the same
# quadrilaterals, which is what makes the two comparable. They go to different
# output directories because they would otherwise overwrite each other.

execute_process(COMMAND ${ELMERGRID_BIN} 1 2 Step.grd -out stepquad)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 Step.grd -triangles -out steptri)

# Only the closing RUN_ELMER_TEST() inspects TEST.PASSED and every run
# overwrites it, so the earlier family is checked here by hand -- which also
# lets a failure name the family that broke. Same shape as ElasticStabilized,
# whose comments explain why the comparison lines are echoed with the message.
EXECUTE_ELMER_SOLVER(tri.sif)
FILE(READ "TEST.PASSED" _res)
IF(NOT _res EQUAL "1")
  SET(_cmp "")
  IF(EXISTS "tri.sif-stdout.log")
    FILE(STRINGS "tri.sif-stdout.log" _lines REGEX "CompareToReferenceSolution")
    STRING(REPLACE ";" "\n  " _cmp "${_lines}")
  ENDIF()
  MESSAGE(FATAL_ERROR
    "pressure stabilization on triangles failed its reference norm\n  ${_cmp}")
ENDIF()

# The quadrilateral goes through RUN_ELMER_TEST(), which reads the sif named in
# ELMERSOLVER_STARTINFO and reports the timing the harness expects.
RUN_ELMER_TEST()
