include(test_macros)

# Five element families, five sifs, one ctest entry. Each sif is the same
# problem on a different mesh, so what is being compared across them is the
# per-family stabilization constant and nothing else.
#
# Only the closing RUN_ELMER_TEST() inspects TEST.PASSED, and every run
# overwrites it -- so the earlier families are checked here by hand, which also
# lets a failure name the family that broke instead of just the test.

execute_process(COMMAND ${ELMERGRID_BIN} 1 2 square.grd)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 trimesh.grd)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 cube.grd)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 prismmesh.grd)
execute_process(COMMAND ${ELMERGRID_BIN} 14 2 tet.msh -out tetmesh)

# EXECUTE_ELMER_SOLVER redirects the run's output to <sif>-stdout.log, which a
# CI log never shows -- so a bare FATAL_ERROR here reports that a family failed
# and nothing about by how much, which is useless for the one thing these norms
# are sensitive to. Echo the comparison lines with the message.
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

EXECUTE_ELMER_SOLVER(quad.sif)
CHECK_FAMILY(quad.sif "quadrilaterals")

EXECUTE_ELMER_SOLVER(tri.sif)
CHECK_FAMILY(tri.sif "triangles")

EXECUTE_ELMER_SOLVER(hex.sif)
CHECK_FAMILY(hex.sif "hexahedra")

EXECUTE_ELMER_SOLVER(prism.sif)
CHECK_FAMILY(prism.sif "prisms")

# The tetrahedron goes through RUN_ELMER_TEST(), which reads the sif named in
# ELMERSOLVER_STARTINFO and reports the timing the harness expects.
RUN_ELMER_TEST()
