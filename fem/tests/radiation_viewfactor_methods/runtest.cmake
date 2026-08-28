include(test_macros)

# The method keywords are read once per run, in GetCartParameters, and the
# view factors are computed by a separate binary before the solver starts.
# So a single Elmer run cannot cross the method matrix; this test is several
# runs instead, which the ctest harness is happy to do.
#
# Note that only the closing RUN_ELMER_TEST() inspects TEST.PASSED, and every
# run overwrites it.  Intermediate runs are therefore checked here by hand,
# which also lets a failure name the corner that broke.

# Each corner's view factors are kept so the corners can be compared against
# each other afterwards, not only against a norm.
MACRO(RUN_CORNER SIF)
  execute_process(COMMAND ${VIEWFACTORS_BIN} ${SIF}.sif
    OUTPUT_FILE ${SIF}-viewfactors.log ERROR_FILE ${SIF}-viewfactors-err.log)
  configure_file(box_in_box/ViewFactors.dat vf_${SIF}.dat COPYONLY)
ENDMACRO()

# A shared reference norm asserts that the corners agree.  On its own that is
# not enough: if a method silently fell back to another, the corners would
# agree perfectly and the test would pass while testing nothing.  So also
# require that every corner produced a *different* view factor matrix, which
# is what proves each one really took the path it asked for.
MACRO(REQUIRE_DIFFERENT A B)
  execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files vf_${A}.dat vf_${B}.dat
    RESULT_VARIABLE _differs)
  IF(_differs EQUAL 0)
    MESSAGE(FATAL_ERROR
      "corners '${A}' and '${B}' produced identical view factors: "
      "one of them is not taking the path it asks for")
  ENDIF()
ENDMACRO()

MACRO(CHECK_CORNER WHAT)
  FILE(READ "TEST.PASSED" _res)
  IF(NOT _res EQUAL "1")
    MESSAGE(FATAL_ERROR "view factor method '${WHAT}' failed its reference norm")
  ENDIF()
ENDMACRO()

execute_process(COMMAND ${ELMERGRID_BIN} 1 2 box_in_box.grd)
execute_process(COMMAND ${RADIATORS_BIN})

# ---------------------------------------------------------------------------
# The shaft cull is a pure speed optimisation: it replaces one tree traversal
# per ray with one shaft build per pair and must not change a single bit of
# the answer.  That is an exact invariant, so check it exactly rather than
# through a norm tolerance.  At eight rays over a six element shadow mesh the
# automatic rule leaves the cull off, so the two runs really do take different
# paths through the code.
# ---------------------------------------------------------------------------
execute_process(COMMAND ${VIEWFACTORS_BIN} culloff.sif
  OUTPUT_FILE culloff-viewfactors.log ERROR_FILE culloff-viewfactors-err.log)
configure_file(box_in_box/ViewFactors.dat vf_culloff.dat COPYONLY)

execute_process(COMMAND ${VIEWFACTORS_BIN} cullon.sif
  OUTPUT_FILE cullon-viewfactors.log ERROR_FILE cullon-viewfactors-err.log)
configure_file(box_in_box/ViewFactors.dat vf_cullon.dat COPYONLY)

execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files vf_culloff.dat vf_cullon.dat
  RESULT_VARIABLE CULL_DIFFERS)
IF(NOT CULL_DIFFERS EQUAL 0)
  MESSAGE(FATAL_ERROR
    "the shaft cull changed the view factors; it is required to be bit identical")
ENDIF()

# ---------------------------------------------------------------------------
# The four corners of the matrix: inner integration (quadrature or the closed
# form contour formula) against shadowing (ray casting or clipping).  They
# share one reference norm, so what is asserted is that the corners agree.
# ---------------------------------------------------------------------------
RUN_CORNER(numrays)
EXECUTE_ELMER_SOLVER(numrays.sif)
CHECK_CORNER("numerical + ray casting")

RUN_CORNER(cfrays)
EXECUTE_ELMER_SOLVER(cfrays.sif)
CHECK_CORNER("closed form + ray casting")

RUN_CORNER(cfclip)
EXECUTE_ELMER_SOLVER(cfclip.sif)
CHECK_CORNER("closed form + clipping")

# The last corner goes through RUN_ELMER_TEST(), which reads the sif named in
# ELMERSOLVER_STARTINFO and reports the timing the harness expects.
RUN_CORNER(numclip)

# Every corner must differ from every other one.
REQUIRE_DIFFERENT(numrays cfrays)
REQUIRE_DIFFERENT(numrays cfclip)
REQUIRE_DIFFERENT(numrays numclip)
REQUIRE_DIFFERENT(cfrays  cfclip)
REQUIRE_DIFFERENT(cfrays  numclip)
REQUIRE_DIFFERENT(cfclip  numclip)

RUN_ELMER_TEST()
