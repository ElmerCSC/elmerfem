include(test_macros)

# The method keywords are read once per run, in GetCartParameters, and the
# view factors are computed by a separate binary before the solver starts.
# So a single Elmer run cannot cross the method matrix; this test is several
# runs instead, which the ctest harness is happy to do.
#
# Note that only the closing RUN_ELMER_TEST() inspects TEST.PASSED, and every
# run overwrites it.  Intermediate runs are therefore checked here by hand,
# which also lets a failure name the corner that broke.

# Run one of the factor binaries, and insist that it actually ran.  Without the
# RESULT_VARIABLE check a binary that never starts -- the Windows loader
# refusing it for a DLL it cannot find is the way this happens -- looks exactly
# like a binary that ran and wrote nothing, and the test then fails several
# lines later on a missing file with no hint of why.  The log is redirected, so
# echo it back on failure or the message would be the only thing anyone sees.
MACRO(RUN_TOOL BIN LOGBASE)
  execute_process(COMMAND ${BIN} ${ARGN}
    OUTPUT_FILE ${LOGBASE}.log ERROR_FILE ${LOGBASE}-err.log
    RESULT_VARIABLE _rc)
  IF(NOT _rc EQUAL 0)
    SET(_out "")
    IF(EXISTS ${LOGBASE}.log)
      FILE(READ ${LOGBASE}.log _out)
    ENDIF()
    SET(_err "")
    IF(EXISTS ${LOGBASE}-err.log)
      FILE(READ ${LOGBASE}-err.log _err)
    ENDIF()
    MESSAGE(FATAL_ERROR
      "'${BIN} ${ARGN}' failed with '${_rc}'\n"
      "--- stdout ---\n${_out}\n--- stderr ---\n${_err}")
  ENDIF()
ENDMACRO()

# Radiators is a thin wrapper: it re-spells its own argv[0] as ViewFactors and
# passes its last argument through, so it does honour a sif name -- but only if
# it is given one.  Called bare it hands ViewFactors nothing but the -radiators
# flag, ViewFactors falls back to ELMERSOLVER_STARTINFO, and every corner of
# the matrix below would get its radiator factors computed with the one method
# STARTINFO happens to name.  So name the sif for each corner explicitly.
#
# Each corner's view factors are kept so the corners can be compared against
# each other afterwards, not only against a norm.
MACRO(RUN_CORNER SIF)
  RUN_TOOL(${RADIATORS_BIN} ${SIF}-radiators ${SIF}.sif)
  RUN_TOOL(${VIEWFACTORS_BIN} ${SIF}-viewfactors ${SIF}.sif)
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

RUN_TOOL(${ELMERGRID_BIN} elmergrid 1 2 box_in_box.grd)

# ---------------------------------------------------------------------------
# The shaft cull is a pure speed optimisation: it replaces one tree traversal
# per ray with one shaft build per pair and must not change a single bit of
# the answer.  That is an exact invariant, so check it exactly rather than
# through a norm tolerance.  At eight rays over a six element shadow mesh the
# automatic rule leaves the cull off, so the two runs really do take different
# paths through the code.
# ---------------------------------------------------------------------------
RUN_CORNER(culloff)
RUN_CORNER(cullon)

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
