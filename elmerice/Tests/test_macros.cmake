# ELMER_HOME and ELMER_LIB are what point the binaries at the build tree being
# tested rather than at whatever happens to be installed.  RUN_ELMER_TEST and
# EXECUTE_ELMER_SOLVER set them, but they only do so when they are called, and
# a test that drives ElmerGrid, ViewFactors or Radiators directly through
# EXECUTE_PROCESS spawns those before any macro has run.  Such a child then
# inherits an empty environment and falls back to the compiled in
# ELMER_SOLVER_HOME, i.e. to the install tree; on Windows getsolverhome()
# ignores that define and looks next to the running exe instead, which in a
# build tree holds no elements.def at all and the run dies before it writes
# anything.  Set them here, at include time, so every child of a runtest.cmake
# sees the same environment the solver gets.  Same values RUN_ELMER_TEST uses.
#
# BINARY_DIR is only defined when this is included from a runtest.cmake driven
# by ctest; at configure time, where the test CMakeLists.txt include it, there
# is nothing to point at and nothing to spawn.
#
# PATH matters just as much, and for a reason that only bites on Windows:
# there is no RPATH there, so a freshly linked ViewFactors.exe finds
# libelmersolver, fhutiter, matc, arpack and the compiler's own runtime DLLs
# through PATH alone.  A child spawned before any macro ran therefore does not
# start at all -- the loader kills it with STATUS_DLL_NOT_FOUND before main(),
# so it writes nothing, not even to the error log the test redirected.  The
# older radiation tests hid this: their ViewFactors call quietly did nothing
# and ElmerSolver, which does run with PATH set, recomputed the missing
# factors itself.  A test that reads ViewFactors.dat before the solver runs
# has nothing to fall back on and just reports a file that does not exist.
IF(DEFINED BINARY_DIR)
  IF(NOT DEFINED ENV{ELMER_HOME} OR "$ENV{ELMER_HOME}" STREQUAL "")
    SET(ENV{ELMER_HOME} "${BINARY_DIR}/fem/src")
  ENDIF()
  IF(NOT DEFINED ENV{ELMER_LIB} OR "$ENV{ELMER_LIB}" STREQUAL "")
    SET(ENV{ELMER_LIB} "${BINARY_DIR}/fem/src/modules")
  ENDIF()
  IF(NOT DEFINED ENV{ELMER_MODULES_PATH} OR "$ENV{ELMER_MODULES_PATH}" STREQUAL "")
    SET(ENV{ELMER_MODULES_PATH} "${BINARY_DIR}/elmerice/Solvers:${BINARY_DIR}/elmerice/Solvers/GridDataReader:${BINARY_DIR}/elmerice/Solvers/ScatteredDataInterpolator:${BINARY_DIR}/elmerice/Solvers/MeshAdaptation_2D:${BINARY_DIR}/elmerice/UserFunctions")
  ENDIF()
  IF(NOT(WIN32))
    SET(ENV{PATH} "${BINARY_DIR}/meshgen2d/src/:${BINARY_DIR}/elmerice/Solvers:${BINARY_DIR}/elmerice/Utils:${BINARY_DIR}/fem/src:$ENV{PATH}")
  ELSE()
    GET_FILENAME_COMPONENT(COMPILER_DIRECTORY ${CMAKE_Fortran_COMPILER} PATH)
    SET(ENV{PATH} "$ENV{ELMER_HOME};$ENV{ELMER_LIB};${BINARY_DIR}/meshgen2d/src/;${BINARY_DIR}/fhutiter/src;${BINARY_DIR}/matc/src;${BINARY_DIR}/mathlibs/src/arpack;${BINARY_DIR}/mathlibs/src/parpack;${COMPILER_DIRECTORY};${BINARY_DIR}/elmerice/Solvers;${BINARY_DIR}/elmerice/Utils;${BINARY_DIR}/fem/src;$ENV{PATH}")
  ENDIF()
ENDIF()

MACRO(ADD_ELMERICE_LABEL test_name label_string)
  SET_PROPERTY(TEST ${test_name} APPEND PROPERTY LABELS ${label_string})
ENDMACRO()

MACRO(ADD_ELMERICE_TEST test_name)
  ADD_TEST(NAME ${test_name}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${CMAKE_COMMAND}
      -DELMERGRID_BIN=${ELMERGRID_BIN}
      -DELMERSOLVER_BIN=${ELMERSOLVER_BIN}
      -DTEST_SOURCE=${CMAKE_CURRENT_SOURCE_DIR}
      -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
      -DBINARY_DIR=${CMAKE_BINARY_DIR}
      -DELMERSOLVER_HOME=${ELMER_SOLVER_HOME}
      -DSHLEXT=${SHL_EXTENSION}
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
      -DMPIEXEC=${MPIEXEC}
      -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
      -DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}
      -DMPIEXEC_POSTFLAGS=${MPIEXEC_POSTFLAGS}
      -DWITH_MPI=${WITH_MPI}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/runTest.cmake)
  SET_TESTS_PROPERTIES(${test_name} PROPERTIES LABELS "elmerice")
  IF(WITH_MPI AND DEFINED NPROCS)
    # Tell ctest how many processors this test requires.
    SET_TESTS_PROPERTIES(${test_name} PROPERTIES
      PROCESSORS ${NPROCS})
  ENDIF()
ENDMACRO()

MACRO(ADD_ELMERICETEST_MODULE test_name module_name file_name)
  IF(APPLE)
    SET(CMAKE_SHARED_MODULE_SUFFIX ".dylib")
  ENDIF(APPLE)
  SET(ELMERICETEST_CMAKE_NAME "${test_name}_${module_name}")
  ADD_LIBRARY(${ELMERICETEST_CMAKE_NAME} MODULE ${file_name})
  SET_TARGET_PROPERTIES(${ELMERICETEST_CMAKE_NAME}
    PROPERTIES PREFIX "")
  SET_TARGET_PROPERTIES(${ELMERICETEST_CMAKE_NAME}
    PROPERTIES OUTPUT_NAME ${module_name} LINKER_LANGUAGE Fortran)
  TARGET_LINK_LIBRARIES(${ELMERICETEST_CMAKE_NAME} elmersolver)
  TARGET_LINK_LIBRARIES(${ELMERICETEST_CMAKE_NAME} ElmerIceUtils)
#  IF(WITH_MPI)
#    ADD_DEPENDENCIES(${ELMERICETEST_CMAKE_NAME} elmersolver)
#  ELSE()
#    ADD_DEPENDENCIES(${ELMERICETEST_CMAKE_NAME} elmersolver)
#  ENDIF()
  UNSET(ELMERICETEST_CMAKE_NAME)
ENDMACRO()

MACRO(RUN_ELMERICE_TEST)
  MESSAGE(STATUS "BINARY_DIR = ${BINARY_DIR}")

  # Query number of physical CPU cores of the host
  cmake_host_system_information(RESULT PHYS_CPU QUERY NUMBER_OF_PHYSICAL_CORES)

  # Optional arguments like WITH_MPI
  SET(LIST_VAR "${ARGN}")
  IF(LIST_VAR STREQUAL "")
    FILE(REMOVE "TEST.PASSED")

    # Limit number of OpenMP threads to a sensible value
    IF(NOT DEFINED ENV{OMP_NUM_THREADS})
      # set limit
      SET(ENV{OMP_NUM_THREADS} ${PHYS_CPU})
    ENDIF()

    IF(WITH_MPI)
      EXECUTE_PROCESS(COMMAND "${MPIEXEC}" ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ${ELMERSOLVER_BIN} ${MPIEXEC_POSTFLAGS}
        OUTPUT_FILE "test-stdout.log"
        ERROR_FILE "test-stderr.log")
    ELSE()
      EXECUTE_PROCESS(COMMAND ${ELMERSOLVER_BIN}
        OUTPUT_FILE "test-stdout.log"
        ERROR_FILE "test-stderr.log")
    ENDIF()
  ELSEIF("${LIST_VAR}" STREQUAL "WITH_MPI" AND WITH_MPI)
    # Macro has been called with WITH_MPI argument and MPI is enabled
    SET(N "${NPROCS}")
    IF("${N}" STREQUAL "")
      MESSAGE(FATAL_ERROR "Test failed:variable <NPROC> not defined. Set <NPROC> in runTest.cmake")
    ENDIF()

    FILE(REMOVE "TEST.PASSED_${N}")

    IF(NOT DEFINED ENV{OMP_NUM_THREADS})
      # Limit number of OpenMP threads to a sensible value
      # Divide by ${NPROCS} and truncate to the nearest lower integer
      MATH(EXPR n_openmp_threads "${PHYS_CPU} / ${NPROCS}")
      IF(${n_openmp_threads} LESS 1)
        # minimum of 1
        SET(ENV{OMP_NUM_THREADS} 1)
      ELSE()
        # set limit
        SET(ENV{OMP_NUM_THREADS} ${n_openmp_threads})
      ENDIF()
    ENDIF()

    EXECUTE_PROCESS(COMMAND "${MPIEXEC}" ${MPIEXEC_NUMPROC_FLAG} ${N} ${MPIEXEC_PREFLAGS} ${ELMERSOLVER_BIN} ${MPIEXEC_POSTFLAGS}
      OUTPUT_FILE "test-stdout.log"
      ERROR_FILE "test-stderr.log")
  ENDIF()

  IF(NPROCS GREATER "1")
    SET(_passed_file "TEST.PASSED_${NPROCS}")
  ELSE()
    SET(_passed_file "TEST.PASSED")
  ENDIF()

  # Same reasoning as the fem copy of this macro: a missing result file means
  # the solver produced nothing, and FILE(READ) would report only that a path
  # does not exist. This macro sends its output to files rather than to a
  # variable, so the norms are recovered from the log.
  IF(NOT EXISTS "${_passed_file}")
    MESSAGE(FATAL_ERROR
      "the solver produced no ${_passed_file} at all -- it did not run to "
      "completion.\n  See test-stdout.log and test-stderr.log in "
      "${CMAKE_CURRENT_BINARY_DIR}.")
  ENDIF()

  FILE(READ "${_passed_file}" RES)
  IF(NOT RES EQUAL "1")
    SET(_cmp "")
    IF(EXISTS "test-stdout.log")
      FILE(STRINGS "test-stdout.log" _cmp_lines
        REGEX "CompareToReferenceSolution")
      IF(_cmp_lines)
        STRING(REPLACE ";" "\n  " _cmp "${_cmp_lines}")
        SET(_cmp "\n  ${_cmp}")
      ENDIF()
    ENDIF()
    MESSAGE(FATAL_ERROR
      "the solver ran but its result did not match the reference${_cmp}\n"
      "  See test-stdout.log and test-stderr.log in "
      "${CMAKE_CURRENT_BINARY_DIR}.")
  ENDIF()
ENDMACRO()

MACRO(EXECUTE_ELMER_SOLVER SIFNAME)
  EXECUTE_PROCESS(COMMAND ${ELMERSOLVER_BIN} ${SIFNAME} OUTPUT_FILE "${SIFNAME}-stdout.log"
    ERROR_FILE "${SIFNAME}-stderr.log")
ENDMACRO()
