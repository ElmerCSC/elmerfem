MACRO(ADD_ELMER_LABEL test_name label_string)
  SET_PROPERTY(TEST ${test_name} APPEND PROPERTY LABELS ${label_string})
ENDMACRO()


MACRO(ADD_ELMER_TEST TestName)
  # Parse optional named arguments, NPROCS and LABELS, which can be lists
  CMAKE_PARSE_ARGUMENTS(_parsedArgs "" "" "NPROCS;LABELS" "${ARGN}")

  IF(_parsedArgs_LABELS)
    # List of task counts was given so this is a parallel test case
    FOREACH(n ${_parsedArgs_NPROCS})
      IF(WITH_MPI)
        # Check the task bounds and add only compatible tests
        IF(${n} GREATER ${MPI_TEST_MAXPROC} OR ${n} LESS ${MPI_TEST_MINPROC})
          MESSAGE(STATUS "Skipping test ${TestName} with ${n} procs")
        ELSE()
          LIST(APPEND tests_list "${TestName}_np${n}")
          LIST(APPEND tasks_list "${n}")
        ENDIF()
      ELSE(WITH_MPI)
        # If there is a single task version in the task list, add
        # it as a test case also for non-MPI builds
        IF(${n} EQUAL 1)
          SET(tests_list "${TestName}")
          SET(tasks_list 1)
        ENDIF()
      ENDIF(WITH_MPI)
    ENDFOREACH()
  ELSE()
    # Serial or purely single task test
    SET(tests_list "${TestName}")
    SET(tasks_list 1)
  ENDIF()

  # Loop over the two lists, which is cumbersome in CMake
  LIST(LENGTH tests_list nt)
  MATH(EXPR ntests "${nt} - 1")
  FOREACH(n RANGE ${ntests})
    LIST(GET tests_list ${n} this_test_name)
    LIST(GET tasks_list ${n} this_test_tasks)
    ADD_TEST(NAME ${this_test_name}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMAND ${CMAKE_COMMAND}
      -DELMERGRID_BIN=${ELMERGRID_BIN}
      -DELMERSOLVER_BIN=${ELMERSOLVER_BIN}
      -DFINDNORM_BIN=${FINDNORM_BIN}
      -DMESH2D_BIN=${MESH2D_BIN}
      -DTEST_SOURCE=${CMAKE_CURRENT_SOURCE_DIR}
      -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
      -DBINARY_DIR=${CMAKE_BINARY_DIR}
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
      -DMPIEXEC=${MPIEXEC}
      -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
      -DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}
      -DMPIEXEC_POSTFLAGS=${MPIEXEC_POSTFLAGS}
      -DWITH_MPI=${WITH_MPI}
      -DMPIEXEC_NTASKS=${this_test_tasks}
      -P ${CMAKE_SOURCE_DIR}/fem/tests/test_macros.cmake
      -P ${CMAKE_CURRENT_SOURCE_DIR}/runtest.cmake)
    IF(_parsedArgs_LABELS)
      FOREACH(lbl ${_parsedArgs_LABELS})
        SET_PROPERTY(TEST ${this_test_name} APPEND PROPERTY LABELS ${lbl})
      ENDFOREACH()
    ENDIF()
  ENDFOREACH()
ENDMACRO(ADD_ELMER_TEST)


MACRO(ADD_ELMERTEST_MODULE test_name module_name file_name)
  IF(APPLE)
    SET(CMAKE_SHARED_MODULE_SUFFIX ".dylib")
  ENDIF(APPLE)
  SET(ELMERTEST_CMAKE_NAME "${test_name}_${module_name}")
  ADD_LIBRARY(${ELMERTEST_CMAKE_NAME} MODULE ${file_name})
  SET_TARGET_PROPERTIES(${ELMERTEST_CMAKE_NAME}
    PROPERTIES PREFIX "")
  TARGET_LINK_LIBRARIES(${ELMERTEST_CMAKE_NAME}
    elmersolver)
  SET_TARGET_PROPERTIES(${ELMERTEST_CMAKE_NAME}
    PROPERTIES OUTPUT_NAME ${module_name} LINKER_LANGUAGE Fortran)
  TARGET_LINK_LIBRARIES(${ELMERTEST_CMAKE_NAME} elmersolver)
  ADD_DEPENDENCIES(${ELMERTEST_CMAKE_NAME} 
    elmersolver Solver_TGT ElmerGrid)
  UNSET(ELMERTEST_CMAKE_NAME)
ENDMACRO()


MACRO(RUN_ELMER_TEST)
  MESSAGE(STATUS "BINARY_DIR = ${BINARY_DIR}")
  SET(ENV{ELMER_HOME} "${BINARY_DIR}/fem/src")
  SET(ENV{ELMER_LIB} "${BINARY_DIR}/fem/src/modules")

  IF(NOT(WIN32))
    SET(ENV{PATH} "$ENV{PATH}:${BINARY_DIR}/meshgen2d/src/:${BINARY_DIR}/fem/src")
  ENDIF(NOT(WIN32))

  # Clean up old result files
  IF(${MPIEXEC_NTASKS} GREATER 1)
    FILE(REMOVE TEST.PASSED_${MPIEXEC_NTASKS})
  ELSE()
    FILE(REMOVE TEST.PASSED)
  ENDIF()

  IF(WIN32)
    SET(ENV{PATH} "$ENV{PATH};${BINARY_DIR}/meshgen2d/src/;${BINARY_DIR}/fem/src")
    GET_FILENAME_COMPONENT(COMPILER_DIRECTORY ${CMAKE_Fortran_COMPILER} PATH)
    SET(ENV{PATH} "$ENV{PATH};${COMPILER_DIRECTORY};$ENV{ELMER_HOME};$ENV{ELMER_LIB};${BINARY_DIR}/fhutiter/src;${BINARY_DIR}/matc/src;${BINARY_DIR}/mathlibs/src/arpack")
  ENDIF(WIN32)

  IF(WITH_MPI)
    EXECUTE_PROCESS(COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_NTASKS} ${MPIEXEC_PREFLAGS} ${ELMERSOLVER_BIN} ${MPIEXEC_POSTFLAGS}
      OUTPUT_FILE "test-stdout_${MPIEXEC_NTASKS}.log"
      ERROR_FILE "test-stderr_${MPIEXEC_NTASKS}.log")
  ELSE()
    EXECUTE_PROCESS(COMMAND ${ELMERSOLVER_BIN}
      OUTPUT_FILE "test-stdout.log"
      ERROR_FILE "test-stderr.log")
  ENDIF()

  # Check the result file (with suffix is more than single task)
  IF(${MPIEXEC_NTASKS} GREATER 1)
    FILE(READ "TEST.PASSED_${MPIEXEC_NTASKS}" RES)
  ELSE()
    FILE(READ "TEST.PASSED" RES)
  ENDIF()
  IF(NOT RES EQUAL "1")
    MESSAGE(FATAL_ERROR "Test failed")
  ENDIF()
ENDMACRO()


MACRO(EXECUTE_ELMER_SOLVER SIFNAME)
  SET(ENV{ELMER_HOME} "${BINARY_DIR}/fem/src")
  SET(ENV{ELMER_LIB} "${BINARY_DIR}/fem/src/modules")

  IF(NOT(WIN32))
    SET(ENV{PATH} "$ENV{PATH}:${BINARY_DIR}/meshgen2d/src/:${BINARY_DIR}/fem/src")
  ENDIF(NOT(WIN32))

  IF(WIN32)
    SET(ENV{PATH} "$ENV{PATH};${BINARY_DIR}/meshgen2d/src/;${BINARY_DIR}/fem/src")
    GET_FILENAME_COMPONENT(COMPILER_DIRECTORY ${CMAKE_Fortran_COMPILER} PATH)
    SET(ENV{PATH} "${COMPILER_DIRECTORY};$ENV{ELMER_HOME};$ENV{ELMER_LIB};${BINARY_DIR}/fhutiter/src;${BINARY_DIR}/matc/src;${BINARY_DIR}/mathlibs/src/arpack")
  ENDIF(WIN32)

  EXECUTE_PROCESS(COMMAND ${ELMERSOLVER_BIN} ${SIFNAME}
    OUTPUT_FILE "${SIFNAME}-stdout.log"
    ERROR_FILE "${SIFNAME}-stderr.log")
ENDMACRO()
