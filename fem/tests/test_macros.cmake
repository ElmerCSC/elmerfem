MACRO(ADD_ELMER_LABEL test_name label_string)
  SET_PROPERTY(TEST ${test_name} APPEND PROPERTY LABELS ${label_string})
ENDMACRO()


MACRO(ADD_ELMER_TEST test_name)
  # Check the number of input arguments
  IF(${ARGC} GREATER 1) # This is a parallel test
    # Construct a list with test names and number of tasks,
    # note the suffix
    FOREACH(n ${ARGN})
      LIST(APPEND test_list "${test_name}_${n}")
      LIST(APPEND tasks_list "${n}")
    ENDFOREACH()
  ELSE()
    # Serial or single task test
    SET(test_list "${test_name}")
    SET(task_list 1)
  ENDIF()

  # Loop over the two lists, which is cumbersome in CMake
  LIST(LENGTH test_list nt)
  MATH(EXPR ntests "${nt} - 1")
  FOREACH(n RANGE ${ntests})
    LIST(GET test_list ${n} this_test_name)
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
  ENDFOREACH()
ENDMACRO()


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

  FILE(REMOVE TEST.PASSED)

  IF(WIN32)
    SET(ENV{PATH} "$ENV{PATH};${BINARY_DIR}/meshgen2d/src/;${BINARY_DIR}/fem/src")
    GET_FILENAME_COMPONENT(COMPILER_DIRECTORY ${CMAKE_Fortran_COMPILER} PATH)
    SET(ENV{PATH} "$ENV{PATH};${COMPILER_DIRECTORY};$ENV{ELMER_HOME};$ENV{ELMER_LIB};${BINARY_DIR}/fhutiter/src;${BINARY_DIR}/matc/src;${BINARY_DIR}/mathlibs/src/arpack")
  ENDIF(WIN32)

  IF(WITH_MPI)
    EXECUTE_PROCESS(COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_NTASKS} ${MPIEXEC_PREFLAGS} ${ELMERSOLVER_BIN} ${MPIEXEC_POSTFLAGS}
      OUTPUT_FILE "test-stdout_${MPIEXEC_NTASKS}.log"
      ERROR_FILE "test-stderr_${MPIEXEC_NTASKS}.log"
      OUTPUT_VARIABLE TESTOUTPUT)
  ELSE()      
    EXECUTE_PROCESS(COMMAND ${ELMERSOLVER_BIN}
      OUTPUT_FILE "test-stdout.log"
      ERROR_FILE "test-stderr.log"
      OUTPUT_VARIABLE TESTOUTPUT)
  ENDIF()

  MESSAGE(STATUS "testoutput.........: ${TESTOUTPUT}")

  FILE(READ "TEST.PASSED" RES)
  IF (NOT RES EQUAL "1")
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

