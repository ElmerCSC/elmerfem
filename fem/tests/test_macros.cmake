

macro(ADD_ELMER_TEST test_name)
  add_test(NAME ${test_name}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${CMAKE_COMMAND}
      -DELMERGRID_BIN=${ELMERGRID_BIN}
      -DELMERSOLVER_BIN=${ELMERSOLVER_BIN}
      -DFINDNORM_BIN=${FINDNORM_BIN}
      -DMESH2D_BIN=${MESH2D_BIN}
      -DTEST_SOURCE=${CMAKE_CURRENT_SOURCE_DIR}
      -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
      -DBINARY_DIR=${CMAKE_BINARY_DIR}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/runtest.cmake)
endmacro()

macro(ADD_ELMERTEST_MODULE test_name module_name file_name)
  SET(ELMERTEST_CMAKE_NAME "${test_name}_${module_name}")
  ADD_LIBRARY(${ELMERTEST_CMAKE_NAME} MODULE ${file_name})
  SET_TARGET_PROPERTIES(${ELMERTEST_CMAKE_NAME}
    PROPERTIES PREFIX "")
  TARGET_LINK_LIBRARIES(${ELMERTEST_CMAKE_NAME}
    elmersolver)
  SET_TARGET_PROPERTIES(${ELMERTEST_CMAKE_NAME}
    PROPERTIES OUTPUT_NAME ${module_name})
  TARGET_LINK_LIBRARIES(${ELMERTEST_CMAKE_NAME} elmersolver)
  IF(WITH_MPI)
    ADD_DEPENDENCIES(${ELMERTEST_CMAKE_NAME} 
      elmersolver ElmerSolver_mpi ElmerGrid)
  ELSE()
    ADD_DEPENDENCIES(${ELMERTEST_CMAKE_NAME} 
      elmersolver ElmerSolver ElmerGrid)
  ENDIF()
  UNSET(ELMERTEST_CMAKE_NAME)
endmacro()

macro(RUN_ELMER_TEST)
  set(ENV{ELMER_HOME} "${BINARY_DIR}/fem/src")
  set(ENV{ELMER_LIB} "${BINARY_DIR}/fem/src/modules")
  set(ENV{PATH} "$ENV{PATH}:${BINARY_DIR}/meshgen2d/src/")
  execute_process(COMMAND ${ELMERSOLVER_BIN} OUTPUT_FILE "test.log"
    ERROR_FILE "test.log")
  execute_process(COMMAND ${FINDNORM_BIN} ${CMAKE_CURRENT_BINARY_DIR}/test.log
    OUTPUT_FILE norm.log ERROR_FILE norm.log RESULT_VARIABLE RES)

  if (RES)
    message(FATAL_ERROR "Test failed")
  endif()
endmacro()

macro(RUN_ELMER_TEST_NEW)
  set(ENV{ELMER_HOME} "${BINARY_DIR}/fem/src")
  set(ENV{ELMER_LIB} "${BINARY_DIR}/fem/src/modules")
  set(ENV{PATH} "$ENV{PATH}:${BINARY_DIR}/meshgen2d/src/")
  execute_process(COMMAND ${ELMERSOLVER_BIN} OUTPUT_FILE "test-stdout.log"
    ERROR_FILE "test-stderr.log")
  file(READ "TEST.PASSED" RES)

  if (NOT RES EQUAL "1")
    message(FATAL_ERROR "Test failed")
  endif()
endmacro()

