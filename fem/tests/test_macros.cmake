

macro(ADD_ELMER_TEST test_name)
  add_test(NAME ${test_name}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${CMAKE_COMMAND}
      -DELMERGRID_BIN=${ELMERGRID_BIN}
      -DELMERSOLVER_BIN=${ELMERSOLVER_BIN}
      -DFINDNORM_BIN=${FINDNORM_BIN}
      -DTEST_SOURCE=${CMAKE_CURRENT_SOURCE_DIR}
      -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/runtest.cmake)
endmacro()


macro(RUN_ELMER_TEST)
  set(ENV{ELMER_LIB} "${PROJECT_SOURCE_DIR}/fem/src")
  set(ENV{ELMER_HOME} "${PROJECT_SOURCE_DIR}/fem/src")
  execute_process(COMMAND ${ELMERSOLVER_BIN} OUTPUT_FILE "test.log"
    ERROR_FILE "test.log")
  execute_process(COMMAND ${FINDNORM_BIN} ${CMAKE_CURRENT_BINARY_DIR}/test.log
    OUTPUT_FILE norm.log ERROR_FILE norm.log RESULT_VARIABLE RES)

  if (RES)
    message(FATAL_ERROR "Test failed")
  endif()
endmacro()
