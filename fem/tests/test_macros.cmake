

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
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
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
  set(ENV{PATH} "$ENV{PATH};${BINARY_DIR}/meshgen2d/src;c:/elmer_stuff/mingw64/bin;$ENV{ELMER_HOME};$ENV{ELMER_LIB};${BINARY_DIR}/fhutiter/src;${BINARY_DIR}/matc/src;${BINARY_DIR}/mathlibs/src/arpack")

  #MESSAGE(STATUS "\$Env:path=\"$ENV{PATH}\"")
  execute_process(COMMAND ${ELMERSOLVER_BIN} OUTPUT_FILE "elmersolver-stdout.log"
    ERROR_FILE "elmersolver-stderr.log" OUTPUT_VARIABLE TESTOUTPUT)
  execute_process(COMMAND ${FINDNORM_BIN} ${CMAKE_CURRENT_BINARY_DIR}/elmersolver-stdout.log
    OUTPUT_FILE "findnorm-stdout.log" ERROR_FILE "findnorm-stderr.log" RESULT_VARIABLE RES OUTPUT_VARIABLE FINDNORM_OUTPUT)

  MESSAGE(STATUS "findnorm res.......: ${RES}")
  MESSAGE(STATUS "testoutput.........: ${TESTOUTPUT}")
  MESSAGE(STATUS "findnorm output....: ${FINDNORM_OUTPUT}")
  if (RES)
    message(FATAL_ERROR "Test failed")
  endif()
endmacro()

macro(RUN_ELMER_TEST_NEW)
  set(ENV{ELMER_HOME} "${BINARY_DIR}/fem/src")
  set(ENV{ELMER_LIB} "${BINARY_DIR}/fem/src/modules")
  set(ENV{PATH} "$ENV{PATH};${BINARY_DIR}/meshgen2d/src/;${BINARY_DIR}/fem/src")
  IF(WIN32)
    GET_FILENAME_COMPONENT(COMPILER_DIRECTORY ${CMAKE_Fortran_COMPILER} PATH)
    set(ENV{PATH} "$ENV{PATH};${COMPILER_DIRECTORY};$ENV{ELMER_HOME};$ENV{ELMER_LIB};${BINARY_DIR}/fhutiter/src;${BINARY_DIR}/matc/src;${BINARY_DIR}/mathlibs/src/arpack")
  ENDIF(WIN32)
  execute_process(COMMAND ${ELMERSOLVER_BIN} OUTPUT_FILE "test-stdout.log"
    ERROR_FILE "test-stderr.log" OUTPUT_VARIABLE TESTOUTPUT)
  MESSAGE(STATUS "testoutput.........: ${TESTOUTPUT}")
  file(READ "TEST.PASSED" RES)
  if (NOT RES EQUAL "1")
    message(FATAL_ERROR "Test failed")
  endif()
endmacro()

macro(EXECUTE_ELMER_SOLVER SIFNAME)
  set(ENV{ELMER_HOME} "${BINARY_DIR}/fem/src")
  set(ENV{ELMER_LIB} "${BINARY_DIR}/fem/src/modules")
  set(ENV{PATH} "$ENV{PATH}:${BINARY_DIR}/meshgen2d/src/:${BINARY_DIR}/fem/src")
  IF(WIN32)
    GET_FILENAME_COMPONENT(COMPILER_DIRECTORY ${CMAKE_Fortran_COMPILER} PATH)
    set(ENV{PATH} "${COMPILER_DIRECTORY};$ENV{ELMER_HOME};$ENV{ELMER_LIB};${BINARY_DIR}/fhutiter/src;${BINARY_DIR}/matc/src;${BINARY_DIR}/mathlibs/src/arpack")
  ENDIF(WIN32)
  execute_process(COMMAND ${ELMERSOLVER_BIN} ${SIFNAME} OUTPUT_FILE "${SIFNAME}-stdout.log"
    ERROR_FILE "${SIFNAME}-stderr.log")
endmacro()
