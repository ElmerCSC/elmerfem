include(CMakeParseArguments)
macro(ADD_ELMER_LABEL test_name label_string)
  set_property(
    TEST ${test_name}
    APPEND
    PROPERTY LABELS ${label_string})
endmacro()

macro(ADD_ELMER_TEST TestName)
  # Parse optional named arguments, NPROCS and LABELS, which can be lists
  cmake_parse_arguments(_parsedArgs "" "" "NPROCS;LABELS" "${ARGN}")

  if(_parsedArgs_NPROCS)
    # List of task counts was given so this is a parallel test case
    foreach(n ${_parsedArgs_NPROCS})
      if(WITH_MPI AND ${n} GREATER 1)
        # Check the task bounds and add only compatible tests
        if(${n} GREATER ${MPI_TEST_MAXPROC} OR ${n} LESS ${MPI_TEST_MINPROC})
          message(STATUS "Skipping test ${TestName} with ${n} procs")
        else()
          list(APPEND tests_list "${TestName}_np${n}")
          list(APPEND label_list "parallel")
          list(APPEND tasks_list "${n}")
        endif()
      else()
        # If there is a single task version in the task list, add it as a test
        # case also for non-MPI builds
        if(${n} EQUAL 1)
          set(tests_list "${TestName}")
          set(label_list "serial")
          set(tasks_list 1)
        endif()
      endif()
    endforeach()
  else()
    # No NPROCS argument, serial test
    set(tests_list "${TestName}")
    set(label_list "serial")
    set(tasks_list 1)
  endif()

  # Loop over the two lists, which is cumbersome in CMake
  list(LENGTH tests_list nt)
  math(EXPR ntests "${nt} - 1")
  foreach(n RANGE ${ntests})
    list(GET tests_list ${n} _this_test_name)
    list(GET tasks_list ${n} _this_test_tasks)
    list(GET label_list ${n} _this_test_label)
    if(_this_test_name)
      add_test(
        NAME ${_this_test_name}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMAND
          ${CMAKE_COMMAND} -DCMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}
          -DELMERGRID_BIN=${ELMERGRID_BIN} -DELMERSOLVER_BIN=${ELMERSOLVER_BIN}
          -DFINDNORM_BIN=${FINDNORM_BIN} -DMESH2D_BIN=${MESH2D_BIN}
          -DTEST_SOURCE=${CMAKE_CURRENT_SOURCE_DIR}
          -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
          -DBINARY_DIR=${CMAKE_BINARY_DIR}
          -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
          -DMPIEXEC=${MPIEXEC} -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
          -DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}
          -DMPIEXEC_POSTFLAGS=${MPIEXEC_POSTFLAGS} -DWITH_MPI=${WITH_MPI}
          -DMPIEXEC_NTASKS=${_this_test_tasks} -P
          ${CMAKE_CURRENT_SOURCE_DIR}/runtest.cmake)
      set_property(TEST ${_this_test_name} APPEND PROPERTY LABELS
                                                           ${_this_test_label})
      # If LABELS argument was given iterate through the given labels and add
      # them to this test
      if(_parsedArgs_LABELS)
        foreach(lbl ${_parsedArgs_LABELS})
          set_property(TEST ${_this_test_name} APPEND PROPERTY LABELS ${lbl})
        endforeach()
      endif()
    endif(_this_test_name)
  endforeach()
endmacro(ADD_ELMER_TEST)

macro(ADD_ELMERTEST_MODULE test_name module_name file_name)
  if(APPLE)
    set(CMAKE_SHARED_MODULE_SUFFIX ".dylib")
  endif(APPLE)
  set(ELMERTEST_CMAKE_NAME "${test_name}_${module_name}")
  add_library(${ELMERTEST_CMAKE_NAME} MODULE ${file_name})
  set_target_properties(${ELMERTEST_CMAKE_NAME} PROPERTIES PREFIX "")
  target_link_libraries(${ELMERTEST_CMAKE_NAME} elmersolver)
  set_target_properties(
    ${ELMERTEST_CMAKE_NAME} PROPERTIES OUTPUT_NAME ${module_name}
                                       LINKER_LANGUAGE Fortran)
  target_link_libraries(${ELMERTEST_CMAKE_NAME} elmersolver)
  add_dependencies(${ELMERTEST_CMAKE_NAME} elmersolver Solver_TGT ElmerGrid)
  unset(ELMERTEST_CMAKE_NAME)
endmacro()

macro(RUN_ELMER_TEST)
  cmake_parse_arguments(_parsedArgs "" "" "ELMER_LIB" "${ARGN}")
  message(STATUS "BINARY_DIR = ${BINARY_DIR}")
  set(ENV{ELMER_HOME} "${BINARY_DIR}/fem/src")
  set(ENV{ELMER_LIB} "${BINARY_DIR}/fem/src/modules")
  if(NOT _parsedArgs_ELMER_LIB STREQUAL "")
    message(STATUS "Extra library directories ${_parsedArgs_ELMER_LIB}")
    foreach(_extra_lib ${_parsedArgs_ELMER_LIB})
      if(NOT (WIN32))
        set(ENV{ELMER_LIB}
            "$ENV{ELMER_LIB}:${BINARY_DIR}/fem/src/modules/${_extra_lib}")
      else()
        set(ENV{ELMER_LIB}
            "$ENV{ELMER_LIB};${BINARY_DIR}/fem/src/modules/${_extra_lib}")
      endif()
    endforeach()
  endif()

  if(NOT (WIN32))
    set(ENV{PATH}
        "$ENV{PATH}:${BINARY_DIR}/meshgen2d/src/:${BINARY_DIR}/fem/src")
  endif(NOT (WIN32))

  # Clean up old result files
  if(${MPIEXEC_NTASKS} GREATER 1)
    file(REMOVE TEST.PASSED_${MPIEXEC_NTASKS})
  else()
    file(REMOVE TEST.PASSED)
  endif()

  if(WIN32)
    set(ENV{PATH}
        "$ENV{PATH};${BINARY_DIR}/meshgen2d/src/;${BINARY_DIR}/fem/src")
    get_filename_component(COMPILER_DIRECTORY ${CMAKE_Fortran_COMPILER} PATH)
    set(ENV{PATH}
        "$ENV{PATH};${COMPILER_DIRECTORY};$ENV{ELMER_HOME};$ENV{ELMER_LIB};${BINARY_DIR}/fhutiter/src;${BINARY_DIR}/matc/src;${BINARY_DIR}/mathlibs/src/arpack"
    )
  endif(WIN32)

  if(WITH_MPI)
    execute_process(
      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_NTASKS}
              ${MPIEXEC_PREFLAGS} ${ELMERSOLVER_BIN} ${MPIEXEC_POSTFLAGS}
      OUTPUT_VARIABLE TEST_STDOUT_VARIABLE
      ERROR_VARIABLE TEST_ERROR_VARIABLE)
    file(WRITE "test-stdout_${MPIEXEC_NTASKS}.log" "${TEST_STDOUT_VARIABLE}")
    file(WRITE "test-stderr_${MPIEXEC_NTASKS}.log" "${TEST_STDERR_VARIABLE}")
  else()
    execute_process(
      COMMAND ${ELMERSOLVER_BIN}
      OUTPUT_VARIABLE TEST_STDOUT_VARIABLE
      ERROR_VARIABLE TEST_ERROR_VARIABLE)
    file(WRITE "test-stdout.log" "${TEST_STDOUT_VARIABLE}")
    file(WRITE "test-stderr.log" "${TEST_STDERR_VARIABLE}")
  endif()
  if($ENV{CTEST_OUTPUT_ON_FAILURE})
    message(STATUS "stdout:")
    message("${TEST_STDOUT_VARIABLE}")
    message(STATUS "stderr:")
    message("${TEST_STDERR_VARIABLE}")
  endif()

  # Check the result file (with suffix is more than single task)
  if(${MPIEXEC_NTASKS} GREATER 1)
    file(READ "TEST.PASSED_${MPIEXEC_NTASKS}" RES)
  else()
    file(READ "TEST.PASSED" RES)
  endif()
  if(NOT RES EQUAL "1")
    message(FATAL_ERROR "Test failed")
  endif()
endmacro()

macro(EXECUTE_ELMER_SOLVER SIFNAME)
  set(ENV{ELMER_HOME} "${BINARY_DIR}/fem/src")
  set(ENV{ELMER_LIB} "${BINARY_DIR}/fem/src/modules")

  if(NOT (WIN32))
    set(ENV{PATH}
        "$ENV{PATH}:${BINARY_DIR}/meshgen2d/src/:${BINARY_DIR}/fem/src")
  endif(NOT (WIN32))

  if(WIN32)
    set(ENV{PATH}
        "$ENV{PATH};${BINARY_DIR}/meshgen2d/src/;${BINARY_DIR}/fem/src")
    get_filename_component(COMPILER_DIRECTORY ${CMAKE_Fortran_COMPILER} PATH)
    set(ENV{PATH}
        "$ENV{PATH};${COMPILER_DIRECTORY};$ENV{ELMER_HOME};$ENV{ELMER_LIB};${BINARY_DIR}/fhutiter/src;${BINARY_DIR}/matc/src;${BINARY_DIR}/mathlibs/src/arpack"
    )
  endif(WIN32)

  execute_process(
    COMMAND ${ELMERSOLVER_BIN} ${SIFNAME}
    OUTPUT_FILE "${SIFNAME}-stdout.log"
    ERROR_FILE "${SIFNAME}-stderr.log")
endmacro()
