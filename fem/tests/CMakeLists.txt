project(ElmerTests Fortran C)
cmake_minimum_required(VERSION 2.8.9)

# ADD_EXECUTABLE(findnorm findnorm_cmake.c)

set(ELMERGRID_BIN "${CMAKE_BINARY_DIR}/elmergrid/src/ElmerGrid")
if(WITH_MPI)
  set(ELMERSOLVER_BIN "${CMAKE_BINARY_DIR}/fem/src/ElmerSolver_mpi")
else()
  set(ELMERSOLVER_BIN "${CMAKE_BINARY_DIR}/fem/src/ElmerSolver")
endif()
set(FINDNORM_BIN "${CMAKE_CURRENT_BINARY_DIR}/findnorm")
set(MESH2D_BIN "${CMAKE_BINARY_DIR}/meshgen2d/src/Mesh2D")

macro(SUBDIRLIST result curdir depth)
  set(glob_pattern "*")
  foreach(D RANGE 1 ${depth})
    file(
      GLOB children
      RELATIVE ${curdir}
      ${glob_pattern})
    foreach(child ${children})
      if(IS_DIRECTORY ${curdir}/${child})
        list(APPEND dirlist ${child})
      endif()
    endforeach()
    set(glob_pattern "${glob_pattern}/*")
  endforeach()
  set(${result} ${dirlist})
endmacro()

subdirlist(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR} 2)

set(_test_counter 0)
foreach(subdir ${SUBDIRS})
  if(NOT ${subdir} MATCHES "^contrib/" OR WITH_CONTRIB)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/CMakeLists.txt")
      math(EXPR _test_counter "${_test_counter}+1")
      add_subdirectory(${subdir})
      list(APPEND TEST_SUBDIRS ${subdir})
    endif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/CMakeLists.txt")
  endif()
endforeach()
message(STATUS "Found ${_test_counter} tests")

set(INSTALL_OLD_TESTS
    FALSE
    CACHE BOOL "(deprecated) Install old test system under build directory.")

foreach(_subdir ${TEST_SUBDIRS})
  file(
    GLOB files
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    "${CMAKE_CURRENT_SOURCE_DIR}/${_subdir}/*")
  foreach(_file ${files})
    if(NOT (${_file} MATCHES "CMakeLists.txt$"))
      if(NOT (${_file} MATCHES "runtest.cmake$"))
        list(APPEND ELMER_TEST_FILES ${_file})
      endif()
    endif()
  endforeach()
endforeach()

if(INSTALL_OLD_TESTS)
  install(
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.
    DESTINATION ${CMAKE_BINARY_DIR}/elmertests
    FILES_MATCHING
    REGEX
      ".+\\.(f90|F90|sif|mif|grd|eg|msh|good|dat|pos|result|in2d|best0)|HB|Makefile|TEST\\.RESULT|ELMERSOLVER_STARTINFO|mesh\\.*"
  )

  # Misc. directories for meshes etc.
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/diffuser_sa/diffuser/.
          DESTINATION ${CMAKE_BINARY_DIR}/elmertests/diffuser_sa/diffuser)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/diffuser_sst/diffuser/.
          DESTINATION ${CMAKE_BINARY_DIR}/elmertests/diffuser_sst/diffuser)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/diffuser_v2f/diffuser/.
          DESTINATION ${CMAKE_BINARY_DIR}/elmertests/diffuser_v2f/diffuser)

  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/findnorm.c
                ${CMAKE_CURRENT_SOURCE_DIR}/ElmerSolver_test_how-to.txt
          DESTINATION ${CMAKE_BINARY_DIR}/elmertests)
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/runtests
          DESTINATION ${CMAKE_BINARY_DIR}/elmertests)
endif(INSTALL_OLD_TESTS)

# ElmerTests_package target
set(ELMER_TESTS_PACKAGE_FILENAME
    "elmerfem-tests-${ELMER_FEM_REVISION}.tar.gz"
    CACHE STRING "Filename of test package")
mark_as_advanced(ELMER_TESTS_PACKAGE_FILENAME)

add_custom_target(ElmerTests "${CMAKE_COMMAND}" "-E" "make_directory"
                             "${CMAKE_BINARY_DIR}/elmerfem-tests")
add_custom_target(
  ElmerTests_package
  DEPENDS ElmerTests
  COMMAND "${CMAKE_COMMAND}" "-E" "tar" "cvfz" "${ELMER_TESTS_PACKAGE_FILENAME}"
          "elmerfem-tests"
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

# Copy the tests in a cross platform compliant manner
foreach(_file ${ELMER_TEST_FILES})
  if(NOT (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${_file}))
    add_custom_command(
      TARGET ElmerTests
      COMMAND
        "${CMAKE_COMMAND}" "-E" "copy" "${CMAKE_CURRENT_SOURCE_DIR}/${_file}"
        "${CMAKE_BINARY_DIR}/elmerfem-tests/${_file}")
  else()
    add_custom_command(
      TARGET ElmerTests
      COMMAND
        "${CMAKE_COMMAND}" "-E" "copy_directory"
        "${CMAKE_CURRENT_SOURCE_DIR}/${_file}"
        "${CMAKE_BINARY_DIR}/elmerfem-tests/${_file}")
  endif()
endforeach()
