PROJECT(ElmerTests Fortran C)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.9)

#ADD_EXECUTABLE(findnorm findnorm_cmake.c)

SET(ELMERGRID_BIN "${CMAKE_BINARY_DIR}/elmergrid/src/ElmerGrid")
IF(WITH_MPI)
  SET(ELMERSOLVER_BIN "${CMAKE_BINARY_DIR}/fem/src/ElmerSolver_mpi")
ELSE()
  SET(ELMERSOLVER_BIN "${CMAKE_BINARY_DIR}/fem/src/ElmerSolver")
ENDIF()
SET(FINDNORM_BIN "${CMAKE_CURRENT_BINARY_DIR}/findnorm")
SET(MESH2D_BIN "${CMAKE_BINARY_DIR}/meshgen2d/src/Mesh2D")

MACRO(SUBDIRLIST result curdir depth)
  set(glob_pattern "*")
  FOREACH(D RANGE 1 depth)
    FILE(GLOB children RELATIVE ${curdir} ${glob_pattern}) 
    FOREACH(child ${children})
      IF(IS_DIRECTORY ${curdir}/${child})
        LIST(APPEND dirlist ${child})
      ENDIF()
    ENDFOREACH()
    set(glob_pattern "${glob_pattern}/*")
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()

SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR} 2)

SET(_test_counter 0)
FOREACH(subdir ${SUBDIRS})
  IF(NOT ${subdir} MATCHES "^contrib/"
      OR WITH_CONTRIB)
    IF(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/CMakeLists.txt")
      MATH(EXPR _test_counter "${_test_counter}+1")
      ADD_SUBDIRECTORY(${subdir})
      LIST(APPEND TEST_SUBDIRS ${subdir})
    ENDIF(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/CMakeLists.txt")
  ENDIF()
ENDFOREACH()
MESSAGE(STATUS "Found ${_test_counter} tests")

SET(INSTALL_OLD_TESTS FALSE CACHE BOOL "(deprecated) Install old test system under build directory.")

FOREACH(_subdir ${TEST_SUBDIRS})
  FILE(GLOB files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/${_subdir}/*")
  FOREACH(_file ${files})
    IF(NOT(${_file} MATCHES "CMakeLists.txt$"))
      IF(NOT(${_file} MATCHES "runtest.cmake$"))
        LIST(APPEND ELMER_TEST_FILES ${_file})
      ENDIF()
    ENDIF()
  ENDFOREACH()
ENDFOREACH()

IF(INSTALL_OLD_TESTS)
  INSTALL(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/. 
    DESTINATION ${CMAKE_BINARY_DIR}/elmertests
    FILES_MATCHING REGEX 
    ".+\\.(f90|F90|sif|mif|grd|eg|msh|good|dat|pos|result|in2d|best0)|HB|Makefile|TEST\\.RESULT|ELMERSOLVER_STARTINFO|mesh\\.*")

  # Misc. directories for meshes etc.
  INSTALL(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/diffuser_sa/diffuser/.
    DESTINATION ${CMAKE_BINARY_DIR}/elmertests/diffuser_sa/diffuser)
  INSTALL(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/diffuser_sst/diffuser/.
    DESTINATION ${CMAKE_BINARY_DIR}/elmertests/diffuser_sst/diffuser)
  INSTALL(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/diffuser_v2f/diffuser/.
    DESTINATION ${CMAKE_BINARY_DIR}/elmertests/diffuser_v2f/diffuser)

  INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/findnorm.c
    ${CMAKE_CURRENT_SOURCE_DIR}/ElmerSolver_test_how-to.txt
    DESTINATION ${CMAKE_BINARY_DIR}/elmertests)
  INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/runtests 
    DESTINATION ${CMAKE_BINARY_DIR}/elmertests)
ENDIF(INSTALL_OLD_TESTS)

# ElmerTests_package target
SET(ELMER_TESTS_PACKAGE_FILENAME "elmerfem-tests-${ELMER_FEM_REVISION}.tar.gz" CACHE STRING "Filename of test package")
MARK_AS_ADVANCED(ELMER_TESTS_PACKAGE_FILENAME)

ADD_CUSTOM_TARGET(ElmerTests
  "${CMAKE_COMMAND}" "-E" "make_directory" "${CMAKE_BINARY_DIR}/elmerfem-tests")
ADD_CUSTOM_TARGET(ElmerTests_package
  DEPENDS ElmerTests
  COMMAND "${CMAKE_COMMAND}" "-E" "tar" "cvfz" "${ELMER_TESTS_PACKAGE_FILENAME}" "elmerfem-tests"
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

# Copy the tests in a cross platform compliant manner
FOREACH(_file ${ELMER_TEST_FILES})
  IF(NOT(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${_file}))
    ADD_CUSTOM_COMMAND(TARGET ElmerTests
      COMMAND "${CMAKE_COMMAND}" "-E" "copy"
      "${CMAKE_CURRENT_SOURCE_DIR}/${_file}"
      "${CMAKE_BINARY_DIR}/elmerfem-tests/${_file}")
  ELSE()
    ADD_CUSTOM_COMMAND(TARGET ElmerTests
      COMMAND "${CMAKE_COMMAND}" "-E" "copy_directory"
      "${CMAKE_CURRENT_SOURCE_DIR}/${_file}"
      "${CMAKE_BINARY_DIR}/elmerfem-tests/${_file}")
  ENDIF()
ENDFOREACH()
