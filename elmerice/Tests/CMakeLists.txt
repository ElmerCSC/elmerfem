project(ElmerIceTests Fortran C)
cmake_minimum_required(VERSION 2.8)

set(ELMERGRID_BIN "${CMAKE_BINARY_DIR}/elmergrid/src/ElmerGrid")
if(WITH_MPI)
  set(ELMERSOLVER_BIN "${CMAKE_BINARY_DIR}/fem/src/ElmerSolver_mpi")
else()
  set(ELMERSOLVER_BIN "${CMAKE_BINARY_DIR}/fem/src/ElmerSolver")
endif()
# SET(FINDNORM_BIN "${CMAKE_CURRENT_BINARY_DIR}/findnorm")

macro(SUBDIRLIST result curdir)
  file(
    GLOB children
    RELATIVE ${curdir}
    ${curdir}/*)
  set(dirlist "")
  foreach(child ${children})
    if(IS_DIRECTORY ${curdir}/${child})
      list(APPEND dirlist ${child})
    endif()
  endforeach()
  set(${result} ${dirlist})
endmacro()

subdirlist(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR})

foreach(subdir ${SUBDIRS})
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/CMakeLists.txt")
    add_subdirectory(${subdir})
    list(APPEND TEST_SUBDIRS ${subdir})
  endif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/CMakeLists.txt")
endforeach()

# SET(INSTALL_OLD_TESTS FALSE CACHE BOOL "(deprecated) Install old test system
# under build directory.")

foreach(_subdir ${TEST_SUBDIRS})
  file(
    GLOB files
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    "${CMAKE_CURRENT_SOURCE_DIR}/${_subdir}/*")
  foreach(_file ${files})
    if(NOT (${_file} MATCHES "CMakeLists.txt$"))
      if(NOT (${_file} MATCHES "runTest.cmake$"))
        list(APPEND ELMERICE_TEST_FILES ${_file})
      endif()
    endif()
  endforeach()
endforeach()

# ElmerTests_package target
set(ELMERICE_TESTS_PACKAGE_FILENAME
    "elmerice-tests-${ELMER_FEM_REVISION}.tar.gz"
    CACHE STRING "Filename of test package")
mark_as_advanced(ELMERICE_TESTS_PACKAGE_FILENAME)

add_custom_target(ElmerIceTests "${CMAKE_COMMAND}" "-E" "make_directory"
                                "${CMAKE_BINARY_DIR}/elmerice-tests")
add_custom_target(
  ElmerIceTests_package
  DEPENDS ElmerIceTests
  COMMAND "${CMAKE_COMMAND}" "-E" "tar" "cvfz"
          "${ELMERICE_TESTS_PACKAGE_FILENAME}" "elmerice-tests"
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

# Copy the tests in a cross platform compliant manner
foreach(_file ${ELMERICE_TEST_FILES})
  if(NOT (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${_file}))
    add_custom_command(
      TARGET ElmerIceTests
      COMMAND
        "${CMAKE_COMMAND}" "-E" "copy" "${CMAKE_CURRENT_SOURCE_DIR}/${_file}"
        "${CMAKE_BINARY_DIR}/elmerice-tests/${_file}")
  else()
    add_custom_command(
      TARGET ElmerIceTests
      COMMAND
        "${CMAKE_COMMAND}" "-E" "copy_directory"
        "${CMAKE_CURRENT_SOURCE_DIR}/${_file}"
        "${CMAKE_BINARY_DIR}/elmerice-tests/${_file}")
  endif()
endforeach()

# IF(NOT(CMAKE_BUILD_TYPE MATCHES "Rel")) INSTALL(DIRECTORY
# ${CMAKE_CURRENT_SOURCE_DIR}/. # DESTINATION ${CMAKE_BINARY_DIR}/elmericetests#
# FILES_MATCHING REGEX ".+\\.(f90|F90|sif|mif|grd|eg|msh|good|dat|pos|result|in2
# d|best0|txt)|HB|Makefile|TEST\\.RESULT|ELMERSOLVER_STARTINFO|mesh\\.*|DATA|PRO
# G")
#
# INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/Compare.f90
# ${CMAKE_CURRENT_SOURCE_DIR}/ElmerIceSolver_test_how-to.txt DESTINATION
# ${CMAKE_BINARY_DIR}/elmericetests) INSTALL(PROGRAMS
# ${CMAKE_CURRENT_BINARY_DIR}/runTest.sh DESTINATION
# ${CMAKE_BINARY_DIR}/elmericetests) ENDIF()
