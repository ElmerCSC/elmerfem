PROJECT(ElmerIceTests Fortran C)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET(ELMERGRID_BIN "${CMAKE_BINARY_DIR}/elmergrid/src/ElmerGrid")
IF(WITH_MPI)
  SET(ELMERSOLVER_BIN "${CMAKE_BINARY_DIR}/fem/src/ElmerSolver_mpi")
ELSE()
  SET(ELMERSOLVER_BIN "${CMAKE_BINARY_DIR}/fem/src/ElmerSolver")
ENDIF()
#SET(FINDNORM_BIN "${CMAKE_CURRENT_BINARY_DIR}/findnorm")

MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
        LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()

SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR})

FOREACH(subdir ${SUBDIRS})
  IF(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/CMakeLists.txt")
    ADD_SUBDIRECTORY(${subdir})
    LIST(APPEND TEST_SUBDIRS ${subdir})
  ENDIF(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/CMakeLists.txt")
ENDFOREACH()

#SET(INSTALL_OLD_TESTS FALSE CACHE BOOL "(deprecated) Install old test system under build directory.")

FOREACH(_subdir ${TEST_SUBDIRS})
  FILE(GLOB files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/${_subdir}/*")
  FOREACH(_file ${files})
    IF(NOT(${_file} MATCHES "CMakeLists.txt$"))
      IF(NOT(${_file} MATCHES "runTest.cmake$"))
        LIST(APPEND ELMERICE_TEST_FILES ${_file})
      ENDIF()
    ENDIF()
  ENDFOREACH()
ENDFOREACH()

# ElmerTests_package target
SET(ELMERICE_TESTS_PACKAGE_FILENAME "elmerice-tests-${ELMER_FEM_REVISION}.tar.gz" CACHE STRING "Filename of test package")
MARK_AS_ADVANCED(ELMERICE_TESTS_PACKAGE_FILENAME)

ADD_CUSTOM_TARGET(ElmerIceTests
  "${CMAKE_COMMAND}" "-E" "make_directory" "${CMAKE_BINARY_DIR}/elmerice-tests")
ADD_CUSTOM_TARGET(ElmerIceTests_package
  DEPENDS ElmerIceTests
  COMMAND "${CMAKE_COMMAND}" "-E" "tar" "cvfz" "${ELMERICE_TESTS_PACKAGE_FILENAME}" "elmerice-tests"
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

# Copy the tests in a cross platform compliant manner
FOREACH(_file ${ELMERICE_TEST_FILES})
  IF(NOT(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${_file}))
    ADD_CUSTOM_COMMAND(TARGET ElmerIceTests
      COMMAND "${CMAKE_COMMAND}" "-E" "copy"
      "${CMAKE_CURRENT_SOURCE_DIR}/${_file}"
      "${CMAKE_BINARY_DIR}/elmerice-tests/${_file}")
  ELSE()
    ADD_CUSTOM_COMMAND(TARGET ElmerIceTests
      COMMAND "${CMAKE_COMMAND}" "-E" "copy_directory"
      "${CMAKE_CURRENT_SOURCE_DIR}/${_file}"
      "${CMAKE_BINARY_DIR}/elmerice-tests/${_file}")
  ENDIF()
ENDFOREACH()

	

#IF(NOT(CMAKE_BUILD_TYPE MATCHES "Rel"))
#INSTALL(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/. #
#		  DESTINATION ${CMAKE_BINARY_DIR}/elmericetests#
#		  FILES_MATCHING REGEX 
#		  ".+\\.(f90|F90|sif|mif|grd|eg|msh|good|dat|pos|result|in2d|best0|txt)|HB|Makefile|TEST\\.RESULT|ELMERSOLVER_STARTINFO|mesh\\.*|DATA|PROG")
#
#INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/Compare.f90
#	      ${CMAKE_CURRENT_SOURCE_DIR}/ElmerIceSolver_test_how-to.txt
#	     	  DESTINATION ${CMAKE_BINARY_DIR}/elmericetests)
#INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/runTest.sh
#		  DESTINATION ${CMAKE_BINARY_DIR}/elmericetests)
#ENDIF()
