
#ADD_EXECUTABLE(findnorm findnorm_cmake.c)

SET(ELMERGRID_BIN "${CMAKE_BINARY_DIR}/elmergrid/src/ElmerGrid")
IF(WITH_MPI)
  SET(ELMERSOLVER_BIN "${CMAKE_BINARY_DIR}/fem/src/ElmerSolver_mpi")
ELSE()
  SET(ELMERSOLVER_BIN "${CMAKE_BINARY_DIR}/fem/src/ElmerSolver")
ENDIF()
SET(FINDNORM_BIN "${CMAKE_CURRENT_BINARY_DIR}/findnorm")
SET(MESH2D_BIN "${CMAKE_BINARY_DIR}/meshgen2d/src/Mesh2D")

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/runtests.in.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/runtests)

# Environment variables for ElmerSolver
# SET(ENV{ELMER_LIB} "${CMAKE_BINARY_DIR}/fem/src")
# SET(ENV{ELMER_HOME} "${CMAKE_BINARY_DIR}/fem/src")

#INCLUDE(test_subdirs.cmake)

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
  ENDIF(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/CMakeLists.txt")
ENDFOREACH()

IF(NOT(CMAKE_BUILD_TYPE MATCHES "Rel"))
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
ENDIF(NOT(CMAKE_BUILD_TYPE MATCHES "Rel"))
