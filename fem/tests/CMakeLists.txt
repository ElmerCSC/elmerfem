
ADD_EXECUTABLE(findnorm findnorm_cmake.c)

SET(ELMERGRID_BIN "${CMAKE_BINARY_DIR}/elmergrid/src/ElmerGrid")
SET(ELMERSOLVER_BIN "${CMAKE_BINARY_DIR}/fem/src/ElmerSolver")
SET(FINDNORM_BIN "${CMAKE_CURRENT_BINARY_DIR}/findnorm")

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/runtests.in.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/runtests)

SET(ENV{ELMER_LIB} "{CMAKE_BINARY_DIR}/fem/src")
SET(ENV{ELMER_HOME} "{CMAKE_BINARY_DIR}/fem/src")

ADD_SUBDIRECTORY(1dtests)
ADD_SUBDIRECTORY(1sttime)
ADD_SUBDIRECTORY(2ndtime)
ADD_SUBDIRECTORY(adaptivity2)

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
