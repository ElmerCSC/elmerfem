INCLUDE(CheckIncludeFiles)

SET(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/elmerice/cmake/Modules")
#TODO: is this proper way of choosing shared library API?
IF(WIN32)
  UNSET(HAVE_DLOPEN_API)
  SET(HAVE_LOADLIBRARY_API TRUE)
ELSEIF(UNIX)
  SET(HAVE_DLOPEN_API TRUE)
  SET(HAVE_DLSYM TRUE)
ENDIF()


INCLUDE(CheckTypeSize)
CHECK_TYPE_SIZE(long LONG_SIZE)
IF(LONG_SIZE EQUAL 4)
  SET(OFF_KIND "selected_int_kind(9)")
ELSEIF(LONG_SIZE EQUAL 8)
  SET(OFF_KIND "selected_int_kind(18)")
ENDIF()

SET(ELMER_SOLVER_HOME "${CMAKE_INSTALL_PREFIX}/share/elmersolver" CACHE PATH "")

SET(SHL_EXTENSION ".so")
IF(WIN32)
  SET(SHL_EXTENSION ".dll")
ENDIF()
IF(APPLE)
  SET(SHL_EXTENSION ".dylib")
ENDIF()

INCLUDE_DIRECTORIES("${PROJECT_BINARY_DIR}")

IF (WITH_PROJ4)
	UNSET(WITH_PROJ4 CACHE)	
 MESSAGE(FATAL_ERROR " DEPRECATED USE WITH_PROJ=TRUE.\n"
                     " See elmerice/Utils/Documentation/ProjUtils.md" )
ENDIF()

IF (WITH_PROJ)
  FIND_PACKAGE(PROJ)

  MESSAGE(STATUS "------------------------------------------------")
  MESSAGE(STATUS "Elmer/Ice package: Looking for Proj library     ")
  IF(PROJ_FOUND)
    SET(HAVE_PROJ TRUE)
    MARK_AS_ADVANCED(HAVE_PROJ)
    ADD_DEFINITIONS(-DHAVE_PROJ)
    MESSAGE(STATUS "  PROJ LIBRARY FOUND=" " ${PROJ_FOUND} ")
    MESSAGE(STATUS "  PROJ_LIBRARIES:   " "${PROJ_LIBRARIES}")
    MESSAGE(STATUS "Compile ProjUtils with Proj support")
  ELSE()
    MESSAGE(STATUS "Library not found: >PROJ_FOUND< = ${PROJ_FOUND}")
    MESSAGE(WARNING "Missing: >PROJ_LIBRARIES< to compile ProjUtils support for Proj")
  ENDIF()
  MESSAGE(STATUS "------------------------------------------------")
ENDIF()

ADD_SUBDIRECTORY(Utils)
ADD_SUBDIRECTORY(UserFunctions)
ADD_SUBDIRECTORY(Solvers)
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(Tests)
ENDIF()

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Solvers/Permafrost/permafrostmaterialdb.dat
  ${CMAKE_CURRENT_BINARY_DIR}/Solvers/Permafrost/permafrostmaterialdb.dat COPYONLY)

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/Solvers/Permafrost/permafrostmaterialdb.dat  DESTINATION "share/elmersolver/lib")
