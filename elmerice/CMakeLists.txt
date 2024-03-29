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
  FIND_PACKAGE(PROJ)
  FIND_PACKAGE(FORTRANGIS)

  IF(PROJ_FOUND AND FORTRANGIS_FOUND)
    SET(HAVE_PROJ TRUE)
    SET(HAVE_FORTRANGIS TRUE)
    MARK_AS_ADVANCED(HAVE_PROJ HAVE_FORTRANGIS)
    ADD_DEFINITIONS(-DHAVE_PROJ)
    MESSAGE(STATUS "  PROJ & FORTRANGIS: " "${PROJ_FOUND} & ${FORTRANGIS_FOUND}")
    MESSAGE(STATUS "  PROJ_INC:   " "${PROJ_INCLUDE_DIR}")
    MESSAGE(STATUS "  PROJ_LIB:   " "${PROJ_LIBRARY}")
    MESSAGE(STATUS "  FORTRANGIS_INC:  " "${FORTRANGIS_INCLUDE_DIR}")
    MESSAGE(STATUS "  FORTRANGIS_LIB:  " "${FORTRANGIS_LIBRARY}")

    MESSAGE(STATUS "Compile ProjUtils with Proj4 support")
  ELSE()
    MESSAGE(STATUS "Library not found: >PROJ_FOUND< = ${PROJ_FOUND} and/or >FORTRANGIS_FOUND< =  ${FORTRANGIS_FOUND}")
    MESSAGE(WARNING "  \n Missing: >PROJ_INCLUDE_DIR< , >PROJ_LIBRARY< , >FORTRANGIS_INCLUDE_DIR< , >FORTRANGIS_LIBRARY< \n to compile ProjUtils support for Proj4")
  ENDIF()
ENDIF()

ADD_SUBDIRECTORY(Utils)
ADD_SUBDIRECTORY(UserFunctions)
ADD_SUBDIRECTORY(Solvers)
ADD_SUBDIRECTORY(Tests)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Solvers/Permafrost/permafrostmaterialdb.dat
  ${CMAKE_CURRENT_BINARY_DIR}/Solvers/Permafrost/permafrostmaterialdb.dat COPYONLY)

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/Solvers/Permafrost/permafrostmaterialdb.dat  DESTINATION "share/elmersolver/lib")
