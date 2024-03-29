#
# this module look for qwt (http://hdf.ncsa.uiuc.edu) support
# it will define the following values
#
# QWT_INCLUDE_DIR  = where qwt.h can be found
# QWT_LIBRARY      = the library to link against qwt
# FOUND_QWT        = set to true after finding the library
#

#INCLUDE(${CMAKE_ROOT}/Modules/FindOpenGL.cmake)
#INCLUDE(${CMAKE_ROOT}/Modules/FindQt.cmake)

IF(EXISTS ${PROJECT_CMAKE}/QwtConfig.cmake)
  INCLUDE(${PROJECT_CMAKE}/QwtConfig.cmake)
ENDIF(EXISTS ${PROJECT_CMAKE}/QwtConfig.cmake)

IF(Qwt_INCLUDE_DIRS)

  FIND_PATH(QWT_INCLUDE_DIR qwt.h ${Qwt_INCLUDE_DIRS})
  FIND_LIBRARY(QWT_LIBRARY qwt HINTS ${Qwt_LIBRARY_DIRS})

ELSE(Qwt_INCLUDE_DIRS)

  FIND_PATH(QWT_INCLUDE_DIR qwt.h 
    PATHS
    /usr/include/qwt
    /usr/include/qwt6
    /usr/local/include/qwt
    /sw/include/qwt
    PATH_SUFFIXES qwt-qt6 qwt-qt5 qwt-qt4 qwt-qt
    HINTS /usr/local/opt/qwt-qt4/lib/
    HINTS /usr/local/Cellar/qwt-qt4/6.1.3_1/lib
    HINTS /usr/include/qt5/qwt6
    HINTS /msys64/mingw64/include/qwt
    HINTS /msys64/mingw64/include/qwt-qt5
    )
  FIND_LIBRARY(QWT_LIBRARY NAMES qwt qwt6-qt5 qwt-qt5 qwt-qt6 qwt-qt4 qwt-qt
    /usr/lib
    /usr/lib64
    /usr/local/lib
    /sw/lib
    HINTS /usr/local/opt/qwt-qt4/lib/
    HINTS /usr/local/Cellar/qwt-qt4/6.1.3_1/lib
    HINTS /usr/lib/qwt-qt5
    HINTS /msys64/mingw64/lib
    )

ENDIF(Qwt_INCLUDE_DIRS)

IF(QWT_INCLUDE_DIR AND QWT_LIBRARY) 
  SET(FOUND_QWT 1 CACHE BOOL "Found qwt library")
  SET(Qwt_FOUND TRUE CACHE BOOL "Found Qwt library")
ELSE(QWT_INCLUDE_DIR AND QWT_LIBRARY) 
  SET(FOUND_QWT 0 CACHE BOOL "Not found qwt library")
  SET(Qwt_FOUND FALSE CACHE BOOL "Not found qwt library")
ENDIF(QWT_INCLUDE_DIR AND QWT_LIBRARY) 

MARK_AS_ADVANCED(
  QWT_INCLUDE_DIR
  QWT_LIBRARY
  FOUND_QWT
  Qwt_FOUND
  )

