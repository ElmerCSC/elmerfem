#
# this module look for qwt (http://hdf.ncsa.uiuc.edu) support
# it will define the following values
#
# QWT_INCLUDE_DIR  = where qwt.h can be found
# QWT_LIBRARY      = the library to link against qwt
# FOUND_QWT        = set to true after finding the library
#

INCLUDE(${CMAKE_ROOT}/Modules/FindOpenGL.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/FindQt.cmake)

IF(EXISTS ${PROJECT_CMAKE}/QwtConfig.cmake)
  INCLUDE(${PROJECT_CMAKE}/QwtConfig.cmake)
ENDIF(EXISTS ${PROJECT_CMAKE}/QwtConfig.cmake)

IF(Qwt_INCLUDE_DIRS)

  FIND_PATH(QWT_INCLUDE_DIR qwt.h ${Qwt_INCLUDE_DIRS})
  FIND_LIBRARY(QWT_LIBRARY qwt ${Qwt_LIBRARY_DIRS})

ELSE(Qwt_INCLUDE_DIRS)

  FIND_PATH(QWT_INCLUDE_DIR qwt.h 
    /usr/include/qwt
    /usr/local/include/qwt
    /sw/include/qwt
    )
  FIND_LIBRARY(QWT_LIBRARY qwt 
    /usr/lib
    /usr/local/lib
    /sw/lib
    )

ENDIF(Qwt_INCLUDE_DIRS)

IF(QWT_INCLUDE_DIR AND QWT_LIBRARY) 
  SET(FOUND_QWT 1 CACHE BOOL "Found qwt library")
ELSE(QWT_INCLUDE_DIR AND QWT_LIBRARY) 
  SET(FOUND_QWT 0 CACHE BOOL "Not found qwt library")
ENDIF(QWT_INCLUDE_DIR AND QWT_LIBRARY) 

MARK_AS_ADVANCED(
  QWT_INCLUDE_DIR
  QWT_LIBRARY
  FOUND_QWT
  )

