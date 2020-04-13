#
# this module look for qwt (http://hdf.ncsa.uiuc.edu) support it will define the
# following values
#
# QWT_INCLUDE_DIR  = where qwt.h can be found QWT_LIBRARY      = the library to
# link against qwt FOUND_QWT        = set to true after finding the library
#

include(${CMAKE_ROOT}/Modules/FindOpenGL.cmake)
include(${CMAKE_ROOT}/Modules/FindQt.cmake)

if(EXISTS ${PROJECT_CMAKE}/QwtConfig.cmake)
  include(${PROJECT_CMAKE}/QwtConfig.cmake)
endif(EXISTS ${PROJECT_CMAKE}/QwtConfig.cmake)

if(Qwt_INCLUDE_DIRS)

  find_path(QWT_INCLUDE_DIR qwt.h ${Qwt_INCLUDE_DIRS})
  find_library(QWT_LIBRARY qwt HINTS ${Qwt_LIBRARY_DIRS})

else(Qwt_INCLUDE_DIRS)

  find_path(
    QWT_INCLUDE_DIR qwt.h
    PATHS /usr/include/qwt /usr/local/include/qwt /sw/include/qwt
    HINTS /usr/local/opt/qwt-qt4/lib/
    HINTS /usr/local/Cellar/qwt-qt4/6.1.3_1/lib)
  find_library(
    QWT_LIBRARY qwt /usr/lib /usr/local/lib /sw/lib
    HINTS /usr/local/opt/qwt-qt4/lib/
    HINTS /usr/local/Cellar/qwt-qt4/6.1.3_1/lib)

endif(Qwt_INCLUDE_DIRS)

if(QWT_INCLUDE_DIR AND QWT_LIBRARY)
  set(FOUND_QWT
      1
      CACHE BOOL "Found qwt library")
  set(Qwt_FOUND
      TRUE
      CACHE BOOL "Found Qwt library")
else(QWT_INCLUDE_DIR AND QWT_LIBRARY)
  set(FOUND_QWT
      0
      CACHE BOOL "Not found qwt library")
  set(Qwt_FOUND
      FALSE
      CACHE BOOL "Not found qwt library")
endif(QWT_INCLUDE_DIR AND QWT_LIBRARY)

mark_as_advanced(QWT_INCLUDE_DIR QWT_LIBRARY FOUND_QWT)
