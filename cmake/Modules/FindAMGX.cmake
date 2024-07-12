# CMake script for finding AMGX
# Thomas Zwinger, CSC - IT Center for Science Ltd.
# 2020/03
#

#  AMGX_INCLUDE_DIR  - user modifiable choice of where to AMGX include dir
#  AMGX_LIBRARY      - user modifiable choice of where AMGX library is

# his module returns these variables for the rest of the project to use.
#
#  AMGX_FOUND              - True if AMGX found 
#  AMGX_LIBRARY           -  AMGX library
#  AMGX_INCLUDE_DIR       - AMGX include dir.
#  CUDA_LIBRARIES         - needed cuda libraries
#  CUDA_LIBDIR             - needed cuda library directory

#INCLUDE(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)

MESSAGE ("----------------------")
MESSAGE ("-- AMGX + CUDA:")
FIND_PACKAGE(CUDA)
MESSAGE("-- Cuda libraries: " ${CUDA_LIBRARIES})

SET(AMGX_FOUND FALSE)

FIND_PATH(AMGX_INCLUDE_DIR amgx_c.h
  HINTS 
  "${AMGXINCLUDE}" "${AMGX_ROOT}/include"
  )




FIND_LIBRARY(AMGX_LIBRARY
  NAMES amgx libamgx.a
  NAMES_PER_DIR
  HINTS "${AMGX_ROOT}/lib" "${AMGXLIB}"
  REQUIRED
  )

IF (AMGX_LIBRARY AND AMGX_INCLUDE_DIR)
  SET (AMGX_FOUND TRUE)
ENDIF()


IF (AMGX_FOUND)
  SET(AMGX_INCLUDE_DIRS "${AMGX_INCLUDE_DIR}")
  SET(CUDA_LIBRARIES ${CUDA_LIBRARIES} ${CUDA_cusparse_LIBRARY} ${CUDA_cublas_LIBRARY} ${CUDA_cusolver_LIBRARY} ${CUDA_nvToolsExt_LIBRARY})
  SET(AMGX_LIBRARIES ${AMGX_LIBRARY} ${CUDA_LIBRARIES})
  MESSAGE ("-- AMGX found")
  MESSAGE ("-- AMGX_ROOT= ${AMGX_ROOT}")
  MESSAGE ("----------------------")
ELSE()
    MESSAGE (FATAL_ERROR, "AMGX not found")
ENDIF()

MARK_AS_ADVANCED(
  AMGX_FOUND
  AMGX_INCLUDE_DIR
  AMGX_LIBRARY
  CUDA_LIBRARIES
  CUDA_LIBDIR)
