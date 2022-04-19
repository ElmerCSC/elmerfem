# CMake script for finding XIOS
#

cmake_minimum_required(VERSION 2.8)

# If XIOS libraries are already defined, do nothing
IF(XIOS_LIBRARIES)
	IF(XIOS_INCLUDE_DIR)
		SET(XIOS_FOUND TRUE)
    RETURN()
  ENDIF()
ENDIF()

