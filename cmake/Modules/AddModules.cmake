# Macros for compiling solver modules
#
# ADD_ELMER_MODULE(<name> [DIRECTORY <subdiretory>])
# Adds solver module <name> having sources in "<name>.F90".
# Optionally treat all .F90 files in <subdirectory> as sources.
#
#
# ADD_ELMER_MODULES([SKIP <name_1> <name_2> ... <skip_n>])
# Treats all directories and .F90 files from current directory as input
# for ADD_ELMER_MODULE. Skips those files and directories that equal
# to <skip_m> for some 0<m<n+1.
#
#
# MAKE_SKIP_LIST(INCLUDE name_1;name_2;...;name_n)
# Find all modules similarly as ADD_ELMER_MODULES and remove those modules from the list
# that match name_1, or name_2, or, ..., name_n. 
#
# Populates a local variable called SKIP_LIST with the resulting list.
#
# Useful when combined with ADD_ELMER_MODULES(SKIP SKIP_LIST)



MACRO(ADD_ELMER_MODULE BASENAME)
  CMAKE_PARSE_ARGUMENTS(_parsedArgs "" "" "SOURCES;DIRECTORY" "${ARGN}")

  SET(_SOURCES "")

  IF(_parsedArgs_DIRECTORY)
    FILE(GLOB DIR_SRC_FILES "${_parsedArgs_DIRECTORY}/*.F90")
    SET(_SOURCES ${DIR_SRC_FILES})
  ELSEIF(_parsedArgs_SOURCES)
    SET(_SOURCES ${_parsedArgs_SOURCES})
  ELSE()
    SET(_SOURCES ${BASENAME}.F90)
  ENDIF()

  #message(STATUS "adding library ${BASENAME} with sources ${_SOURCES}")
  ADD_LIBRARY(${BASENAME} SHARED ${_SOURCES})

  SET_TARGET_PROPERTIES(${BASENAME} PROPERTIES PREFIX "")
  SET_TARGET_PROPERTIES(${BASENAME} PROPERTIES LINKER_LANGUAGE Fortran)
  TARGET_LINK_LIBRARIES(${BASENAME} elmersolver)
  ADD_DEPENDENCIES(${BASENAME} elmersolver)
  INSTALL(TARGETS ${BASENAME} 
    LIBRARY DESTINATION "share/elmersolver/lib"
    RUNTIME DESTINATION "share/elmersolver/lib")
    #ARCHIVE DESTINATION "share/elmersolver/lib")
ENDMACRO()

MACRO(ADD_ELMER_MODULES)
  set(num_modules 0)
  CMAKE_PARSE_ARGUMENTS(_parsedArgs "" "" "SKIP" "${ARGN}")

  # Find files
  FILE(GLOB SRC_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "*.F90")

  FOREACH(FNAME ${SRC_FILES})
    GET_FILENAME_COMPONENT(BASENAME ${FNAME} NAME_WE)
    LIST(FIND _parsedArgs_SKIP ${BASENAME} FILE_INDEX)
    IF(FILE_INDEX EQUAL -1)
      IF(NOT BASENAME STREQUAL "")
        ADD_ELMER_MODULE(${BASENAME})
      ENDIF()
      MATH(EXPR num_modules "${num_modules}+1")
    ENDIF()
  ENDFOREACH()

  # Find subdirs
  FILE(GLOB SRC_DIRS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "*")

  FOREACH(DIRNAME ${SRC_DIRS})
    IF(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${DIRNAME}")
      LIST(FIND _parsedArgs_SKIP ${DIRNAME} DIR_INDEX)
      IF(DIR_INDEX EQUAL -1)
        ADD_ELMER_MODULE(${DIRNAME} DIRECTORY ${DIRNAME}) 
        MATH(EXPR num_modules "${num_modules}+1")
      ENDIF()
    ENDIF()
  ENDFOREACH()

  message(STATUS "Found ${num_modules} modules from ${CMAKE_CURRENT_SOURCE_DIR}")
ENDMACRO()

MACRO(MAKE_SKIP_LIST) # Provide list of modules to be included
  CMAKE_PARSE_ARGUMENTS(_parsedArgs "" "" "INCLUDE" "${ARGN}")

  FILE(GLOB SRC_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "*.F90")

  SET(SKIP_LIST "")

  FOREACH(FNAME ${SRC_FILES})
    GET_FILENAME_COMPONENT(BASENAME ${FNAME} NAME_WE)

    IF(NOT ${BASENAME} IN_LIST _parsedArgs_INCLUDE)
      LIST(APPEND SKIP_LIST ${BASENAME})
    ELSE()
      message(STATUS "Including ${BASENAME} module")
    ENDIF()
  ENDFOREACH()


  FILE(GLOB SRC_DIRS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "*")

  FOREACH(DIRNAME ${SRC_DIRS})
    IF(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${DIRNAME}")
      LIST(FIND ARGN ${DIRNAME} DIR_INDEX)
      IF(NOT ${DIRNAME} IN_LIST _parsedArgs_INCLUDE)
        LIST(APPEND SKIP_LIST ${DIRNAME})
      ELSE()
        message(STATUS "Including ${DIRNAME} module")
      ENDIF()
    ENDIF()
  ENDFOREACH()

ENDMACRO()
