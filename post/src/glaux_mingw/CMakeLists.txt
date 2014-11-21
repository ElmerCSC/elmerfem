#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../../)

add_subdirectory(src)


########### install files ###############




#original Makefile.am contents follow:

#SUBDIRS=src
