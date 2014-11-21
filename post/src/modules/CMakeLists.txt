
########### next target ###############
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(module_STAT_SRCS
   matctcl.c
   misc.c
   readfile.c
   transforms.c
)

add_library(module STATIC ${module_STAT_SRCS})


########### install files ###############




#original Makefile.am contents follow:

#noinst_LIBRARIES=libmodule.a
#
#libmodule_a_SOURCES=\
#	matctcl.c \
#	misc.c \
#	readfile.c \
#	transforms.c 
