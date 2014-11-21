
########### next target ###############
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(objects_STAT_SRCS
   objects.c
   objects.h
   transformations.c
)

add_library(objects STATIC ${objects_STAT_SRCS})


########### install files ###############




#original Makefile.am contents follow:

#noinst_LIBRARIES=libobjects.a
#
#libobjects_a_SOURCES=\
#	objects.c \
#	objects.h \
#	transformations.c
