
########### next target ###############

SET(glaux_STAT_SRCS
   3d.h
   font.c
   glaux.c
   image.c
   shapes.c
   teapot.c
   vect3d.c
   glaux.h
   xform.c
)

add_library(glaux STATIC ${glaux_STAT_SRCS})


########### install files ###############




#original Makefile.am contents follow:

#
#if NOTMINGW32
#
#noinst_LIBRARIES=libglaux.a
#
#libglaux_a_SOURCES=\
#	3d.h \
#	font.c \
#	glaux.c \
#	image.c \
#	shapes.c \
#	teapot.c \
#	vect3d.c \
#	glaux.h \
#	xform.c 
#
#endif
