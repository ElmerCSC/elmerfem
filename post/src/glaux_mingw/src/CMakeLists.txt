
########### install files ###############

#install(FILES  GL/glaux.h DESTINATION )

SET(libglaux_a_SOURCES font.c  image.c teapot.c tkdib.c tkimage.c xform.c glaux.c shapes.c tk.c tkfont.c vect3d.c)

add_library(glaux STATIC ${libglaux_a_SOURCES})

#original Makefile.am contents follow:

#nobase_include_HEADERS=GL/glaux.h
#
#lib_LIBRARIES=libglaux.a
#
#libglaux_a_SOURCES=font.c \
#                   image.c \
#                   teapot.c \
#                   tkdib.c \
#                   tkimage.c \
#                   xform.c \
#                   glaux.c \
#                   shapes.c \
#                   tk.c \
#                   tkfont.c \
#                   vect3d.c
