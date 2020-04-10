
########### next target ###############
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

# find_package(X11)

# if (X11_FOUND)
# 	if (APPLE)
# 		find_path (X11_GL_INCLUDE_PATH GL/glx.h ${X11_INC_SEARCH_PATH})
# 		if (NOT X11_GL_INCLUDE_PATH)
# 			message (WARNING "Could not find GL/glx.h")
# 			set (X11_FOUND FALSE)
# 		endif (NOT X11_GL_INCLUDE_PATH)
# 		set (X11_INCLUDE_DIR ${X11_INCLUDE_DIR} ${X11_GL_INCLUDE_PATH})

# 		find_library (X11_GL_LIB GL ${X11_LIB_SEARCH_PATH})
# 		if (NOT X11_GL_LIB)
# 			message (WARNING "Could not find libGL.dylib")
# 			set (X11_FOUND FALSE)
# 		endif (NOT X11_GL_LIB)

# 		find_library (X11_GLU_LIB GLU ${X11_LIB_SEARCH_PATH})
# 		if (NOT X11_GLU_LIB)
# 			message (WARNING "Could not find libGLU.dylib")
# 			set (X11_FOUND FALSE)
# 		endif (NOT X11_GLU_LIB)
# 	else ()
# 		set (X11_INCLUDE_DIR ${X11_INCLUDE_DIR} ${OPENGL_INCLUDE_DIR})
# 		set (X11_GL_LIB ${OPENGL_gl_LIBRARY})
# 		set (X11_GLU_LIB ${OPENGL_glu_LIBRARY})
# 	endif ()
# endif (X11_FOUND)

# if(APPLE)
#    find_path (X11_GL_INCLUDE_PATH GL/glx.h
#       HINTS /opt/X11/include/
#    )
# endif(APPLE)




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
