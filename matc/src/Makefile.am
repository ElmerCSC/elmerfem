#  
# - which are needed? 
# - should it be <matc.h> or <matc/matc.h>
# 
nobase_include_HEADERS=elmer/matc.h elmer/fnames.h elmer/gra.h 
noinst_HEADERS=str.h

COMMON_SRCS =  \
	c3d.c \
	clip.c \
	dri_ps.c \
	eig.c \
	error.c \
	eval.c \
	files.c \
	funcs.c \
	gra.c \
	gra_com.c \
	jacobi.c \
	lists.c \
	lu.c \
	matc.c \
	matrix.c \
	oper.c \
	optim.c \
	parser.c \
	str.c \
	urand.c \
	variable.c 

bin_PROGRAMS = matc

matc_SOURCES = main.c

matc_LDADD = libmatc.a

lib_LIBRARIES = libmatc.a
libmatc_a_SOURCES = $(COMMON_SRCS)
