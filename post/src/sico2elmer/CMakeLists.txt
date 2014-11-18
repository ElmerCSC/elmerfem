
########### next target ###############

SET(sico2elmer_SRCS
   sico2elmerc.c
   sico2elmer.f90
   sico2elmer.h
)

add_executable(sico2elmer ${sico2elmer_SRCS})

target_link_libraries(sico2elmer)
SET_TARGET_PROPERTIES(sico2elmer PROPERTIES LINKER_LANGUAGE Fortran)
install(TARGETS sico2elmer DESTINATION bin)


########### install files ###############




#original Makefile.am contents follow:

#if NOTWINDOWS
#
#bin_PROGRAMS=sico2elmer
#
#sico2elmer_SOURCES=\
#	sico2elmerc.c \
#	sico2elmer.f90 \
#	sico2elmer.h 
#
#endif
