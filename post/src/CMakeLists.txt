include_directories(${CMAKE_CURRENT_BINARY_DIR}/graphics)
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR}/camera)


INCLUDE_DIRECTORIES(${TCL_INCLUDE_PATH} ${TK_INCLUDE_PATH})


IF(NOT(WIN32))
add_subdirectory(tk)
ENDIF(NOT(WIN32))

IF(NOT(MINGW))
  add_subdirectory(glaux)
ELSE()
  add_subdirectory(glaux_mingw)
ENDIF(NOT(MINGW))
add_subdirectory(camera)
add_subdirectory(graphics)
add_subdirectory(modules)
add_subdirectory(sico2elmer)
add_subdirectory(objects)
add_subdirectory(visuals)
add_subdirectory(tcl)
add_subdirectory(elements)


########### next target ###############

IF(NOT(WIN32))
  SET(ElmerPost_SRCS
    ElmerPost.c
    elmerpost.h
    geometry.c
    geometry.h
    LUDecomp.c
    second.c
    cplusplus.cpp
    fttext.cpp
    teksti.c)
ELSEIF(MINGW AND WIN32)
  SET(ElmerPost_SRCS
    ElmerPost.c
    elmerpost.h
    geometry.c
    geometry.h
    LUDecomp.c
    second.c
    cplusplus.cpp
    fttext.cpp
    ElmerPost-res.rc
    )
ENDIF()


IF(MINGW AND WIN32)
  SET(CMAKE_RC_COMPILER_INIT windres)
  ENABLE_LANGUAGE(RC)
  SET(CMAKE_RC_COMPILE_OBJECT
    "<CMAKE_RC_COMPILER> <FLAGS> -O coff <DEFINES> -i <SOURCE> -o <OBJECT>")
ENDIF()

add_executable(ElmerPost ${ElmerPost_SRCS})

target_link_libraries(ElmerPost camera elements module objects visuals graphics
  glaux ${OPENGL_LIBRARIES} ${TCL_LIBRARY} ${TK_LIBRARY} matc m)

IF(NOT(WIN32))
  SET_TARGET_PROPERTIES(ElmerPost PROPERTIES INSTALL_RPATH ${ELMERSOLVER_RPATH_STRING})
ENDIF()

IF(WIN32)
  target_link_libraries(ElmerPost gdi32 user32 advapi32)
ENDIF(WIN32)

install(TARGETS ElmerPost DESTINATION bin COMPONENT ElmerPost)


########### next target ###############
IF(NOT(WIN32))
  SET(QueryGLXExt_SRCS
    QueryGLXExt.c
    )

  add_executable(QueryGLXExt ${QueryGLXExt_SRCS})

  target_link_libraries(QueryGLXExt ${OPENGL_LIBRARIES})

  install(TARGETS QueryGLXExt DESTINATION bin)
ENDIF(NOT(WIN32))

install(DIRECTORY "lib" DESTINATION ${ELMER_POST_DATADIR} COMPONENT ElmerPost)


########### next target ###############





#original Makefile.am contents follow:

#SUBDIRS=tk glaux camera graphics modules sico2elmer objects visuals tcl elements 
#
#if NOTWINDOWS
#
#bin_PROGRAMS=ElmerPost QueryGLXExt
#
#ElmerPost_SOURCES=\
#	ElmerPost.c \
#	elmerpost.h \
#	geometry.c \
#	geometry.h \
#	LUDecomp.c \
#	second.c \
#	cplusplus.cpp \
#	fttext.cpp \
#	teksti.c 
#
#ElmerPost_LDADD=\
#	./camera/libcamera.a \
#	./elements/libelements.a \
#	./modules/libmodule.a \
#	./objects/libobjects.a \
#	./visuals/libvisuals.a \
#	./graphics/libgraphics.a \
#	./glaux/libglaux.a \
#	./tk/libtk.a
#
#QueryGLXExt_SOURCES=QueryGLXExt.c 
#
#all-local: screensave.so
#clean-local:
#	rm -f *.so
#
#EXTRA_DIST=screensave.c
#
#screensave.$(OBJEXT): screensave.c
#	$(COMPILE) -c $< -fPIC
#
#screensave.so: screensave.$(OBJEXT)
#	$(CC) $(CFLAGS) $(LDFLAGS) $(SH_LDFLAGS) -o $@ $< $(LIBS)
#
#else
#
## WINDOWS
#
#.rc.o:
#	windres $^ -o $@
#%.o : %.rc
#	windres $^ -o $@
#
#bin_PROGRAMS=ElmerPost 
#
#ElmerPost_SOURCES=\
#	ElmerPost.c \
#	elmerpost.h \
#	geometry.c \
#	geometry.h \
#	LUDecomp.c \
#	second.c \
#	cplusplus.cpp \
#	fttext.cpp \
#	ElmerPost-res.rc
#
#if NOTMINGW32
#ElmerPost_LDADD=\
#	./camera/libcamera.a \
#	./elements/libelements.a \
#	./modules/libmodule.a \
#	./objects/libobjects.a \
#	./visuals/libvisuals.a \
#	./graphics/libgraphics.a \
#	./glaux/libglaux.a \
#	./tk/libtk.a
#else
#ElmerPost_LDADD=\
#	./camera/libcamera.a \
#	./elements/libelements.a \
#	./modules/libmodule.a \
#	./objects/libobjects.a \
#	./visuals/libvisuals.a \
#	./graphics/libgraphics.a \
#	-lglaux  \
#	-lgdi32 \
#	-luser32 \
#	-ladvapi32
#endif
#endif
#
#ElmerPost_LDFLAGS=$(SH_EXPALL_FLAG)
#
#ELMER_POST_DATADIR=$(prefix)/share/elmerpost
#
#install-data-local:
#	@echo "Installing modules to ${ELMER_POST_DATADIR}/modules"
#
#	mkdir -p $(DESTDIR)$(ELMER_POST_DATADIR)/modules
#	-cp screensave.so $(DESTDIR)$(ELMER_POST_DATADIR)/modules 
#
#	@echo "Installing help to ${ELMER_POST_DATADIR}/help"
#	mkdir -p $(DESTDIR)$(ELMER_POST_DATADIR)/help
#	cp -r help/* $(DESTDIR)$(ELMER_POST_DATADIR)/help
#	rm -rf `find $(DESTDIR)$(ELMER_POST_DATADIR)/help -name .svn`
#
#	@echo "Installing tcl scripts to ${ELMER_POST_DATADIR}/tcl"
#	mkdir -p $(DESTDIR)$(ELMER_POST_DATADIR)/tcl
#	cp -r tcl/* $(DESTDIR)$(ELMER_POST_DATADIR)/tcl
#	rm -rf `find $(DESTDIR)$(ELMER_POST_DATADIR)/tcl -name .svn`
#
#	@echo "Installing initialization scripts to ${ELMER_POST_DATADIR}/lib"
#	mkdir -p $(DESTDIR)$(ELMER_POST_DATADIR)/lib
#	cp -r lib/* $(DESTDIR)$(ELMER_POST_DATADIR)/lib
#	rm -rf `find $(DESTDIR)$(ELMER_POST_DATADIR)/lib -name .svn`
#
#	@echo "Installing fonts to $(ELMER_POST_DATADIR)/fonts"
#	mkdir -p $(DESTDIR)$(ELMER_POST_DATADIR)/fonts
#	cp -r fonts/* $(DESTDIR)$(ELMER_POST_DATADIR)/fonts
#	rm -rf `find $(DESTDIR)$(ELMER_POST_DATADIR)/fonts -name .svn`
##
## Copy the help and tcl
##
#dist-hook:
#	cp -r tcl $(distdir)
#	cp -r help $(distdir)
#	cp -r lib $(distdir)
#	cp -r fonts $(distdir)
#	rm -rf `find $(distdir)/ -name .svn`
#	rm -rf `find $(distdir)/ -name CVS`
#	rm -rf `find $(distdir)/ -name \#*`
#	rm -rf `find $(distdir)/ -name *~`
#
#
#
#
#
#
