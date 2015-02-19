
########### install files ###############

FILE(GLOB TCLFILES *.tcl)
INSTALL(FILES ${TCLFILES} loadsingle DESTINATION "${ELMER_POST_DATADIR}/tcl" COMPONENT ElmerPost)


#original Makefile.am contents follow:

## nothing
