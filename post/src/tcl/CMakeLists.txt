
########### install files ###############

FILE(GLOB TCLFILES *.tcl)
INSTALL(FILES ${TCLFILES} loadsingle DESTINATION "${ELMER_POST_DATADIR}/tcl")


#original Makefile.am contents follow:

## nothing
