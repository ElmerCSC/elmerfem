ELMER_HOME=$$(ELMER_HOME)

isEmpty(ELMER_HOME) {
   error("ELMER_HOME is undefined")
}

TEMPLATE = lib
TARGET = tetplugin
DEPENDPATH += .
INCLUDEPATH += .
CONFIG += release dll warn_off
QMAKE_CXXFLAGS += -DTETLIBRARY

unix: {
   SOURCES_NOOPTIMIZE = predicates.cxx
   nooptimize.name = nooptimize
   nooptimize.input = SOURCES_NOOPTIMIZE
   nooptimize.dependency_type = TYPE_C
   nooptimize.variable_out = OBJECTS
   nooptimize.output = ${QMAKE_VAR_OBJECTS_DIR}${QMAKE_FILE_IN_BASE}$${first(QMAKE_EXT_OBJ)}
   nooptimize.commands = $${QMAKE_CXX} -fPIC -O0 $(INCPATH) -c ${QMAKE_FILE_IN} -o ${QMAKE_FILE_OUT}
   QMAKE_EXTRA_COMPILERS += nooptimize   
}

HEADERS += tetgen.h
SOURCES += predicates.cxx tetgen.cxx ElmerAPI.cpp
unix: SOURCES -= predicates.cxx

target.path = $${ELMER_HOME}/lib
INSTALLS += target
