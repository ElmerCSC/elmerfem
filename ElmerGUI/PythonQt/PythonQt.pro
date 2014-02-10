#----------------------------------------------------------------------
#                   qmake project file for PythonQt
#----------------------------------------------------------------------
include(../ElmerGUI.pri)

TARGET = PythonQt
TEMPLATE = lib
CONFIG += qt release warn_off

INCLUDEPATH += $${PY_INCLUDEPATH} src
LIBPATH += $${PY_LIBPATH}
LIBS += $${PY_LIBS}

DESTDIR = lib
MOC_DIR = tmp
OBJECTS_DIR = tmp

DEFINES += PYTHONQT_EXPORTS

unix {
   CONFIG += staticlib 
}

win32 {
   DEFINES += EG_MODS
   pythonqtdll.path = $${ELMERGUI_HOME}
   pythonqtdll.files = lib\PythonQt.dll
   INSTALLS += pythonqtdll
}

HEADERS += src/PythonQt.h \
           src/PythonQtStdDecorators.h \
           src/PythonQtClassInfo.h \
           src/PythonQtImporter.h \
           src/PythonQtObjectPtr.h \
           src/PythonQtSlot.h \
           src/PythonQtStdOut.h \
           src/PythonQtMisc.h \
           src/PythonQtMethodInfo.h \
           src/PythonQtImportFileInterface.h \
           src/PythonQtConversion.h \
           src/PythonQtSignalReceiver.h \
           src/PythonQtWrapper.h \
           src/PythonQtMetaObjectWrapper.h \
           src/PythonQtCppWrapperFactory.h \
           src/PythonQtVariants.h \
           src/PythonQtVariantWrapper.h \
           src/wrapper/PythonQtWrappedVariants.h \
           src/gui/PythonQtScriptingConsole.h    \
           src/PythonQtSystem.h
  
SOURCES += src/PythonQtStdDecorators.cpp \
           src/PythonQt.cpp \
           src/PythonQtClassInfo.cpp \
           src/PythonQtImporter.cpp \
           src/PythonQtObjectPtr.cpp \
           src/PythonQtStdOut.cpp \
           src/PythonQtSlot.cpp \
           src/PythonQtMisc.cpp \
           src/PythonQtMethodInfo.cpp \
           src/PythonQtConversion.cpp \
           src/PythonQtSignalReceiver.cpp \
           src/PythonQtVariants.cpp \
           src/PythonQtVariantWrapper.cpp \
           src/PythonQtWrapper.cpp \
           src/PythonQtMetaObjectWrapper.cpp \
           src/gui/PythonQtScriptingConsole.cpp
