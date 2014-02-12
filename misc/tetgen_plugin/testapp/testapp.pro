TEMPLATE = app
TARGET = testmain
DEPENDPATH += .
INCLUDEPATH += . ../plugin
QMAKE_CXXFLAGS += -DTETLIBRARY
CONFIG += release console
SOURCES += main.cpp
win32: {
  INSTALLS += target
  target.path = .
}
