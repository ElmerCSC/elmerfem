TEMPLATE = app
TARGET = ElmerClips
DEPENDPATH += . src
INCLUDEPATH += . src

win32 {
   INCLUDEPATH += D:/ffmpeg/include src/win32
   QMAKE_LIBDIR += D:/ffmpeg/bin
   LIBS += -lavcodec -lavutil -lswscale
   DESTDIR = ElmerClips
}

unix {
  LIBS += -lavcodec -lavutil -lswscale
}

CONFIG += release

HEADERS += src/preview.h src/encoder.h
SOURCES += src/main.cpp src/preview.cpp src/encoder.cpp
RESOURCES += ElmerClips.qrc
RC_FILE += ElmerClips.rc
