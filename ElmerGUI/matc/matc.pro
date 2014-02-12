#----------------------------------------------------------------------
#                 qmake project file for MATC library
#----------------------------------------------------------------------
include(../ElmerGUI.pri)

TARGET = matc
TEMPLATE = lib
CONFIG -= qt debug
CONFIG += staticlib release warn_off
DESTDIR = lib
OBJECTS_DIR = obj
win32 {
   DEFINES = WIN32 _CRT_SECURE_NO_WARNINGS 
   SOURCES = src\*.c
} else {
   SOURCES = src/*.c
}
