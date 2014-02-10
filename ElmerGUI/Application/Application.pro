#==============================================================================
#
#      ElmerGUI: qmake project file for Unix, Win32, and MacX
#
#==============================================================================

include(../ElmerGUI.pri)

#------------------------------------------------------------------------------
# Target:
#------------------------------------------------------------------------------
TARGET = ElmerGUI
TEMPLATE = app
CONFIG += release

#------------------------------------------------------------------------------
# Installation directory and files:
#------------------------------------------------------------------------------
target.path = $${ELMERGUI_HOME}
INSTALLS += target
edf.path = $${ELMERGUI_HOME}/edf
edf.files = edf/*
INSTALLS += edf
edf-extra.path = $${ELMERGUI_HOME}/edf-extra
edf-extra.files = edf-extra/*
INSTALLS += edf-extra

#------------------------------------------------------------------------------
# Compiler flags:
#------------------------------------------------------------------------------
CONFIG += warn_off

# QMAKE_LFLAGS += -Wl,-rpath,$$(VTKHOME)/lib/vtk-5.2

win32 {
  QMAKE_LFLAGS += /NODEFAULTLIB:library
}

#------------------------------------------------------------------------------
# Directories:
#------------------------------------------------------------------------------
DEPENDPATH += . src forms plugins vtkpost cad twod
INCLUDEPATH += .
MOC_DIR = tmp
OBJECTS_DIR = tmp
RCC_DIR = tmp
UI_DIR = tmp
DESTDIR = .

#------------------------------------------------------------------------------
# QT:
#------------------------------------------------------------------------------
QT += opengl xml script
CONFIG += uitools

#------------------------------------------------------------------------------
# MATC (see ../matc/README for more details):
#------------------------------------------------------------------------------
contains(DEFINES, EG_MATC) {
   LIBPATH += ../matc/lib
   LIBS += -lmatc
}

#------------------------------------------------------------------------------
# NETGEN (see ../netgen/README for more details):
#------------------------------------------------------------------------------
INCLUDEPATH += ../netgen/libsrc/interface
LIBPATH += ../netgen/ngcore
LIBS += -lng

#------------------------------------------------------------------------------
# QWT:
#------------------------------------------------------------------------------
contains(DEFINES, EG_QWT) {
   INCLUDEPATH += $${QWT_INCLUDEPATH}
   LIBPATH += $${QWT_LIBPATH}
   LIBS += $${QWT_LIBS}
}

#------------------------------------------------------------------------------
# VTK:
#------------------------------------------------------------------------------
contains(DEFINES, EG_VTK) {
   INCLUDEPATH += $${VTK_INCLUDEPATH}
   LIBPATH += $${VTK_LIBPATH}
   LIBS += $${VTK_LIBS}
}

#------------------------------------------------------------------------------
# OpenCASCADE:
#------------------------------------------------------------------------------
contains(DEFINES, EG_OCC) {
   contains(BITS, 64):  DEFINES += _OCC64

   unix: DEFINES += HAVE_CONFIG_H HAVE_IOSTREAM HAVE_FSTREAM HAVE_LIMITS_H
   win32: DEFINES += WNT CSFDB
   macx: DEFINED -= EG_OCC         # not supported at the moment

   INCLUDEPATH += $${OCC_INCLUDEPATH}
   LIBPATH += $${OCC_LIBPATH}
   LIBS += $${OCC_LIBS}
}

#------------------------------------------------------------------------------
# PYTHONQT (see ../PythonQt/README for more details):
#------------------------------------------------------------------------------
contains(DEFINES, EG_PYTHONQT) {
   INCLUDEPATH += $${PY_INCLUDEPATH} ../PythonQt/src
   LIBPATH += $${PY_LIBPATH} ../PythonQt/lib
   LIBS += $${PY_LIBS} -lPythonQt
}

#------------------------------------------------------------------------------
# Process info query on win32:
#------------------------------------------------------------------------------
win32: LIBS += -lpsapi

#------------------------------------------------------------------------------
# OpenGL GLU
#------------------------------------------------------------------------------
unix:  LIBS += -lGLU
#------------------------------------------------------------------------------
# Input files:
#------------------------------------------------------------------------------
HEADERS += src/bodypropertyeditor.h \
           src/boundarydivision.h \
           src/boundarypropertyeditor.h \
           src/checkmpi.h \
           src/dynamiceditor.h \
           src/edfeditor.h \
           src/egini.h \
           src/generalsetup.h \
           src/glcontrol.h \
           src/glwidget.h \
           src/helpers.h \
           src/mainwindow.h \
           src/materiallibrary.h \
           src/maxlimits.h \
           src/meshcontrol.h \
           src/meshingthread.h \
           src/meshtype.h \
           src/meshutils.h \
           src/operation.h \
           src/parallel.h \
           src/projectio.h \
           src/sifgenerator.h \
           src/sifwindow.h \
           src/solverparameters.h \
           src/summaryeditor.h \
           plugins/egconvert.h \
           plugins/egdef.h \
           plugins/egmain.h \
           plugins/egmesh.h \
           plugins/egnative.h \
           plugins/egtypes.h \
           plugins/egutils.h \
           plugins/elmergrid_api.h \
           plugins/nglib_api.h \
           plugins/tetgen.h \
           plugins/tetlib_api.h \
           twod/renderarea.h \
           twod/twodview.h \
           twod/curveeditor.h

FORMS += forms/bodypropertyeditor.ui \
         forms/boundarydivision.ui \
         forms/boundarypropertyeditor.ui \
         forms/generalsetup.ui \
         forms/glcontrol.ui \
         forms/materiallibrary.ui \
         forms/meshcontrol.ui \
         forms/parallel.ui \
         forms/solverparameters.ui \
         forms/summaryeditor.ui

SOURCES += src/bodypropertyeditor.cpp \
           src/boundarydivision.cpp \
           src/boundarypropertyeditor.cpp \
           src/checkmpi.cpp \
           src/dynamiceditor.cpp \
           src/edfeditor.cpp \
           src/egini.cpp \
           src/generalsetup.cpp \
           src/glcontrol.cpp \
           src/glwidget.cpp \
           src/helpers.cpp \
           src/main.cpp \
           src/mainwindow.cpp \
           src/materiallibrary.cpp \
           src/maxlimits.cpp \
           src/meshcontrol.cpp \
           src/meshingthread.cpp \
           src/meshtype.cpp \
           src/meshutils.cpp \
           src/operation.cpp \
           src/parallel.cpp \
           src/projectio.cpp \
           src/sifgenerator.cpp \
           src/sifwindow.cpp \
           src/solverparameters.cpp \
           src/summaryeditor.cpp \
           plugins/egconvert.cpp \
           plugins/egmain.cpp \
           plugins/egmesh.cpp \
           plugins/egnative.cpp \
           plugins/egutils.cpp \
           plugins/elmergrid_api.cpp \
           plugins/nglib_api.cpp \
           plugins/tetlib_api.cpp \
           twod/renderarea.cpp \
           twod/twodview.cpp \
           twod/curveeditor.cpp

#------------------------------------------------------------------------------
# Optional input files:
#------------------------------------------------------------------------------
contains(DEFINES, EG_QWT) {
   HEADERS += src/convergenceview.h
   SOURCES += src/convergenceview.cpp
}

contains(DEFINES, EG_VTK) {
   HEADERS += vtkpost/axes.h \
              vtkpost/featureedge.h \
              vtkpost/vtkpost.h \
              vtkpost/isosurface.h \
              vtkpost/isocontour.h \
              vtkpost/epmesh.h \
              vtkpost/colorbar.h \
              vtkpost/meshpoint.h \
              vtkpost/meshedge.h \
              vtkpost/surface.h \
              vtkpost/preferences.h \
              vtkpost/vector.h \
              vtkpost/readepfile.h \
              vtkpost/streamline.h \
              vtkpost/timestep.h \
              vtkpost/ecmaconsole.h \
              vtkpost/text.h

   FORMS += vtkpost/axes.ui \
            vtkpost/featureedge.ui \
            vtkpost/isosurface.ui \
            vtkpost/isocontour.ui \
            vtkpost/colorbar.ui \
            vtkpost/surface.ui \
            vtkpost/meshpoint.ui \
            vtkpost/meshedge.ui \
            vtkpost/preferences.ui \
            vtkpost/vector.ui \
            vtkpost/readepfile.ui \
            vtkpost/streamline.ui \
            vtkpost/timestep.ui \
            vtkpost/text.ui

   SOURCES += vtkpost/axes.cpp \
              vtkpost/featureedge.cpp \
              vtkpost/vtkpost.cpp \
              vtkpost/isosurface.cpp \
              vtkpost/isocontour.cpp \
              vtkpost/epmesh.cpp \
              vtkpost/colorbar.cpp \
              vtkpost/meshpoint.cpp \
              vtkpost/meshedge.cpp \
              vtkpost/surface.cpp \
              vtkpost/preferences.cpp \
              vtkpost/vector.cpp \
              vtkpost/readepfile.cpp \
              vtkpost/streamline.cpp \
              vtkpost/timestep.cpp \
              vtkpost/ecmaconsole.cpp \
              vtkpost/text.cpp

   contains(DEFINES, EG_MATC) {
      HEADERS += vtkpost/matc.h \
                 vtkpost/mc.h

      FORMS += vtkpost/matc.ui

      SOURCES += vtkpost/matc.cpp
   }
}

contains(DEFINES, EG_OCC) {
   HEADERS += cad/cadview.h \
              cad/cadpreferences.h

   FORMS += cad/cadpreferences.ui

   SOURCES += cad/cadview.cpp \
              cad/cadpreferences.cpp
}

#------------------------------------------------------------------------------
# Resource files:
#------------------------------------------------------------------------------
RESOURCES += ElmerGUI.qrc
win32: RC_FILE += ElmerGUI.rc
macx: RC_FILE = M3Dicon.icns

#------------------------------------------------------------------------------
# END OF FILE
#------------------------------------------------------------------------------
