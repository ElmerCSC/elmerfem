#==============================================================================
#
#                       ElmerGUI: configuration file
#
#==============================================================================

#------------------------------------------------------------------------------
# Optional components (undefine or comment out to exclude from compilation):
#------------------------------------------------------------------------------
DEFINES += EG_QWT      # Use QWT for convergence monitor?
DEFINES += EG_VTK      # Use VTK for postprocessing?
DEFINES += EG_PARAVIEW # Use ParaView for postprocessing?
DEFINES += EG_MATC     # Use MATC for internal operations in postprocessing?
DEFINES += EG_OCC      # Use OpenCASCADE 6.3 for importing CAD files? Needs VTK.
DEFINES -= EG_PYTHONQT # Use PythonQt for scripting in post processor?

#------------------------------------------------------------------------------
# 64 bit system?
#------------------------------------------------------------------------------
BITS = 32


#------------------------------------------------------------------------------
# Installation directory:
#------------------------------------------------------------------------------
ELMERGUI_HOME = $$(ELMERGUI_HOME)
isEmpty(ELMERGUI_HOME) {
   ELMER_HOME = $$(ELMER_HOME)
   isEmpty(ELMER_HOME) {
      unix: ELMER_HOME = /usr/local
      win32: ELMER_HOME = c:/Elmer7
      macx: ELMER_HOME = /usr/local
   }
   ELMERGUI_HOME = $${ELMER_HOME}/bin
}

#------------------------------------------------------------------------------
# Python library:
#------------------------------------------------------------------------------
unix {
   PY_INCLUDEPATH = /usr/include/python2.5
   PY_LIBPATH = /usr/lib
   PY_LIBS = -lpython2.5
}

win32 {
   PY_INCLUDEPATH = c:/PYTHON/Python-2.6.1/Include
   PY_LIBPATH = c:/PYTHON/Python-2.6.1/PCbuild
   PY_LIBS = -lpython26
}

macx {
}

#------------------------------------------------------------------------------
# QWT library:
#------------------------------------------------------------------------------
unix {
  QWT_INCLUDEPATH = /usr/include/qwt-qt4
  QWT_LIBPATH = /usr/lib
  QWT_LIBS = -lqwt-qt4
}

win32 {
  QWT_INCLUDEPATH = c:/Source/Qwt/include
  QWT_LIBPATH = c:/Source/Qwt/lib
  QWT_LIBS = -lqwt5
}

macx {
  QWT_INCLUDEPATH = /usr/local/qwt-5.0.2/include
  QWT_LIBPATH = /usr/local/qwt-5.0.2/lib
  QWT_LIBS =  -lqwt5
}

#------------------------------------------------------------------------------
# VTK library:
#------------------------------------------------------------------------------
unix {
   VTK_INCLUDEPATH = /usr/include/vtk-5.2
   VTK_LIBPATH = /usr/lib
   VTK_LIBS = -lQVTK \
              -lvtkCommon \
              -lvtkDICOMParser \
              -lvtkFiltering \
              -lvtkGenericFiltering \
              -lvtkGraphics \
              -lvtkHybrid \
              -lvtkIO \
              -lvtkImaging \
              -lvtkInfovis \
              -lvtkNetCDF \
              -lvtkRendering \
              -lvtkViews \
              -lvtkVolumeRendering \
              -lvtkWidgets
}

win32 {
   VTK_INCLUDEPATH = c:/Source/VTK/include/vtk-5.4
   VTK_LIBPATH = c:/Source/VTK/lib/vtk-5.4
   VTK_LIBS = -lQVTK \
              -lvtkCommon \
              -lvtkDICOMParser \
              -lvtkFiltering \
              -lvtkGenericFiltering \
              -lvtkGraphics \
              -lvtkHybrid \
              -lvtkIO \
              -lvtkImaging \
              -lvtkInfovis \
              -lvtkNetCDF \
              -lvtkRendering \
              -lvtkViews \
              -lvtkVolumeRendering \
              -lvtkWidgets \
              -lvtkexoIIc \
              -lvtkexpat \
              -lvtkfreetype \
              -lvtkftgl \
              -lvtkjpeg \
              -lvtklibxml2 \
              -lvtkmetaio \
              -lvtkpng \
              -lvtksys \
              -lvtktiff \
              -lvtkverdict \
              -lvtkzlib \
              -ladvapi32
}

macx {
   VTK_INCLUDEPATH = /usr/local/include/vtk-5.0
   VTK_LIBPATH = /usr/lib
   VTK_LIBS = -lvtkHybrid \
              -lvtkWidgets \
	      -lQVTK
}

#------------------------------------------------------------------------------
# OpenCASCADE library:
#------------------------------------------------------------------------------
unix {
   OCC_INCLUDEPATH = /usr/include/opencascade
   OCC_LIBPATH = /usr/lib
   OCC_LIBS = -lTKSTL \
              -lTKBRep \
              -lTKernel \
              -lTKG2d \
              -lTKG3d \
              -lTKGeomAlgo \
              -lTKGeomBase \
              -lTKMath \
              -lTKMesh \
              -lTKShHealing \
              -lTKSTEP \
              -lTKSTEP209 \
              -lTKSTEPAttr \
              -lTKSTEPBase \
              -lTKIGES \
              -lTKTopAlgo \
              -lTKXSBase
}

win32 {
   OCC_INCLUDEPATH = $(CASROOT)/inc
   OCC_LIBPATH = $(CASROOT)/win32/lib
   OCC_LIBS = $(CASROOT)/win32/lib/TKBRep.lib \
              $(CASROOT)/win32/lib/TKernel.lib \
              $(CASROOT)/win32/lib/TKG2d.lib \
              $(CASROOT)/win32/lib/TKG3d.lib \
              $(CASROOT)/win32/lib/TKGeomAlgo.lib \
              $(CASROOT)/win32/lib/TKGeomBase.lib \
              $(CASROOT)/win32/lib/TKMath.lib \
              $(CASROOT)/win32/lib/TKMesh.lib \
              $(CASROOT)/win32/lib/TKShHealing.lib \
              $(CASROOT)/win32/lib/TKSTEP.lib \
              $(CASROOT)/win32/lib/TKSTEP209.lib \
              $(CASROOT)/win32/lib/TKSTEPAttr.lib \
              $(CASROOT)/win32/lib/TKSTEPBase.lib \
              $(CASROOT)/win32/lib/TKIGES.lib \
              $(CASROOT)/win32/lib/TKTopAlgo.lib \
              $(CASROOT)/win32/lib/TKXSBase.lib
}

macx {
   OCC_INCLUDEPATH = 
   OCC_LIBPATH = 
   OCC_LIBS = 
}
