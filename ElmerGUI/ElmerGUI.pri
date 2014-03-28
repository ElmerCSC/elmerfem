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
DEFINES += EG_OCC      # Use OpenCASCADE 6.5 for importing CAD files? Needs VTK.
DEFINES -= EG_PYTHONQT # Use PythonQt for scripting in post processor?

#------------------------------------------------------------------------------
# 64 bit system?
#------------------------------------------------------------------------------
BITS = 64


#------------------------------------------------------------------------------
# Installation directory:
#------------------------------------------------------------------------------
#ELMERGUI_HOME = /home/ssillanp/elmer/ElmerGUI/Installation/ElmerGUI/bin
#ELMERGUI_HOME = $$(ELMERGUI_HOME)
ELMERGUI_HOME = /opt/elmer/bin
isEmpty(ELMERGUI_HOME) {
   ELMER_HOME = $$(ELMER_HOME)
   isEmpty(ELMER_HOME) {
      unix: ELMER_HOME = /opt/elmer
      win32: ELMER_HOME = c:/Elmer7
      macx: ELMER_HOME = /usr/local
   }
#  ELMERGUI_HOME = $${ELMER_HOME}/bin
}


#------------------------------------------------------------------------------
# Paraview directory
#------------------------------------------------------------------------------
PARAVIEW_HOME = /opt/ParaView
isEmpty(PARAVIEW_HOME) {
   PARAVIEW_HOME = $$(PARAVIEW_HOME)
   isEmpty(PARAVIEW_HOME) {
      unix: PARAVIEW_HOME = /opt/ParaView
      win32: PARAVIEW_HOME = c:/ParaView
      macx: PARAVIEW_HOME = /usr/local/ParaView
   }
}


#------------------------------------------------------------------------------
# Python library:
#------------------------------------------------------------------------------
unix {
   PY_INCLUDEPATH = /home/ssillanp/elmer/ElmerGUI/Installation/Python-2.7.3/include/python2.7
   PY_LIBPATH = /home/ssillanp/elmer/ElmerGUI/Installation/Python-2.7.3/lib
   PY_LIBS = -lpython2.7
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
  QWT_INCLUDEPATH = /home/ssillanp/elmer/ElmerGUI/Installation/qwt/include
  QWT_LIBPATH = /home/ssillanp/elmer/ElmerGUI/Installation/qwt/lib
  QWT_LIBS = -lqwt
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
   VTK_INCLUDEPATH = /home/ssillanp/elmer/ElmerGUI/Installation/vtk/include/vtk-5.10
   VTK_LIBPATH = /home/ssillanp/elmer/ElmerGUI/Installation/vtk/lib/vtk-5.10
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
   OCC_INCLUDEPATH = /home/ssillanp/elmer/ElmerGUI/Installation/occt/inc
   OCC_LIBPATH = /home/ssillanp/elmer/ElmerGUI/Installation/occt/lib
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
