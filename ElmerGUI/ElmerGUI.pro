#==============================================================================
#
#        ElmerGUI: qmake project file for Unix, Win32, and MacX
#
#     For more details, see the project files in the subdirectories
#
#==============================================================================

include(ElmerGUI.pri)

#------------------------------------------------------------------------------
# Test the configuration for some headers:
#------------------------------------------------------------------------------
contains(DEFINES,EG_QWT) {
   !exists($${QWT_INCLUDEPATH}/qwt.h) {
      message("EG_QWT has been defined, but qwt.h was not found")
      message("Check QWT_INCLUDEPATH or undefine EG_QWT in ElmerGUI.pri")
      error("Detected inconsistent configuration. Unable to continue.")
   }
}

contains(DEFINES,EG_VTK) {
   !exists($${VTK_INCLUDEPATH}/QVTKWidget.h) {
      message("EG_VTK has been defined, but QVTKWidget.h was not found")
      message("Check VTK_INCLUDEPATH or undefine EG_VTK in ElmerGUI.pri")
      error("Detected inconsistent configuration. Unable to continue.")
   }
}

contains(DEFINES,EG_OCC) {
   !exists($${OCC_INCLUDEPATH}/BRepTools.hxx) {
      message("EG_OCC has been defined, but BRepTools.hxx was not found")
      message("Check OCC_INCLUDEPATH or undefine EG_OCC in ElmerGUI.pri")
      error("Detected inconsistent configuration. Unable to continue.")
   }
}

contains(DEFINES,EG_PYTHONQT) {
   !exists($${PY_INCLUDEPATH}/Python.h) {
      message("EG_PYTHONQT has been defined, but Python.h was not found")
      message("Check PY_INCLUDEPATH or undefine EG_PYTHONQT in ElmerGUI.pri")
      error("Detected inconsistent configuration. Unable to continue.")
   }
}

message(ELMERGUI_HOME=$${ELMERGUI_HOME})

#------------------------------------------------------------------------------
# Build in all subdirectories:
#------------------------------------------------------------------------------
TEMPLATE = subdirs
SUBDIRS = matc netgen
contains(DEFINES, EG_PYTHONQT): SUBDIRS += PythonQt
SUBDIRS += Application
