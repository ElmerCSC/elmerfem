CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
SET(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
MESSAGE(STATUS "------------------------------------------------")


IF(WITH_QT5)
  MESSAGE(STATUS "------------------------------------------------")
  SET(QT5_PKG_LIST Qt5OpenGL Qt5Xml Qt5Script Qt5Gui Qt5Core Qt5PrintSupport)
  FOREACH(_pkg ${QT5_PKG_LIST})
    FIND_PACKAGE(${_pkg} PATHS ${QT5_PATH})
  ENDFOREACH()
  ADD_DEFINITIONS(-DWITH_QT5)
  MESSAGE(STATUS "  [ElmerGUIlogger] Qt5:               " ${Qt5_FOUND})
  MESSAGE(STATUS "  [ElmerGUIlogger] Qt5 Libraries: ${Qt5OpenGL_LIBRARIES} ${Qt5Xml_LIBRARIES} ${Qt5Script_LIBRARIES} ${Qt5Gui_LIBRARIES} ${Qt5Core_LIBRARIES}")
  MESSAGE(STATUS "------------------------------------------------")
ELSE()  
  MESSAGE(STATUS "------------------------------------------------")
  FIND_PACKAGE(Qt4 COMPONENTS QtOpenGL QtXml QtScript QtGui QtCore REQUIRED)
  MESSAGE(STATUS "  [ElmerGUIlogger] Qt4:               " ${Qt4_FOUND})
  MESSAGE(STATUS "  [ElmerGUIlogger] Qt4_LIBRARIES:     " ${Qt4_LIBRARY})
  MESSAGE(STATUS "------------------------------------------------")
ENDIF()


IF(ElmerGUIloggerSTANDALONE)
  PROJECT(ElmerGUIlogger CXX C)
  IF(WIN32)
    INCLUDE(cmake/windows_bundle.cmake)
  ENDIF(WIN32)
ENDIF(ElmerGUIloggerSTANDALONE)


SET(CMAKE_INCLUDE_CURRENT_DIR ON)
SET(CMAKE_AUTOMOC ON)
SET(CMAKE_AUTORCC ON)
SET(CMAKE_AUTOUIC ON)

SET(HEADERS src/mainwindow.h)
SET(SOURCES src/main.cpp src/mainwindow.cpp)
SET(RESOURCES ElmerGUIlogger.qrc)
SET(RC_FILE ElmerGUIlogger.rc)


IF(WITH_QT5)
  QT5_ADD_RESOURCES(UI_RESOURCES ElmerGUIlogger.qrc)
  ADD_DEFINITIONS(-DWITH_QT5)
ELSE()
  QT4_ADD_RESOURCES(UI_RESOURCES ElmerGUIlogger.qrc)
  MESSAGE(WARNING "QT_USE_FILE: ${QT_USE_FILE}")
  INCLUDE(${QT_USE_FILE})
  ADD_DEFINITIONS(${QT_DEFINITIONS})
ENDIF()



IF(MINGW AND WIN32)
  SET(CMAKE_RC_COMPILER_INIT windres)
  ENABLE_LANGUAGE(RC)
  SET(CMAKE_RC_COMPILE_OBJECT
    "<CMAKE_RC_COMPILER> <FLAGS> -O coff <DEFINES> -i <SOURCE> -o <OBJECT>")
  SET(SOURCES ${SOURCES} ElmerGUIlogger.rc)
ENDIF()


ADD_EXECUTABLE(ElmerGUIlogger WIN32 ${SOURCES} ${UI_RESOURCES})
TARGET_LINK_LIBRARIES(ElmerGUIlogger ${QT_LIBRARIES})

INSTALL(TARGETS ElmerGUIlogger RUNTIME DESTINATION "bin" COMPONENT "elmergui")
