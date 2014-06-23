find_package(PythonLibs 2.7 REQUIRED)
find_package(Qt4 COMPONENTS QtCore QtGui REQUIRED)

SET(CMAKE_AUTOMOC ON)
SET(CMAKE_INCLUDE_CURRENT_DIR ON)

SET(HEADERS src/PythonQt.h 
  src/PythonQtStdDecorators.h 
  src/PythonQtClassInfo.h 
  src/PythonQtImporter.h 
  src/PythonQtObjectPtr.h 
  src/PythonQtSlot.h 
  src/PythonQtStdOut.h 
  src/PythonQtMisc.h 
  src/PythonQtMethodInfo.h 
  src/PythonQtImportFileInterface.h 
  src/PythonQtConversion.h 
  src/PythonQtSignalReceiver.h 
  src/PythonQtWrapper.h 
  src/PythonQtMetaObjectWrapper.h 
  src/PythonQtCppWrapperFactory.h 
  src/PythonQtVariants.h 
  src/PythonQtVariantWrapper.h 
  src/wrapper/PythonQtWrappedVariants.h 
  src/gui/PythonQtScriptingConsole.h    
  src/PythonQtSystem.h)


SET(SOURCES 
  src/PythonQtStdDecorators.cpp 
  src/PythonQt.cpp 
  src/PythonQtClassInfo.cpp 
  src/PythonQtImporter.cpp 
  src/PythonQtObjectPtr.cpp 
  src/PythonQtStdOut.cpp 
  src/PythonQtSlot.cpp 
  src/PythonQtMisc.cpp 
  src/PythonQtMethodInfo.cpp 
  src/PythonQtConversion.cpp 
  src/PythonQtSignalReceiver.cpp 
  src/PythonQtVariants.cpp 
  src/PythonQtVariantWrapper.cpp 
  src/PythonQtWrapper.cpp 
  src/PythonQtMetaObjectWrapper.cpp 
  src/gui/PythonQtScriptingConsole.cpp)

INCLUDE(${QT_USE_FILE})
ADD_DEFINITIONS(${QT_DEFINITIONS})
ADD_LIBRARY(PythonQt ${SOURCES} ${HEADERS})
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/src)
TARGET_LINK_LIBRARIES(PythonQt ${PYTHON_LIBRARIES} ${QT_LIBRARIES})
MESSAGE(STATUS "PYTHON_LIBRARIES ${PYTHON_LIBRARIES}")
INSTALL(TARGETS PythonQt RUNTIME DESTINATION "bin" LIBRARY DESTINATION "lib"
	ARCHIVE DESTINATION "lib" COMPONENT "elmergui")

