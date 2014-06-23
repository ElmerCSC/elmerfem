# Use CPack only if it really found
IF(EXISTS "${CMAKE_ROOT}/Modules/CPack.cmake")

  SET(CPACK_PACKAGE_NAME "Elmer")
  SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Open Source Finite Element Software for Multiphysical Problems")
  
  SET(CPACK_PACKAGE_DESCRIPTION "Elmer is an open source multiphysical
  simulation software mainly developed by CSC - IT Center for Science (CSC).
  Elmer development was started 1995 in collaboration with Finnish
  Universities, research institutes and industry. After it's open source
  publication in 2005, the use and development of Elmer has become
  international.

  Elmer includes physical models of fluid dynamics, structural mechanics,
  electromagnetics, heat transfer and acoustics, for example. These are
  described by partial differential equations which Elmer solves by the Finite
  Element Method (FEM).")

  SET(CPACK_PACKAGE_VERSION_MAJOR "${ELMER_FEM_MAJOR_VERSION}")
  SET(CPACK_PACKAGE_VERSION_MINOR "${ELMER_FEM_MINOR_VERSION}")

  #SET(CPACK_PACKAGE_FILE_NAME "elmerfem-${ELMER_FEM_MAJOR_VERSION}.${ELMER_FEM_MINOR_VERSION}_${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
  STRING(TIMESTAMP DATE "%Y%m%d")
  SET(CPACK_PACKAGE_FILE_NAME "elmerfem-snapshot-${DATE}_${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
  SET(CPACK_PACKAGE_VENDOR "CSC")
  SET(CPACK_PACAKGE_VERSION "${ELMER_FEM_MAJOR_VERSION}.${ELMER_FEM_MINOR_VERSION}")
  SET(CPACK_PACKAGE_CONTACT "elmeradm@csc.fi")

  SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/license_texts/LICENSES_GPL") # @TODO: License for gfortran?

  IF(CMAKE_SYSTEM_NAME MATCHES "Linux")
    MARK_AS_ADVANCED(MAKE_DEB_PACKAGE MAKE_RPM_PACKAGE)
    SET(MAKE_DEB_PACKAGE TRUE CACHE BOOL "Create DEB package with cpack")
    SET(MAKE_RPM_PACKAGE TRUE CACHE BOOL "Create RPM package with cpack")
    SET(CPACK_GENERATOR "TGZ")
    IF(MAKE_DEB_PACKAGE)
      SET(CPACK_GENERATOR "${CPACK_GENERATOR};DEB")
    ENDIF()
    IF(MAKE_RPM_PACKAGE)  # @TODO: untested
      SET(CPACK_GENERATOR "${CPACK_GENERATOR};RPM")
    ENDIF()
  ENDIF()

  IF(CMAKE_SYSTEM_NAME MATCHES "Windows")
    MARK_AS_ADVANCED(MAKE_NSIS_PACKAGE MAKE_ZIP_PACKAGE CPACK_BUNDLE_EXTRA_WINDOWS_DLLS)
    SET(MAKE_ZIP_PACKAGE TRUE CACHE BOOL "Create windows .zip file")
    SET(MAKE_NSIS_PACKAGE TRUE CACHE BOOL "Create windows installer executable")
    SET(CPACK_BUNDLE_EXTRA_WINDOWS_DLLS TRUE CACHE BOOL "Bundle dlls in windows install.")

    IF(CPACK_BUNDLE_EXTRA_WINDOWS_DLLS)
      INSTALL(FILES ${LAPACK_LIBRARIES} DESTINATION "bin")
      IF(NOT(LAPACK_LIB))
        FIND_FILE(LAPACK_LIB liblapack.dll PATH_SUFFIXES "bin")
      ENDIF()
      IF(NOT(BLAS_LIB))
        FIND_FILE(BLAS_LIB libblas.dll PATH_SUFFIXES "bin")
      ENDIF()

      # mingw runtime dynamic link libraries
      FIND_FILE(MINGW_GFORT_LIB libgfortran-3.dll)
      FIND_FILE(QUADMATH_LIB libquadmath-0.dll)
      FIND_FILE(WINPTHREAD_LIB libwinpthread-1.dll)
      FIND_FILE(GCC_LIB libgcc_s_sjlj-1.dll)
      FIND_FILE(STDCPP_LIB libstdc++-6.dll)

      # Qt4 dynamic link libraries for ElmerGUI
      IF(WITH_ELMERGUI)
        FIND_FILE(QTCORE_DLL QtCore4.dll PATH_SUFFIXES "bin")
        FIND_FILE(QTGUI_DLL QtGui4.dll PATH_SUFFIXES "bin")
        FIND_FILE(QTOPENGL_DLL QtOpenGL4.dll PATH_SUFFIXES "bin")
        FIND_FILE(QTSCRIPT_DLL QtScript4.dll PATH_SUFFIXES "bin")
        FIND_FILE(QTXML_DLL QtXml4.dll PATH_SUFFIXES "bin")
        INSTALL(FILES ${QTCORE_DLL} ${QTGUI_DLL} ${QTOPENGL_DLL} ${QTSCRIPT_DLL} ${QTXML_DLL} DESTINATION "bin" COMPONENT "elmergui")
      ENDIF()

      INSTALL(FILES ${MINGW_GFORT_LIB} ${QUADMATH_LIB} ${WINPTHREAD_LIB} ${GCC_LIB} ${STDCPP_LIB} ${BLAS_LIB} ${LAPACK_LIB} DESTINATION "bin")

      IF(BUNDLE_STRIPPED_GFORTRAN)
        INSTALL(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../stripped_gfortran" DESTINATION "." COMPONENT "stripped_gfortran")
        SET(CPACK_COMPONENT_STRIPPED_GFORTRAN_DESCRIPTION "A stripped version of x86_64-w64-mingw32-gfortran 4.8.3 (sjlj) compiler for compiling Elmer modules.")
        SET(CPACK_COMPONENT_STRIPPED_GFORTRAN_DISPLAY_NAME "gfortran 4.8.3")
      ENDIF()
    ENDIF()

    IF(MAKE_NSIS_PACKAGE)
      SET(CPACK_GENERATOR "NSIS")
    ENDIF()
    IF(MAKE_ZIP_PACKAGE)
      SET(CPACK_GENERATOR "${CPACK_GENERATOR};ZIP")
    ENDIF()
  ENDIF()

  IF(MAKE_NSIS_PACKAGE)
    INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/cpack/NSISCPack.cmake)
  ENDIF()
  SET(CPACK_INSTALL_CMAKE_PROJECTS "${CMAKE_BINARY_DIR}" "Elmer" "ALL" "/")
  INCLUDE(CPack)
ENDIF()
