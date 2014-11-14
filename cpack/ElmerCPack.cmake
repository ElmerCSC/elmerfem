# Use CPack only if its cmake script exists
IF(NOT EXISTS "${CMAKE_ROOT}/Modules/CPack.cmake")
  MESSAGE(WARNING "${CMAKE_ROOT}/Modules/CPack.cmake does not exist")
  RETURN()
ENDIF()

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
SET(CPACK_PACKAGE_VERSION_PATCH "${ELMER_FEM_REVISION}")

#SET(CPACK_PACKAGE_FILE_NAME "elmerfem-${ELMER_FEM_MAJOR_VERSION}.${ELMER_FEM_MINOR_VERSION}_${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
IF(${CMAKE_VERSION} VERSION_GREATER 2.8.11)
  STRING(TIMESTAMP DATE "%Y%m%d")
ELSE()
  MESSAGE(WARNING "cmake ${CMAKE_VERSION} does not support STRING(TIMESTAMP ...)")
ENDIF()

SET(CPACK_PACKAGE_BASE_FILE_NAME "elmerfem" CACHE STRING "")
MARK_AS_ADVANCED(CPACK_PACKAGE_BASE_FILE_NAME)
SET(CPACK_PACKAGE_VENDOR "CSC")
SET(CPACK_PACKAGE_VERSION "${ELMER_FEM_MAJOR_VERSION}.${ELMER_FEM_MINOR_VERSION}-${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_PACKAGE_CONTACT "elmeradm@csc.fi")
SET(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_BASE_FILE_NAME}-${CPACK_PACKAGE_VERSION}-${DATE}_${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}" CACHE STRING "")

SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/license_texts/LICENSES_GPL.txt") 

MESSAGE(STATUS "------------------------------------------------")
MESSAGE(STATUS "  Package filename: ${CPACK_PACKAGE_FILE_NAME} ")
MESSAGE(STATUS "  Patch version: ${CPACK_PACKAGE_VERSION} ")

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

    INSTALL(FILES ${MINGW_GFORT_LIB} ${QUADMATH_LIB} ${WINPTHREAD_LIB} ${GCC_LIB} ${STDCPP_LIB} ${BLAS_LIB} ${LAPACK_LIB} DESTINATION "bin")

    IF(BUNDLE_STRIPPED_GFORTRAN)
      # TODO: This will make the windows package to be GPL3
      INSTALL(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../stripped_gfortran" DESTINATION "." COMPONENT "stripped_gfortran")
      SET(CPACK_COMPONENT_STRIPPED_GFORTRAN_DESCRIPTION "A stripped version of x86_64-w64-mingw32-gfortran 4.8.3 (sjlj) compiler for compiling Elmer modules.")
      SET(CPACK_COMPONENT_STRIPPED_GFORTRAN_DISPLAY_NAME "gfortran 4.8.3")
    ENDIF()

    IF(WITH_MPI)
      IF(BUNDLE_MSMPI_REDIST)
        INSTALL(FILES "${CMAKE_CURRENT_SOURCE_DIR}/../msmpi_redist/MSMpiSetup.exe" DESTINATION "redist" COMPONENT "MS_MPI_Redistributable")
        INSTALL(FILES "${CMAKE_CURRENT_SOURCE_DIR}/../msmpi_redist/vcredist_x64.exe" DESTINATION "redist" COMPONENT "MS_MPI_Redistributable")
        INSTALL(FILES "${CMAKE_CURRENT_SOURCE_DIR}/../msmpi_redist/vcredist_x86.exe" DESTINATION "redist" COMPONENT "MS_MPI_Redistributable")
        SET(CPACK_COMPONENT_MS_MPI_REDISTRIBUTABLE_DESCRIPTION "Install HPC Pack 2012 R2 MS-MPI Redistributable Package")
        SET(CPACK_COMPONENT_MS_MPI_REDISTRIBUTABLE_DISPLAY_NAME "MS-MPI")
        LIST(APPEND CPACK_NSIS_EXTRA_INSTALL_COMMANDS "
        IfFileExists '$INSTDIR\\\\redist\\\\MSMpiSetup.exe' MSMpiSetupExists MsMpiSetupNotExist
        MsMpiSetupExists:
        ExecWait '$INSTDIR\\\\redist\\\\vcredist_x64.exe'
        ExecWait '$INSTDIR\\\\redist\\\\vcredist_x86.exe'
        ExecWait '$INSTDIR\\\\redist\\\\MSMpiSetup.exe'
        MsMpiSetupNotExist:
        ")
      ENDIF()
    ENDIF()
  ENDIF()

  IF(MAKE_NSIS_PACKAGE)
    SET(CPACK_GENERATOR "NSIS")
  ENDIF()
  IF(MAKE_ZIP_PACKAGE)
    SET(CPACK_GENERATOR "${CPACK_GENERATOR};ZIP")
  ENDIF()

  IF(MAKE_NSIS_PACKAGE)
    INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/cpack/NSISCPack.cmake)
  ENDIF()
ENDIF()


SET(CPACK_INSTALL_CMAKE_PROJECTS "${CMAKE_BINARY_DIR}" "Elmer" "ALL" "/")

IF(WITH_ELMERGUI)
  SET(CPACK_PACKAGE_EXECUTABLES "ElmerGUI" "ElmerGUI")
  SET(CPACK_CREATE_DESKTOP_LINKS "ElmerGUI")
ENDIF(WITH_ELMERGUI)

IF(WITH_ELMERGUITESTER)
  SET(CPACK_PACKAGE_EXECUTABLES ${CPACK_PACKAGE_EXECUTABLES} "ElmerGUItester" "ElmerGUItester")
ENDIF(WITH_ELMERGUITESTER)

INCLUDE(CPack)
