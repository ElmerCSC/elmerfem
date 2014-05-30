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
  SET(CPACK_PACKAGE_VENDOR "TODO: CSC")
  SET(CPACK_PACAKGE_VERSION "${ELMER_FEM_MAJOR_VERSION}.${ELMER_FEM_MINOR_VERSION}")
  SET(CPACK_PACKAGE_CONTACT "TODO:elmeradm@csc.fi")

  #SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/Copyright.txt")
  #SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/Copyright.txt")
  SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/licenses/LICENSES_GPL")

  IF(CMAKE_SYSTEM_NAME MATCHES "Linux")
    SET(CPACK_GENERATOR "TGZ")
    IF(MAKE_DEB_PACKAGE)
      SET(CPACK_GENERATOR "${CPACK_GENERATOR};DEB")
    ENDIF()
    IF(MAKE_RPM_PACKAGE)  # TODO: untested
      SET(CPACK_GENERATOR "${CPACK_GENERATOR};RPM")
    ENDIF()
  ENDIF()

  IF(CMAKE_SYSTEM_NAME MATCHES "Windows")
    SET(CPACK_GENERATOR "NSIS;ZIP")
  ENDIF()

  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/cpack/NSISCPack.cmake)
  SET(CPACK_INSTALL_CMAKE_PROJECTS "${CMAKE_BINARY_DIR}" "Elmer" "ALL" "/")
  INCLUDE(CPack)
ENDIF()
