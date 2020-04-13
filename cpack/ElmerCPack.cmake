# Use CPack only if its cmake script exists
if(NOT EXISTS "${CMAKE_ROOT}/Modules/CPack.cmake")
  message(WARNING "${CMAKE_ROOT}/Modules/CPack.cmake does not exist")
  return()
endif()

set(CPACK_PACKAGE_NAME "Elmer")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "Open Source Finite Element Software for Multiphysical Problems")

set(CPACK_PACKAGE_DESCRIPTION
    "Elmer is an open source multiphysical
simulation software mainly developed by CSC - IT Center for Science (CSC).
Elmer development was started 1995 in collaboration with Finnish
Universities, research institutes and industry. After it's open source
publication in 2005, the use and development of Elmer has become
international.

Elmer includes physical models of fluid dynamics, structural mechanics,
electromagnetics, heat transfer and acoustics, for example. These are
described by partial differential equations which Elmer solves by the Finite
Element Method (FEM).")

set(CPACK_PACKAGE_VERSION_MAJOR "${ELMER_FEM_MAJOR_VERSION}")
set(CPACK_PACKAGE_VERSION_MINOR "${ELMER_FEM_MINOR_VERSION}")
set(CPACK_PACKAGE_VERSION_PATCH "${ELMER_FEM_REVISION}")

# SET(CPACK_PACKAGE_FILE_NAME "elmerfem-${ELMER_FEM_MAJOR_VERSION}.${ELMER_FEM_M
# INOR_VERSION}_${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
if(${CMAKE_VERSION} VERSION_GREATER 2.8.11)
  string(TIMESTAMP DATE "%Y%m%d")
else()
  message(
    WARNING "cmake ${CMAKE_VERSION} does not support STRING(TIMESTAMP ...)")
endif()

set(CPACK_PACKAGE_BASE_FILE_NAME
    "elmerfem"
    CACHE STRING "")
mark_as_advanced(CPACK_PACKAGE_BASE_FILE_NAME)
set(CPACK_PACKAGE_VENDOR "CSC")
set(CPACK_PACKAGE_VERSION
    "${ELMER_FEM_MAJOR_VERSION}.${ELMER_FEM_MINOR_VERSION}-${CPACK_PACKAGE_VERSION_PATCH}"
)
set(CPACK_PACKAGE_CONTACT "elmeradm@csc.fi")
if(CPACK_PACKAGE_FILE_NAME STREQUAL "")
  set(CPACK_PACKAGE_FILE_NAME
      "${CPACK_PACKAGE_BASE_FILE_NAME}-${CPACK_PACKAGE_VERSION}-${DATE}_${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}"
      CACHE STRING "" FORCE)
else()
  set(CPACK_PACKAGE_FILE_NAME
      "${CPACK_PACKAGE_BASE_FILE_NAME}-${CPACK_PACKAGE_VERSION}-${DATE}_${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}"
      CACHE STRING "")
endif(CPACK_PACKAGE_FILE_NAME STREQUAL "")

set(CPACK_RESOURCE_FILE_LICENSE
    "${CMAKE_CURRENT_SOURCE_DIR}/license_texts/LICENSES_GPL.txt")

message(STATUS "------------------------------------------------")
message(STATUS "  Package filename: ${CPACK_PACKAGE_FILE_NAME} ")
message(STATUS "  Patch version: ${CPACK_PACKAGE_VERSION} ")

if(NOT (BYPASS_DEB_DEPENDENCIES))
  set(CPACK_DEBIAN_PACKAGE_DEPENDS "libblas-dev, liblapack-dev")

  macro(ADD_DEBIAN_DEPENDENCY WITH_RULE DEPS)
    if(${WITH_RULE})
      list(APPEND DEP_LIST ${DEPS})
    endif(${WITH_RULE})
  endmacro()

  add_debian_dependency(WITH_MPI "openmpi-bin")
  add_debian_dependency(WITH_Mumps "libmumps-4.10.0")
  add_debian_dependency(WITH_Hypre "libhypre-2.8.0b")
  add_debian_dependency(WITH_ELMERGUI "libqt4-opengl")
  add_debian_dependency(WITH_ELMERGUILOGGER "libqt4-core")
  add_debian_dependency(WITH_ELMERGUITESTER "libqt4-core")
  add_debian_dependency(WITH_OCC "liboce-foundation" "liboce-modeling8")
  add_debian_dependency(WITH_PARAVIEW "paraview")
  add_debian_dependency(WITH_VTK "libvtk5.8-qt4" "libvtk5.8")
  add_debian_dependency(WITH_QWT "libqwt6")

  foreach(arg ${DEP_LIST})
    set(CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, ${arg}")
  endforeach()
endif()

if(CMAKE_SYSTEM_NAME MATCHES "Linux")
  mark_as_advanced(MAKE_DEB_PACKAGE MAKE_RPM_PACKAGE MAKE_TGZ_PACKAGE)
  # MESSAGE(STATUS "DEB package dependencies ${CPACK_DEBIAN_PACKAGE_DEPENDS}")
  set(MAKE_DEB_PACKAGE
      TRUE
      CACHE BOOL "Create DEB package with cpack")
  set(MAKE_RPM_PACKAGE
      TRUE
      CACHE BOOL "Create RPM package with cpack")
  set(MAKE_TGZ_PACKAGE
      TRUE
      CACHE BOOL "Create TGZ package with cpack")
  if(MAKE_TGZ_PACKAGE)
    list(APPEND CPACK_GENERATOR TGZ)
  endif()
  if(MAKE_DEB_PACKAGE)
    list(APPEND CPACK_GENERATOR DEB)
  endif()
  if(MAKE_RPM_PACKAGE) # @TODO: untested
    set(CPACK_GENERATOR "${CPACK_GENERATOR};RPM")
  endif()
endif()

if(CMAKE_SYSTEM_NAME MATCHES "Windows")
  mark_as_advanced(MAKE_NSIS_PACKAGE MAKE_ZIP_PACKAGE
                   CPACK_BUNDLE_EXTRA_WINDOWS_DLLS)
  set(MAKE_ZIP_PACKAGE
      TRUE
      CACHE BOOL "Create windows .zip file")
  set(MAKE_NSIS_PACKAGE
      TRUE
      CACHE BOOL "Create windows installer executable")
  set(CPACK_BUNDLE_EXTRA_WINDOWS_DLLS
      TRUE
      CACHE BOOL "Bundle dlls in windows install.")

  if(CPACK_BUNDLE_EXTRA_WINDOWS_DLLS)
    install(FILES ${LAPACK_LIBRARIES} DESTINATION "bin")
    if(NOT (LAPACK_LIB))
      find_file(LAPACK_LIB liblapack.dll PATH_SUFFIXES "bin")
    endif()
    if(NOT (BLAS_LIB))
      find_file(BLAS_LIB libblas.dll PATH_SUFFIXES "bin")
    endif()

    # mingw runtime dynamic link libraries
    find_file(QUADMATH_LIB libquadmath-0.dll)
    find_file(WINPTHREAD_LIB libwinpthread-1.dll)
    find_file(STDCPP_LIB libstdc++-6.dll)
    # if 1
    find_file(MINGW_GFORT_LIB libgfortran-5.dll)
    find_file(GCC_LIB libgcc_s_seh-1.dll)
    find_file(DBL_LIB libdouble-conversion.dll)
    find_file(GMP_LIB libgmp-10.dll)
    find_file(Z1_LIB zlib1.dll)
    install(
      FILES ${MINGW_GFORT_LIB}
            ${QUADMATH_LIB}
            ${WINPTHREAD_LIB}
            ${GCC_LIB}
            ${STDCPP_LIB}
            ${BLAS_LIB}
            ${LAPACK_LIB}
            ${DBL_LIB}
            ${GMP_LIB}
            ${Z1_LIB}
      DESTINATION "bin")
    # else
    find_file(MINGW_GFORT_LIB libgfortran-3.dll)
    find_file(GCC_LIB libgcc_s_sjlj-1.dll)
    install(
      FILES ${MINGW_GFORT_LIB}
            ${QUADMATH_LIB}
            ${WINPTHREAD_LIB}
            ${GCC_LIB}
            ${STDCPP_LIB}
            ${BLAS_LIB}
            ${LAPACK_LIB}
      DESTINATION "bin")
    # endif

    # Here we augment the installation by some needed dll's that should be
    # included with QT5. This is a quick and dirty remedy. I'm sure there is a
    # prettier way too.
    if(WITH_QT5)
      find_file(QTF0 tbb.dll)
      find_file(QTF1 libbz2-1.dll)
      find_file(QTF2 libfreetype-6.dll)
      find_file(QTF3 libglib-2.0-0.dll)
      find_file(QTF4 libgraphite2.dll)
      find_file(QTF5 libharfbuzz-0.dll)
      find_file(QTF6 libiconv-2.dll)
      find_file(QTF7 libicudt65.dll)
      find_file(QTF8 libicuin65.dll)
      find_file(QTF9 libicuuc65.dll)
      find_file(QTF10 libintl-8.dll)
      find_file(QTF11 libpcre-1.dll)
      find_file(QTF12 libpcre2-16-0.dll)
      find_file(QTF13 libpng16-16.dll)
      find_file(QTF14 libzstd.dll)
      install(
        FILES ${QTF0}
              ${QTF1}
              ${QTF2}
              ${QTF3}
              ${QTF4}
              ${QTF5}
              ${QTF6}
              ${QTF7}
              ${QTF8}
              ${QTF9}
              ${QTF10}
              ${QTF11}
              ${QTF12}
              ${QTF13}
              ${QTF14}
        DESTINATION "bin")
      install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../platforms"
              DESTINATION "bin")
    endif()

    if(BUNDLE_STRIPPED_GFORTRAN)
      # TODO: This will make the windows package to be GPL3
      install(
        DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../stripped_gfortran"
        DESTINATION "."
        COMPONENT "stripped_gfortran")
      set(CPACK_COMPONENT_STRIPPED_GFORTRAN_DESCRIPTION
          "A stripped version of x86_64-w64-mingw32-gfortran 9.2.0 compiler for compiling Elmer modules."
      )
      set(CPACK_COMPONENT_STRIPPED_GFORTRAN_DISPLAY_NAME "gfortran 9.2.0")
    endif()

    if(WITH_MPI)
      if(BUNDLE_MSMPI_REDIST)
        install(
          FILES "${CMAKE_CURRENT_SOURCE_DIR}/../msmpi_redist/msmpisetup.exe"
          DESTINATION "redist"
          COMPONENT "MS_MPI_Redistributable")
        # these are for Microsoft C++ (not needed for gcc/gfortan) INSTALL(FILES
        # "${CMAKE_CURRENT_SOURCE_DIR}/../msmpi_redist/vcredist_x64.exe"
        # DESTINATION "redist" COMPONENT "MS_MPI_Redistributable") INSTALL(FILES
        # "${CMAKE_CURRENT_SOURCE_DIR}/../msmpi_redist/vcredist_x86.exe"
        # DESTINATION "redist" COMPONENT "MS_MPI_Redistributable")
        set(CPACK_COMPONENT_MS_MPI_REDISTRIBUTABLE_DESCRIPTION
            "Install MS-MPI 10.1.1. Redistributable Package")
        set(CPACK_COMPONENT_MS_MPI_REDISTRIBUTABLE_DISPLAY_NAME "MS-MPI")
        list(
          APPEND
          CPACK_NSIS_EXTRA_INSTALL_COMMANDS
          "
        IfFileExists '$INSTDIR\\\\redist\\\\msmpisetup.exe' MSMpiSetupExists MsMpiSetupNotExist
        MsMpiSetupExists:
#        ExecWait '$INSTDIR\\\\redist\\\\vcredist_x64.exe'
#        ExecWait '$INSTDIR\\\\redist\\\\vcredist_x86.exe'
        ExecWait '$INSTDIR\\\\redist\\\\msmpisetup.exe'
        MsMpiSetupNotExist:
        ")
      endif()
    endif()
  endif()

  if(MAKE_NSIS_PACKAGE)
    set(CPACK_GENERATOR "NSIS")
  endif()
  if(MAKE_ZIP_PACKAGE)
    set(CPACK_GENERATOR "${CPACK_GENERATOR};ZIP")
  endif()

  if(MAKE_NSIS_PACKAGE)
    include(${CMAKE_CURRENT_SOURCE_DIR}/cpack/NSISCPack.cmake)
  endif()
endif()

set(CPACK_INSTALL_CMAKE_PROJECTS "${CMAKE_BINARY_DIR}" "Elmer" "ALL" "/")

if(WITH_ELMERGUI)
  set(CPACK_PACKAGE_EXECUTABLES "ElmerGUI" "ElmerGUI")
  set(CPACK_CREATE_DESKTOP_LINKS "ElmerGUI")
endif(WITH_ELMERGUI)

if(WITH_ELMERPOST)
  set(CPACK_PACKAGE_EXECUTABLES ${CPACK_PACKAGE_EXECUTABLES} "ElmerPost"
                                "ElmerPost")
endif()

if(WITH_ELMERGUITESTER)
  set(CPACK_PACKAGE_EXECUTABLES ${CPACK_PACKAGE_EXECUTABLES} "ElmerGUItester"
                                "ElmerGUItester")
endif(WITH_ELMERGUITESTER)

if(WITH_ELMERGUILOGGER)
  set(CPACK_PACKAGE_EXECUTABLES ${CPACK_PACKAGE_EXECUTABLES} "ElmerGUIlogger"
                                "ElmerGUIlogger")
endif(WITH_ELMERGUILOGGER)

include(CPack)
