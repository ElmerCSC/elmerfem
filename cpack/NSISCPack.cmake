# CPack configuration for NSIS installer

IF(WIN32)
  SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} Elmer")
  SET(CPACK_NSIS_HELP_LINK "TODO: http://www.elmerfem.org")
  SET(CPACK_NSIS_CONTACT "TODO: elmeradm@csc.fi")
  SET(CPACK_NSIS_EXTRA_INSTALL_COMMANDS
"   !include \\\"winmessages.nsh\\\"
   ; HKLM (all users) vs HKCU (current user) defines
   !define env_hklm 'HKLM \\\"SYSTEM\\\\CurrentControlSet\\\\Control\\\\Session Manager\\\\Environment\\\"'
   !define env_hkcu 'HKCU \\\"Environment\\\"'
   StrCmp \\\$ADD_TO_PATH_ALL_USERS \\\"1\\\" WriteAllElmerHomeKey
     DetailPrint \\\"Selected environment for current user\\\"
     WriteRegExpandStr \\\${env_hkcu} ELMER_HOME \\\$INSTDIR
     WriteRegExpandStr \\\${env_hkcu} ELMERGUI_HOME \\\$INSTDIR
     Goto DoSendElmerHome
   WriteAllElmerHomeKey:
     DetailPrint \\\"Selected environment for all users\\\"
     WriteRegExpandStr \\\${env_hklm} ELMER_HOME \\\$INSTDIR
     WriteRegExpandStr \\\${env_hklm} ELMERGUI_HOME \\\$INSTDIR  
     DoSendElmerHome:
   SendMessage \\\${HWND_BROADCAST} \\\${WM_WININICHANGE} 0 \\\"STR:Environment\\\" /TIMEOUT=5000 ")
  SET(CPACK_NSIS_EXTRA_UNINSTALL_COMMANDS
"   ; delete variable
   StrCmp \\\${ADD_TO_PATH_ALL_USERS} \\\"1\\\" unWriteAllElmerHome
     DeleteRegValue \\\${env_hkcu} ELMER_HOME 
     DeleteRegValue \\\${env_hkcu} ELMERGUI_HOME 
     Goto unDoSendElmerHome
   unWriteAllElmerHome:
     DeleteRegValue \\\${env_hklm} ELMER_HOME
     DeleteRegValue \\\${env_hklm} ELMERGUI_HOME
   unDoSendElmerHome:
     SendMessage \\\${HWND_BROADCAST} \\\${WM_WININICHANGE} 0 \\\"STR:Environment\\\" /TIMEOUT=5000")
  SET(CPACK_NSIS_MODIFY_PATH "ON")
  # @TODO: This is build-bot specific hack to find runtime libraries
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

  SET(CPACK_COMPONENT_UNSPECIFIED_DISPLAY_NAME "Elmerfem solver")
  SET(CPACK_COMPONENT_UNSPECIFIED_DESCRIPTION "The main application: ElmerSolver, ElmerGrid, matc and runtime binaries.")

  SET(CPACK_NSIS_COMPONENT_INSTALL TRUE)
  SET(CPACK_COMPONENT_INSTALL_ALL "elmergui gfortran_minimal Unspecified")
ENDIF(WIN32)
