# CPack configuration for NSIS installer

IF(WIN32)
  SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} Elmer")
  SET(CPACK_NSIS_HELP_LINK "http://www.elmerfem.org")
  #SET(CPACK_NSIS_CONTACT "TODO: elmeradm@csc.fi")
  SET(CPACK_NSIS_CONTACT "")
  LIST(APPEND CPACK_NSIS_EXTRA_INSTALL_COMMANDS
"   !include \\\"winmessages.nsh\\\"
   ; HKLM (all users) vs HKCU (current user) defines
   !define env_hklm 'HKLM \\\"SYSTEM\\\\CurrentControlSet\\\\Control\\\\Session Manager\\\\Environment\\\"'
   !define env_hkcu 'HKCU \\\"Environment\\\"'
   StrCmp \\\$ADD_TO_PATH_ALL_USERS \\\"1\\\" WriteAllElmerHomeKey
     DetailPrint \\\"Selected environment for current user\\\"
     WriteRegExpandStr \\\${env_hkcu} ELMER_HOME \\\$INSTDIR
     WriteRegExpandStr \\\${env_hkcu} ELMERGUI_HOME \\\$INSTDIR\\\\share\\\\ElmerGUI
     Goto DoSendElmerHome
   WriteAllElmerHomeKey:
     DetailPrint \\\"Selected environment for all users\\\"
     WriteRegExpandStr \\\${env_hklm} ELMER_HOME \\\$INSTDIR
     WriteRegExpandStr \\\${env_hklm} ELMERGUI_HOME \\\$INSTDIR\\\\share\\\\ElmerGUI 
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

  SET(CPACK_COMPONENT_UNSPECIFIED_DISPLAY_NAME "Elmerfem solver")
  SET(CPACK_COMPONENT_UNSPECIFIED_DESCRIPTION "The main application: ElmerSolver, ElmerGrid, matc and runtime binaries.")
  SET(CPACK_COMPONENT_ELMERGUI_DISPLAY_NAME "ElmerGUI")
  SET(CPACK_COMPONENT_ELMERGUI_SAMPLES_DISPLAY_NAME "ElmerGUI samples")
  SET(CPACK_COMPONENT_ELMERGUI_SAMPLES_DESCRIPTION "Geometry samples for ElmerGUI")

  SET(CPACK_COMPONENT_ELMERPOST_DISPLAY_NAME "ElmerPost")
  SET(CPACK_COMPONENT_ELMERPOST_DESCRIPTION "A post processor for Elmer")

  SET(CPACK_NSIS_COMPONENT_INSTALL TRUE)
  SET(CPACK_COMPONENT_INSTALL_ALL "elmergui gfortran_minimal Unspecified elmergui_samples ElmerPost")
ENDIF(WIN32)
