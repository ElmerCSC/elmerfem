; elmer.nsi
; 
; Elmer installler generating script
;
!define ALL_USERS

!include WriteEnvStr.nsh 
!include AddToPath.nsh 

; The name of the installer
Name "Elmer 5.0"

; The file to write
OutFile "Elmer.exe"

; Icon c:\msys\1.0\home\admin\buildtools\elmer.ico

; The default installation directory
InstallDir c:\elmer

; The text to prompt the user to enter a directory
DirText "This will install the Elmer on your computer. Choose a directory. Hit next to begin installation."

; The stuff to install
Section "Main"
  ; Set output path to the installation directory.
  SetOutPath $INSTDIR
  ; Put file there
  File /r "c:\msys\1.0\home\admin\win32\*"
  SetOutPath "$INSTDIR\bin"
  ;  CreateShortCut "$DESKTOP\ElmerPost.lnk" "$INSTDIR\bin\elmerpost.bat" "" "$INSTDIR\bin\elmerpost.ico" 0
  SetOutPath "$INSTDIR"
  ; Write the installation path into the registry
  WriteRegStr HKLM SOFTWARE\Elmer5x "Install_Dir" "$INSTDIR"
  ; Write the uninstall keys for Windows
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Elmer5x" "DisplayName" "Elmer (remove only)"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Elmer5x" "UninstallString" '"$INSTDIR\uninstall.exe"'

  ; Global environment variables for all users across boots
  Push ELMER_HOME
  Push "$INSTDIR"
  Call WriteEnvStr

  Push ELMER_LIB
  Push "$INSTDIR\share\elmersolver\lib"
  Call WriteEnvStr

  Push ELMER_POST_HOME
  Push "$INSTDIR\share\elmerpost"
  Call WriteEnvStr


  ; PATH
  Push "$INSTDIR\share\elmersolver\lib"
  Call AddToPath

  Push "$INSTDIR\bin"
  Call AddToPath

  WriteUninstaller "$INSTDIR\Uninstall.exe"

SectionEnd ; end the section

;  uninstall section.
Section "Uninstall"
  ; remove registry keys
  DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Elmer5x"
  DeleteRegKey HKLM SOFTWARE\Elmer5x
  ; remove files
  ; MUST REMOVE UNINSTALLER, too
  Delete $INSTDIR\uninstall.exe
  RMDir /r "$INSTDIR"
  ; remove shortcuts, if any.
  ; remove directories used.
  Delete "$DESKTOP\ElmerPost.lnk"
  Delete "$DESKTOP\ElmerPost.pif"


  ; remove the env variables
  Push ELMER_HOME
  Call un.DeleteEnvStr

  Push ELMER_LIB
  Call un.DeleteEnvStr

  Push ELMER_POST_HOME
  Call un.DeleteEnvStr

  ; remove from path
  Push "$INSTDIR\bin"
  Call un.RemoveFromPath

  Push "$INSTDIR\share\elmersolver\lib"
  Call un.RemoveFromPath


SectionEnd

; eof
