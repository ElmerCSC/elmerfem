@echo off
Rem Elmer Front call 

Rem Martti Verho, 21.09.99

set ELMER_HOME=\elmer\dist
set ELMER_HOME_BS=\elmer\dist

set ELMER_FRONT_HOME=%ELMER_HOME%/Front

set PATH=%ELMER_HOME_BS%\bin;%ELMER_HOME_BS%\Front\bin;%PATH%

set TK_LIBRARY=%ELMER_HOME%/lib/tk8.3
set TCL_LIBRARY=%ELMER_HOME%/lib/tcl8.3

Rem NOTE: start ... calling needed!!!
Rem Otherwise calling other Elmer module batch fails!!
start %ELMER_HOME_BS%\bin\Front.exe %*
rem exit
