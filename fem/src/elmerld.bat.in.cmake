ECHO OFF
SET LIBDIR=%ELMER_HOME%/share/elmersolver/lib
SET INCLUDE=%LIBDIR%/../include

SET LD=%ELMER_HOME%/stripped_gfortran/bin/@NSIS_GFORTRAN_RUNTIME_FILENAME@

SET cmd=%LD% -shared %* -L%ELMER_HOME%/bin -L%LIBDIR% -lelmersolver
echo %cmd%
%LD% -shared %* -L%ELMER_HOME%/bin -L%LIBDIR% -lelmersolver
