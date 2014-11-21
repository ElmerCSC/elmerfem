@ECHO OFF
SET LIBDIR=%ELMER_HOME%/bin
SET INCLUDE=%ELMER_HOME%/share/elmersolver/include

SET FC="%ELMER_HOME%/stripped_gfortran/bin/@NSIS_GFORTRAN_RUNTIME_FILENAME@"

SET cmd=%FC% %* @CMAKE_Fortran_FLAGS@ @ELMER_F90FLAGS@ @CMAKE_SHARED_LIBRARY_Fortran_FLAGS@ @CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS@ -I"%INCLUDE%" -L"%LIBDIR%" -shared -lelmersolver
echo %cmd%
%cmd%
