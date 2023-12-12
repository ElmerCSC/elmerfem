@echo off
setlocal
set PATH=C:\tools\msys64\mingw64\bin;%PATH%

if "%~1" == "clean" (
    echo Cleaning build directory...
    rmdir /s /q build
    exit /b 0
)

echo Creating build directory...
if not exist build mkdir build

cd build
@REM cmake -G "MSYS Makefiles" ..
@REM cmake -G "MSYS Makefiles" .. -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON
@REM cmake --build .
C:\tools\msys64\usr\bin\cmake.exe .. -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -Wno-dev
@REM "C:\Program Files\CMake\bin\cmake.exe" .. -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON
C:\tools\msys64\usr\bin\cmake.exe --build . -- VERBOSE=1
@REM C:\tools\msys64\usr\bin\cmake.exe --build . -- -j -- VERBOSE=1
@REM C:\tools\msys64\usr\bin\cmake.exe --build . -- -j -DCMAKE_C_FLAGS="-O3" -DCMAKE_CXX_FLAGS="-O3" -DCMAKE_Fortran_FLAGS="-O3"
cd ..

endlocal
exit /b 0
