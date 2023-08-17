@echo off
setlocal

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
C:\tools\msys64\usr\bin\cmake.exe .. -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON
@REM C:\tools\msys64\usr\bin\cmake.exe --build .
cd ..

endlocal
exit /b 0
