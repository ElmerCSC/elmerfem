==== Windows under msys2

 * Install MSYS from https://www.msys2.org/
 * Launch MSYS via the "MSYS2 MinGW x64" shortcut.
 * Download fresh MSYS package databases and upgrade installed packages by running `pacman -Syu` twice.
 * Install Elmer MSYS dependencies:
    ** `pacman -S --noconfirm --needed base-devel mingw-w64-x86_64-toolchain mingw64/mingw-w64-x86_64-cmake mingw64/mingw-w64-x86_64-openblas mingw64/mingw-w64-x86_64-qt5 mingw64/mingw-w64-x86_64-qwt-qt5 mingw64/mingw-w64-x86_64-nsis git`
 * Get the Elmer source code:
    ** `git clone https://github.com/ElmerCSC/elmerfem`
 * Create directories required for building a local Elmer install
    ** `mkdir -p bundle_msys2/bin bundle_qt5/bin platforms`
 * Create a build directory for build artifacts
    ** `mkdir -p build`
 * Run CMake to prepare the build with executable binaries in an "install" directory. Note that the QWT_INCLUDE_DIR needs to be correctly set to match the MSYS installation location.
    ** `cd build`
    ** `cmake -G "MSYS Makefiles" -DWITH_ELMERGUI:BOOL=TRUE -DWITH_MPI:BOOL=FALSE -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_Fortran_COMPILER=/mingw64/bin/gfortran.exe -DQWT_INCLUDE_DIR=C:/msys64/mingw64/include/qwt-qt5/ -DWIN32:BOOL=TRUE -DCPACK_BUNDLE_EXTRA_WINDOWS_DLLS:BOOL=TRUE ../elmerfem`
 * Build the source code and create a local installation
    ** `make install`
 * Copy additional dependencies from /mingw64/bin/ to the "install" directory:
    ** libgfortran-5.dll libgcc_s_seh-1.dll libopenblas.dll libquadmath-0.dll libwinpthread-1.dll libstdc++-6.dll qwt-qt5.dll libdouble-conversion.dll libicuin69.dll libicuuc69.dll libpcre2-16-0.dll libharfbuzz-0.dll libmd4c.dll libpng16-16.dll zlib1.dll libzstd.dll libicudt69.dll libfreetype-6.dll libglib-2.0-0.dll libgraphite2.dll libintl-8.dll libbz2-1.dll libbrotlidec.dll libpcre-1.dll libiconv-2.dll libbrotlicommon.dll
 * Copy Qt platform dependencies into the "install/bin/platforms" directory
    ** `cp /mingw64/share/qt5/plugins/platforms/qwindows.dll ../install/bin/platforms`
 * Binaries like ElmerSolver.exe or ElmerGUI.exe can now be run from the ../install/bin directory.





//=== Package managers

//[.text-center]
//image::https://repology.org/badge/vertical-allrepos/elmerfem.svg["Packaging status", link=https://repology.org/project/elmerfem/versions]

//==== Chocolatey

//[.text-center]
//image:https://img.shields.io/chocolatey/dt/elmer-mpi["Chocolatey", link=https://chocolatey.org/packages/elmer-mpi]

