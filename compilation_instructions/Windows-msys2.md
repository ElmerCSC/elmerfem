Build Elmer on Windows using MSYS2
==================================

* If you'd like to build with MPI, install the latest version of MSMPI:
  https://learn.microsoft.com/en-us/message-passing-interface/microsoft-mpi-release-notes
* Install MSYS2 from https://www.msys2.org/
* Open the MSYS2 `bash` via the "MSYS2 UCRT64" (or "MSYS2 MINGW64") shortcut.
* Download the latest MSYS2 package databases and upgrade installed packages by
  running `pacman -Syu` twice.
* Install dependencies of Elmer from MSYS2:
```
pacman -S --needed base-devel git pactoys
pacboy -S --needed toolchain:p cc:p fc:p cmake:p openblas:p parmetis:p mumps:p suitesparse:p qt6-declarative:p qwt-qt6:p opencascade:p
```
* If you'd like to build with VTK, also install the following packages:
```
pacboy -S --needed vtk:p adios2:p anari-sdk:p boost:p cgns:p cli11:p eigen3:p fast_float:p ffmpeg:p gl2ps:p liblas:p libmariadbclient:p netcdf:p openslide:p openvdb:p scnlib:p unixodbc:p utf8cpp:p
```
* If you'd like to build with MPI, also install the following package and set
the environment variable `MSMPI_BIN` to the path where you installed MSMPI:
```
pacboy -S --needed msmpi:p
export MSMPI_BIN="/c/Program Files/Microsoft MPI/Bin"
```
* Get the Elmer source code:
```
git clone --recurse-submodules https://github.com/ElmerCSC/elmerfem
```
* Create a build directory and a directory for the installation:
```
cd elmerfem
mkdir -p build
mkdir -p install
```
* Configure using CMake (depending on your choices change the respective flags
to `-DWITH_MPI=OFF` and/or `-DWITH_VTK=OFF`:
```
cd build
mumps_args=()
mumps_args+=("-DWITH_Mumps=ON")
mumps_args+=("-DMumps_LIBRARIES=$(pkg-config --libs mumps-dmo mumps-zmo mumps-smo mumps-cmo)")
mumps_args+=("-DMumps_INCLUDE_DIR=$(pkg-config --variable=includedir mumps-dmo mumps-zmo mumps-smo mumps-cmo)")
cmake \
  -DCMAKE_BUILD_TYPE="Release" \
  -DCMAKE_INSTALL_PREFIX="../install" \
  -DCPACK_BUNDLE_EXTRA_WINDOWS_DLLS=OFF \
  -DBLA_VENDOR="OpenBLAS" \
  -DWITH_OpenMP=ON \
  -DWITH_LUA=ON \
  -DWITH_MPI=ON \
  -DMPI_TEST_MAXPROC=$(nproc) \
  -DMPIEXEC_EXECUTABLE="$(cygpath -m "${MSMPI_BIN}")/mpiexec.exe" \
  -DWITH_Zoltan=ON \
  -DCMAKE_C_FLAGS="-Wno-error=incompatible-pointer-types" \
  -DParMetis_LIBRARIES="$(pkg-config --libs parmetis)" \
  -DParMetis_INCLUDE_DIR="$(pkg-config --cflags parmetis)" \
  "${mumps_args[@]}" \
  -DWITH_ElmerIce=ON \
  -DWITH_ELMERGUI=ON \
  -DWITH_QT6=ON \
  -DWITH_VTK=ON \
  -DWITH_OCC=ON \
  -DWITH_MATC=ON \
  -DWITH_PARAVIEW=ON \
  -DCREATE_PKGCONFIG_FILE=ON \
  ..
```
* Still inside the same directory, build the source code and install at the
  previously created prefix:
```
cmake --build .
cmake --install .
```
* Create a batch file `ElmerGUI.bat` with the following content adapting the
  values of the variables `MSY2_ENV_BIN` and `ELMERFEM_INST_BIN` as needed:
```
REM Set paths and start ElmerGUI

REM Adapt the following to the location of your MSYS2 installation
SET MSYS2_ENV_BIN=C:\msys64\ucrt64\bin

REM Adapt the following to the install prefix you created above
SET ELMERFEM_INST_BIN=C:\repo\elmerfem\install\bin

REM Do not change anything below this line

REM Adjust PATH variable
SET PATH=%ELMERFEM_INST_BIN%;%MSYS2_ENV_BIN%;%PATH%

REM Start ElmerGUI
ElmerGUI.exe
```
* To open ElmerGUI, double-click the batch file.


Install Elmer on Windows using MSYS2
====================================

Alternatively, you can install the version of ElmerFEM that is distributed by
MSYS2. This version might be (significantly) older than the latest head of the
development repository (see above). To do that:

* Install MSYS2 from https://www.msys2.org/
* Open the MSYS2 `bash` via the "MSYS2 UCRT64" shortcut.
* Download latest MSYS2 package databases and upgrade installed packages by
  running `pacman -Syu` twice.
* Install ElmerFEM from MSYS2:
```
pacman -S --needed pactoys
pacboy -S --needed elmerfem:p
```
* To use ElmerFEM, open the MSYS2 `bash` via the "MSYS2 UCRT64" shortcut.
  Then, start the GUI using the following command:
```
ElmerGUI
```


<!---
=== Package managers

[.text-center]
image::https://repology.org/badge/vertical-allrepos/elmerfem.svg["Packaging status", link=https://repology.org/project/elmerfem/versions]

==== Chocolatey

[.text-center]
image:https://img.shields.io/chocolatey/dt/elmer-mpi["Chocolatey", link=https://chocolatey.org/packages/elmer-mpi]
-->
