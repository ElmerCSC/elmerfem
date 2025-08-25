Compilation of Elmer under Ubuntu
=================================

(note: instructions may be out of date)


 * Download the source code and create `build` directory as above
 * Install the dependencies `sudo apt install git cmake build-essential gfortran libopenmpi-dev libblas-dev liblapack-dev`
 * Without MPI:
    ** `cmake .. -DWITH_OpenMP:BOOLEAN=TRUE`
 * With MPI:
    ** `cmake .. -DWITH_OpenMP:BOOLEAN=TRUE -DWITH_MPI:BOOLEAN=TRUE`
 * With ElmerGUI:
    ** `sudo apt install libqt4-dev libqwt-dev`
    ** `cmake .. -DWITH_OpenMP:BOOLEAN=TRUE -DWITH_MPI:BOOLEAN=TRUE -DWITH_ELMERGUI:BOOLEAN=TRUE`
 * With Elmer/Ice (enabling netcdf and MUMPS):
    ** `sudo apt install libnetcdf-dev libnetcdff-dev libmumps-dev libparmetis-dev`
    ** `cmake .. -DWITH_OpenMP:BOOLEAN=TRUE -DWITH_MPI:BOOLEAN=TRUE -DWITH_ElmerIce:BOOLEAN=TRUE -DWITH_Mumps:BOOL=TRUE` 
 * `make`
 * `sudo make install`
 * The executables are in `/usr/local/bin` folder, you may add this to your PATH if you will


These are more recent instructions on discussion forum for Ubuntu 22.04:
* http://www.elmerfem.org/forum/viewtopic.php?p=28018#p28018

Following, detailed and working instructions to compile and install on Ubuntu 24.04
* Clone the repository and cd into it
<br> `git clone https://github.com/ElmerCSC/elmerfem.git` <br>
`cd elmerfem`

* Install dependencies
<br> `sudo apt install git cmake build-essential gfortran libopenmpi-dev libblas-dev liblapack-dev`

As an alternative
* Install all dependencies for complete installation:
<br> `sudo apt install git cmake build-essential gfortran libopenmpi-dev libblas-dev liblapack-dev qtscript5-dev libqt5svg5-dev libqwt-qt5-dev libnetcdf-dev libnetcdff-dev libmumps-dev libparmetis-dev`

* Create the build directory inside the repo "elmerfem" and cd into it
<br> `mkdir build`
<br> `cd build` 

* Command to create the Makefile, with the basis. No GUI, No MPI, nothing...
<br> `cmake -DWITH_ELMERGUI:BOOL=FALSE -DWITH_MPI:BOOL=FALSE -DCMAKE_INSTALL_PREFIX=<the-path-where-you-want-to-install> ../../elmerfem`

* With MPI add the following to the cmake command: ** `cmake .. -DWITH_OpenMP:BOOLEAN=TRUE -DWITH_MPI:BOOLEAN=TRUE`

* With ElmerGUI, you first need to install the following packages: ** `sudo apt install qtscript5-dev libqt5svg5-dev libqwt-qt5-dev`
then add the following options to the cmake command: `cmake .. -DWITH_OpenMP:BOOLEAN=TRUE -DWITH_MPI:BOOLEAN=TRUE -DWITH_ELMERGUI:BOOLEAN=TRUE`

* With Elmer/Ice (enabling netcdf and MUMPS), you need to install the following packages: ** `sudo apt install libnetcdf-dev libnetcdff-dev libmumps-dev libparmetis-dev`
Then add the followinf options to the cmake command: `cmake .. -DWITH_OpenMP:BOOLEAN=TRUE -DWITH_MPI:BOOLEAN=TRUE -DWITH_ElmerIce:BOOLEAN=TRUE -DWITH_Mumps:BOOL=TRUE`

* Do you want to Compile Everything?
<br> `cmake -DWITH_ELMERGUI:BOOL=TRUE -DWITH_OpenMP:BOOLEAN=TRUE -DWITH_MPI:BOOL=TRUE -DWITH_ElmerIce:BOOLEAN=TRUE -DWITH_Mumps:BOOL=TRUE -DCMAKE_INSTALL_PREFIX=<the-path-where-you-want-to-install> ../../elmerfem`

* Compile
<br> `make`

* Install
<br> `(sudo) make install`
The executables are in the folder you chose before <the-path-where-you-want-to-install> , you may add this to your PATH if you will, by modifying the file /etc/environment

* Issues with windows without borders? Modify /etc/environment by adding the following line `export QT_QPA_PLATFORM="xcb"`
