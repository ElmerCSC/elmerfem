Compilation for macOS
=====================

(note: instructions may be out of date)


 * Download this repository either az a zip file via GitHub or using `git clone https://github.com/ElmerCSC/elmerfem.git`
 * Go to the downloaded directory `mkdir build` and `cd build`
 * Install Homebrew
 * Install GCC `brew install gcc`
 * Install CMake `brew install cmake`
 * Without MPI: 
    ** `cmake .. -D WITH_OpenMP:BOOLEAN=TRUE`
 * With MPI:
    ** Install OpenMPI `brew install open-mpi`
    ** `cmake .. -D WITH_OpenMP:BOOLEAN=TRUE -D WITH_MPI:BOOLEAN=TRUE`
 * With ElmerGUI:
    ** install qt4 with `brew install cartr/qt4/qt@4` 
    ** install qwt with `brew install cartr/qt4/qwt-qt4`
    ** `cmake .. -D WITH_OpenMP:BOOLEAN=TRUE -D WITH_MPI:BOOLEAN=TRUE -D WITH_ELMERGUI:BOOLEAN=TRUE`
 * With ElmerPost:
    ** `brew cask install xquartz`
    ** ....
 * `make`
 * `make install`

