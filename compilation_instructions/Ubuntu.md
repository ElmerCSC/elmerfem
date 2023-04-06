==== Ubuntu

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
