This folder contains an experimental Xdmf/HDF5 result writer module for Elmer

Only parallel version (ElmerSolver_mpi) is supported at the moment

Sample usage (for more details, refer to the solver-blocks in *.sif):

     ./compile.sh

     ElmerGrid 1 2 angle3d
     ElmerGrid 1 2 angle3d
     ElmerGrid 2 2 angle3d -metis 8

     rm results.*
     mpirun -np 8 ElmerSolver_mpi

The results are stored in the files

     results.xmf
     results.h5

The file "results.xmf" can be opened with in Paraview.

Keywords for the solver block (here N in an integer between 1-1000):

     Base File Name = String "filename"      (default: "results")
     Single Precision = Logical true/false   (default: false)
     Vector Field N = String "vectorname"    (default: none)
     Scalar Field N = String "scalarname"    (default: none)
