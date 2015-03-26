# CMake generated Testfile for 
# Source directory: /home/jpr/wrk/ff/gitelmer/fem/tests/CurvedBndryPFEM2
# Build directory: /home/jpr/wrk/ff/gitelmer/fem/tests/CurvedBndryPFEM2
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
ADD_TEST(CurvedBndryPFEM2 "/usr/bin/cmake" "-DELMERGRID_BIN=/home/jpr/wrk/ff/gitelmer/elmergrid/src/ElmerGrid" "-DELMERSOLVER_BIN=/home/jpr/wrk/ff/gitelmer/fem/src/ElmerSolver_mpi" "-DFINDNORM_BIN=/home/jpr/wrk/ff/gitelmer/fem/tests/findnorm" "-DMESH2D_BIN=/home/jpr/wrk/ff/gitelmer/meshgen2d/src/Mesh2D" "-DTEST_SOURCE=/home/jpr/wrk/ff/gitelmer/fem/tests/CurvedBndryPFEM2" "-DPROJECT_SOURCE_DIR=/home/jpr/wrk/ff/gitelmer/fem/tests" "-DBINARY_DIR=/home/jpr/wrk/ff/gitelmer" "-DCMAKE_Fortran_COMPILER=/home/jpr/wrk/gcc-4.5/bin/gfortran" "-DMPIEXEC=/home/jpr/wrk/openmpi/bin/mpiexec" "-DMPIEXEC_NUMPROC_FLAG=-np" "-DMPIEXEC_PREFLAGS=" "-DMPIEXEC_POSTFLAGS=" "-DWITH_MPI=TRUE" "-P" "/home/jpr/wrk/ff/gitelmer/fem/tests/CurvedBndryPFEM2/runtest.cmake")
SET_TESTS_PROPERTIES(CurvedBndryPFEM2 PROPERTIES  WORKING_DIRECTORY "/home/jpr/wrk/ff/gitelmer/fem/tests/CurvedBndryPFEM2")
