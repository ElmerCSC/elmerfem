Name: ElmerFEM library
Version: @ELMER_FEM_VERSION@-@ELMER_FEM_REVISION@
Description: Finite element solver library Elmer
URL: elmerfem.org
prefix=@CMAKE_INSTALL_PREFIX@
libdir=${prefix}/lib/elmersolver
FC=@CMAKE_Fortran_COMPILER@
MPIF90=@MPI_Fortran_COMPILER@
includedir=${prefix}/share/elmersolver/include
Libs: -L${libdir} -lelmersolver
Cflags: -Wl,-rpath=${libdir} @CMAKE_Fortran_FLAGS@ @ELMER_F90FLAGS@ -I${includedir}

