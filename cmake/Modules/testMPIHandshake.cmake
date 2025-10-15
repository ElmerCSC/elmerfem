MESSAGE(STATUS "Checking whether XIOS supports MPI Handshake")

# Define the test source code (Fortran code) for the MPI Handshake check
FILE(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/test_xios_mpi_handshake.f90"
"
PROGRAM test_mpi_handshake
  ! Try importing mpi_handshake modules
  USE xios, ONLY: xios_mpi_handshake, xios_MAX_GROUPNAME_LEN
  IMPLICIT NONE
END PROGRAM test_mpi_handshake
")

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I${XIOS_INCLUDE_DIR}")
TRY_COMPILE(TEST_XIOS_HAS_MPI_HANDSHAKE
    ${CMAKE_BINARY_DIR}
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/test_xios_mpi_handshake.f90
    OUTPUT_VARIABLE OUTPUT
)

# Check if the compilation was successful
IF(TEST_XIOS_HAS_MPI_HANDSHAKE)
  message(STATUS "Checking whether XIOS supports MPI Handshake -- yes")
  set(XIOS_HAS_MPI_HANDSHAKE 1 CACHE BOOL "")
ELSE(TEST_XIOS_HAS_MPI_HANDSHAKE)
  message(STATUS "Checking whether XIOS supports MPI Handshake -- no")
  set(XIOS_HAS_MPI_HANDSHAKE 0 CACHE BOOL "")
ENDIF(TEST_XIOS_HAS_MPI_HANDSHAKE)
MARK_AS_ADVANCED(XIOS_HAS_MPI_HANDSHAKE)
