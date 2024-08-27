
# Some implementations of MPI don't support all features with all compilers
# on all platforms (e.g., 'MPI_IN_PLACE' with MSMPI and GFortran on Windows).
# Check for those features using some tests.

if(CMAKE_CROSSCOMPILING)
  # assume it is working
  message(STATUS "Checking whether MPI_IN_PLACE is supported with ${CMAKE_Fortran_COMPILER} -- assuming yes")
  set(CHECK_MPI_IN_PLACE_RUN_ERROR OFF)
else()

  set(save_CMAKE_Fortran_COMPILER ${CMAKE_Fortran_COMPILER})
  set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
  message(STATUS "Checking whether MPI_IN_PLACE is supported with ${CMAKE_Fortran_COMPILER}")

  file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testMPI_IN_PLACE.f90
    "
    PROGRAM TEST_MPI_IN_PLACE
      IMPLICIT NONE
      INTEGER :: ierr
      REAL(8) :: test1(3), test2(3)

      INCLUDE \"mpif.h\"

      test1(:) = 1
      test2 = test1
      PRINT *, \"test1 =\", test1

      CALL MPI_Init(ierr)
      CALL MPI_Allreduce(MPI_IN_PLACE, test1, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MPI_Finalize(ierr)
      PRINT *, \"test1 =\", test1

      IF (.NOT.ALL(test1 == test2)) THEN
        STOP 1
      END IF

    END PROGRAM TEST_MPI_IN_PLACE
    ")
  try_run(CHECK_MPI_IN_PLACE_RUN_ERROR CHECK_MPI_IN_PLACE_COMPILE ${CMAKE_BINARY_DIR}
    SOURCES ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testMPI_IN_PLACE.f90)
endif()

if(CHECK_MPI_IN_PLACE_RUN_ERROR)
  message(STATUS "Checking whether MPI_IN_PLACE is supported with ${CMAKE_Fortran_COMPILER} -- no")
  set(ELMER_BROKEN_MPI_IN_PLACE ON CACHE INTERNAL "")
else()
  message(STATUS "Checking whether MPI_IN_PLACE is supported with ${CMAKE_Fortran_COMPILER} -- yes")
  set(ELMER_BROKEN_MPI_IN_PLACE OFF CACHE INTERNAL "")
endif()

if(NOT CMAKE_CROSSCOMPILING)
  set(CMAKE_Fortran_COMPILER ${save_CMAKE_Fortran_COMPILER})
endif()
