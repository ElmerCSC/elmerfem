
# Some implementations of MPI don't support all features with all compilers
# on all platforms (e.g., 'MPI_IN_PLACE' with MSMPI and GFortran on Windows
# without specific "glue").
# Check for those features using some tests.

if(CMAKE_CROSSCOMPILING)
  # assume it is working
  message(STATUS "Checking whether MPI_IN_PLACE is supported with ${MPI_Fortran_COMPILER} -- assuming yes")
  set(CHECK_MPI_IN_PLACE_RUN_ERROR OFF)
  set(CHECK_MPI_IN_PLACE_COMPILE ON)
else()

  message(STATUS "Checking whether MPI_IN_PLACE is supported with ${MPI_Fortran_COMPILER}")

  file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testMPI_IN_PLACE.F90
    "
    PROGRAM TEST_MPI_IN_PLACE

#if defined(ELMER_HAVE_MPI_MODULE)
      USE MPI
#endif

      IMPLICIT NONE
      INTEGER :: ierr
      REAL(8) :: test1(3), test2(3)

#if defined(ELMER_HAVE_MPIF_HEADER)
      INCLUDE \"mpif.h\"
#endif

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

  if(MPI_Fortran_HAVE_F90_MODULE)
    set(ELMER_MPI_DEFINITION "-DELMER_HAVE_MPI_MODULE")
  else()
    set(ELMER_MPI_DEFINITION "-DELMER_HAVE_MPIF_HEADER")
  endif()

  # Compile and run separately (instead of using try_run) so that the run can
  # be limited by a timeout. A broken MPI stack can hang in MPI_Init when the
  # program is started without mpiexec (e.g., Open MPI 5.0.10 with the UCX 'uct'
  # BTL from UCX 1.20), which would otherwise stall the configuration forever.
  set(ELMER_MPI_IN_PLACE_EXE
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/testMPI_IN_PLACE${CMAKE_EXECUTABLE_SUFFIX})
  file(REMOVE ${ELMER_MPI_IN_PLACE_EXE})

  try_compile(CHECK_MPI_IN_PLACE_COMPILE ${CMAKE_BINARY_DIR}
    SOURCES ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testMPI_IN_PLACE.F90
    COMPILE_DEFINITIONS ${ELMER_MPI_DEFINITION}
    LINK_LIBRARIES MPI::MPI_Fortran
    COPY_FILE ${ELMER_MPI_IN_PLACE_EXE})

  if(CHECK_MPI_IN_PLACE_COMPILE)
    execute_process(COMMAND ${ELMER_MPI_IN_PLACE_EXE}
      TIMEOUT 60
      RESULT_VARIABLE CHECK_MPI_IN_PLACE_RUN_RESULT
      OUTPUT_QUIET ERROR_QUIET)

    if(CHECK_MPI_IN_PLACE_RUN_RESULT STREQUAL "0")
      set(CHECK_MPI_IN_PLACE_RUN_ERROR OFF)
    elseif(CHECK_MPI_IN_PLACE_RUN_RESULT MATCHES "timeout")
      message(WARNING
        "The MPI_IN_PLACE test program did not finish within 60 seconds and was "
        "killed. MPI programs started without mpiexec appear to hang in MPI_Init, "
        "so ElmerSolver might hang as well when started that way. With Open MPI 5 "
        "this can be caused by the UCX 'uct' BTL; try adding the line "
        "'btl = ^uct' to ~/.openmpi/mca-params.conf (or setting the environment "
        "variable OMPI_MCA_btl=^uct). Assuming that MPI_IN_PLACE is supported.")
      set(CHECK_MPI_IN_PLACE_RUN_ERROR OFF)
    else()
      set(CHECK_MPI_IN_PLACE_RUN_ERROR ON)
    endif()
  endif()
endif()

if(CHECK_MPI_IN_PLACE_RUN_ERROR OR NOT CHECK_MPI_IN_PLACE_COMPILE)
  message(STATUS "Checking whether MPI_IN_PLACE is supported with ${MPI_Fortran_COMPILER} -- no")
  set(ELMER_BROKEN_MPI_IN_PLACE ON CACHE INTERNAL "")
else()
  message(STATUS "Checking whether MPI_IN_PLACE is supported with ${MPI_Fortran_COMPILER} -- yes")
  set(ELMER_BROKEN_MPI_IN_PLACE OFF CACHE INTERNAL "")
endif()
