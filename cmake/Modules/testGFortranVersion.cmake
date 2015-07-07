
# Quick hack to test if the GNU Fortran version is >= 4.8

IF (NOT(CMAKE_CROSSCOMPILING))
FILE(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testGFortranVersion.F90
"
program versiontest
  character(len=2) :: major, minor
#if defined(__GFORTRAN__) && defined(__GNUC__) && defined(__GNUC_MINOR__)
  write(major, '(I2)') __GNUC__
  write(minor, '(I2)') __GNUC_MINOR__
  write(*,'(A)')  trim(adjustl(major)) // '.' // trim(adjustl(minor))
#else
  write(*,'(A)') '1.0'
#endif  

end program versiontest
")

TRY_RUN(GFORTRAN_VERSIONTEST_RUN_RESULT
  GFORTRAN_VERSIONTEST_COMPILE_RESULT
  ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testGFortranVersion.F90
  COMPILE_OUTPUT_VARIABLE GFORTRAN_VERSIONTEST_COMPILED
  RUN_OUTPUT_VARIABLE     GFORTRAN_VERSION)

IF(GFORTRAN_VERSIONTEST_COMPILE_RESULT AND
   GFORTRAN_VERSIONTEST_RUN_RESULT EQUAL 0)
  IF(GFORTRAN_VERSION VERSION_GREATER 4.7)
    MESSAGE(STATUS "Checking whether GFortran version >= 4.8 -- yes")
    SET(CMAKE_Fortran_COMPILER_GNU_VERSION_OK 1 CACHE BOOL "")
  ELSE()
    MESSAGE(STATUS "Checking whether GFortran version >= 4.8 -- no")
    SET(CMAKE_Fortran_COMPILER_GNU_VERSION_OK 0 CACHE BOOL "")
  ENDIF()
ELSE()
  MESSAGE(STATUS "Could not determine the Fortran compiler version")
  MESSAGE(STATUS "${COMPILE_RESULT}")
  MESSAGE(STATUS "${RUN_RESULT}")
ENDIF()
ELSE()
  # Handle cross-compilation without user intervention
  FILE(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testGFortranVersion.F90   
"
program versiontest
#if defined(__GFORTRAN__) && defined(__GNUC__) && defined(__GNUC_MINOR__)
#if __GNUC__<5
#if __GNUC_MINOR__<8
! Gfortran versions <4.8 are known not to fully work with Elmer
#error 
#endif
#endif
#endif
end program versiontest
")

  TRY_COMPILE(GFORTRAN_VERSIONTEST_COMPILE_RESULT ${CMAKE_BINARY_DIR}
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testGFortranVersion.F90
    OUTPUT_VARIABLE OUTPUT)

  IF(GFORTRAN_VERSIONTEST_COMPILE_RESULT)
    MESSAGE(STATUS "Checking whether GFortran version >= 4.8 -- yes")
    SET(CMAKE_Fortran_COMPILER_GNU_VERSION_OK 1 CACHE BOOL "")
  ELSE()
    MESSAGE(STATUS "Checking whether GFortran version >= 4.8 -- no")
    SET(CMAKE_Fortran_COMPILER_GNU_VERSION_OK 0 CACHE BOOL "")
  ENDIF()
ENDIF()

MARK_AS_ADVANCED(CMAKE_Fortran_COMPILER_GNU_VERSION_OK)
