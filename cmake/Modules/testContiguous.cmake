
# Quick hack to test if the Fortran compiler supports the CONTIGUOUS attribute

message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports CONTIGUOUS")
file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranCompilerContiguous.f90
  "
      PROGRAM TESTFortranCONT
      REAL, POINTER, CONTIGUOUS :: foo(:)
      REAL, ALLOCATABLE, TARGET :: a(:)
      ALLOCATE(a(100))
      foo => a
      foo(1) = 10
      END PROGRAM TESTFortranCONT
  ")
try_compile(FC_HAS_CONTIGUOUS ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranCompilerContiguous.f90
  OUTPUT_VARIABLE OUTPUT)
if(FC_HAS_CONTIGUOUS)
  message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports CONTIGUOUS -- yes")
  set(FC_HAS_CONTIGUOUS 1 CACHE INTERNAL "")
else(FC_HAS_CONTIGUOUS)
  message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports CONTIGUOUS -- no")
  set(FC_HAS_CONTIGUOUS 0 CACHE INTERNAL "")
endif(FC_HAS_CONTIGUOUS)
