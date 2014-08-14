
# Quick hack to test if the Fortran compiler supports the CONTIGUOUS attribute

message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports CONTIGUOUS")
file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranCompilerContiguous.f90
"
  PROGRAM TESTFortranCONT	
  TYPE conttype			
     REAL, POINTER, CONTIGUOUS :: x(:) => NULL()
  END TYPE conttype
  REAL, POINTER, CONTIGUOUS :: foo(:)
  REAL, ALLOCATABLE, TARGET :: a(:)
  TYPE(conttype) :: ct
  ALLOCATE(a(100), ct % x(100))
  foo => a
  foo(1) = 10
  foo => ct % x
  foo => NULL()
  DEALLOCATE(ct % x)
  END PROGRAM TESTFortranCONT   
  ")
try_compile(FC_HAS_CONTIGUOUS ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranCompilerContiguous.f90
  OUTPUT_VARIABLE OUTPUT)
if(FC_HAS_CONTIGUOUS)
  message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports CONTIGUOUS -- yes")
  set(CMAKE_Fortran_COMPILER_SUPPORTS_CONTIGUOUS 1 CACHE BOOL "")
else(FC_HAS_CONTIGUOUS)
  message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports CONTIGUOUS -- no")
  set(CMAKE_Fortran_COMPILER_SUPPORTS_CONTIGUOUS 0 CACHE BOOL "")
endif(FC_HAS_CONTIGUOUS)
MARK_AS_ADVANCED(CMAKE_Fortran_COMPILER_SUPPORTS_CONTIGUOUS)
