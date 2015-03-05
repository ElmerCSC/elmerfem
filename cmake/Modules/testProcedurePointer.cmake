# Quick hack to test if the Fortran compiler supports procedure pointer

message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports PROCEDURE POINTER")
file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranProcedurePointer.f90
"
FUNCTION addrof(fn) RESULT(faddr)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        SUBROUTINE dummysubr(x,y) bind(C)
          IMPORT
          IMPLICIT NONE
          REAL(KIND=C_FLOAT) :: x, y
        END SUBROUTINE dummysubr
    END INTERFACE
    PROCEDURE(dummysubr) :: fn  
    INTEGER(KIND=8) :: faddr
    TYPE(C_FUNPTR) :: cptr
    cptr = C_FUNLOC(fn)
    faddr = TRANSFER(cptr, faddr)
END FUNCTION addrof
SUBROUTINE mysin(x,y) BIND(C)
  USE ISO_C_BINDING
  IMPLICIT NONE
  REAL(KIND=C_FLOAT) :: x, y
  INTRINSIC SIN
  y=SIN(x)
END SUBROUTINE mysin
PROGRAM TESTFortranProcPtr
  USE ISO_C_BINDING
  IMPLICIT NONE
  ABSTRACT INTERFACE
    SUBROUTINE trig(x,y) BIND(C)
      IMPORT
      IMPLICIT NONE
      REAL(KIND=C_FLOAT) :: x, y
    END SUBROUTINE trig
  END INTERFACE
  INTERFACE
    SUBROUTINE mysin(x,y) BIND(C)
      IMPORT
      IMPLICIT NONE
      REAL(KIND=C_FLOAT) :: x, y
    END SUBROUTINE mysin
  END INTERFACE
  TYPE(C_FUNPTR) :: cfptr
  PROCEDURE(trig), POINTER :: trigf
  REAL(KIND=C_FLOAT) :: x, y, yc
  INTEGER(KIND=8) :: addrof, addr
  EXTERNAL :: addrof
  ! Mimic Elmer behaviour by transferring a function pointer to a 
  ! C pointer with TRANSFER, transfer back to Fortran pointer and 
  ! then call the function
  addr = addrof(mysin)
  cfptr = TRANSFER(addr,cfptr) 
  CALL C_F_PROCPOINTER(cfptr, trigf)
  x = real(3.14/2,C_FLOAT)
  CALL trigf(x,y)
  CALL mysin(x,yc)
  write (*,*) y, yc
END PROGRAM TESTFortranProcPtr
  ")
try_compile(FC_HAS_PROCEDUREPOINTER ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranProcedurePointer.f90
  OUTPUT_VARIABLE OUTPUT)
if(FC_HAS_PROCEDUREPOINTER)
  message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports PROCEDURE POINTER -- yes")
  set(CMAKE_Fortran_COMPILER_SUPPORTS_PROCEDUREPOINTER 1 CACHE BOOL "")
else()
  message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports PROCEDURE POINTER -- no")
  set(CMAKE_Fortran_COMPILER_SUPPORTS_PROCEDUREPOINTER 0 CACHE BOOL "")
endif()
MARK_AS_ADVANCED(CMAKE_Fortran_COMPILER_SUPPORTS_PROCEDUREPOINTER)
