!-------------------------------------------------------------------------------
FUNCTION a(Model, n, X) RESULT(Y)
!-------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: X(2)
  REAL(KIND=dp) :: Y
  real(kind=dp) :: a_mult
  TYPE(Valuelist_t), POINTER :: Params
  logical :: found
!-------------------------------------------------------------------------------
  Params => ListGetSolverParams()
  a_mult = ListGetCReal(Params, 'a_mult', found)
  y = 0._dp  
  if (.not. found) return
  Y = X(1)*a_mult
!-------------------------------------------------------------------------------
END FUNCTION a
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
FUNCTION m(Model, n, X) RESULT(Y)
!-------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: X(3)
  REAL(KIND=dp) :: Y
  real(kind=dp) :: a_mult
  TYPE(Valuelist_t), POINTER :: Params
  logical :: found
!-------------------------------------------------------------------------------
  Params => ListGetSolverParams()
  a_mult = ListGetCReal(Params, 'a_mult', found)
  y = 0._dp  
  if (.not. found) return
  Y = a_mult*1.26d6
!-------------------------------------------------------------------------------
END FUNCTION m
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
FUNCTION jplus(Model, n, X) RESULT(Y) ! {{{
!-------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: X(3)
  REAL(KIND=dp) :: Y
  real(kind=dp) :: a_mult
  TYPE(Valuelist_t), POINTER :: Params
  logical :: found
!-------------------------------------------------------------------------------

  Params => ListGetSolverParams()
  a_mult = ListGetCReal(Params, 'a_mult', found)
  y = 0._dp
  if (.not. found) return
  Y = a_mult
!-------------------------------------------------------------------------------
END FUNCTION ! }}}
!-------------------------------------------------------------------------------
FUNCTION jminus(Model, n, X) RESULT(Y) ! {{{
!-------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: X(3)
  REAL(KIND=dp) :: Y
  real(kind=dp) :: a_mult
  TYPE(Valuelist_t), POINTER :: Params
  logical :: found
!-------------------------------------------------------------------------------

  Params => ListGetSolverParams()
  a_mult = ListGetCReal(Params, 'a_mult', found)
  y = 0._dp
  if (.not. found) return
  Y = -a_mult
!-------------------------------------------------------------------------------
END FUNCTION ! }}}
!-------------------------------------------------------------------------------
