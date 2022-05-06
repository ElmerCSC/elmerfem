!-------------------------------------------------------------------------------
SUBROUTINE reluctfun(Model, n, X, Y)
!-------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER :: n,i
  REAL(KIND=dp) :: X(*)
  REAL(KIND=dp), POINTER CONTIG :: Y(:,:)

  REAL(KIND=dp) :: B(3), Norm
!-------------------------------------------------------------------------------

! Here the three first values of the input X are automatically set to be
! the components of the magnetic induction B. In this example they are not
! used in the computation of the reluctivity

  B(1:3) = X(1:3)
  Norm = SQRT(SUM(B**2))
  ! PRINT *, 'B NORM', norm

  ! The offset corresponding to the automated input containing B:
  i = 3
  
  Y = 0._dp
  Y(1,1) = 0.001_dp
  Y(2,2) = 1/(SQRT(x(i+1)**2 + x(i+3)**2)*1e3 + 1e5)
  Y(3,3) = 0.001_dp
!-------------------------------------------------------------------------------
END SUBROUTINE reluctfun
!-------------------------------------------------------------------------------
