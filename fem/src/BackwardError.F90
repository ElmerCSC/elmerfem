#include "huti_fdefs.h"

!------------------------------------------------------------------------------
!> The normwise relative backward error err = ||b-Ax||/(||A|| ||x|| + ||b||) 
!> where ||.|| is the 2-norm
!------------------------------------------------------------------------------
FUNCTION NormwiseBackwardError2( x,b,r,ipar,dpar ) RESULT(err)
!------------------------------------------------------------------------------
  USE ParallelUtils

  INTEGER :: ipar(*),n
  DOUBLE PRECISION :: x(HUTI_NDIM),b(HUTI_NDIM),r(HUTI_NDIM),dpar(*),err
  DOUBLE PRECISION :: res(HUTI_NDIM)

  n = HUTI_NDIM

  IF(ParEnv % PEs>1) THEN
    CALL SParMatrixVector(x,res,ipar)
  ELSE
    CALL CRS_MatrixVectorMultiply(GlobalMatrix,x,res)
  END IF
  res = res - b(1:n)

  err = SQRT(ParallelReduction(SUM( res(1:n)**2) )) /  &
      SQRT(ParallelReduction(SUM(GlobalMatrix % Values**2))) * &
      SQRT(ParallelReduction(SUM(x(1:n)**2))) + &
      SQRT(ParallelReduction((SUM(b(1:n)**2))) )

!------------------------------------------------------------------------------
END FUNCTION NormwiseBackwardError2
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> The normwise relative backward error err = ||b-Ax||/(||A|| ||x|| + ||b||) 
!> where ||.|| is the supremum norm and A is assumed to be scaled such that its 
!> norm is the unity (setting Linear System Row Equilibration = Logical True).
!------------------------------------------------------------------------------
FUNCTION NormwiseBackwardError( x,b,r,ipar,dpar ) RESULT(err)
!------------------------------------------------------------------------------
  USE ParallelUtils

  INTEGER :: ipar(*),n
  DOUBLE PRECISION :: x(HUTI_NDIM),b(HUTI_NDIM),r(HUTI_NDIM),dpar(*),err
  DOUBLE PRECISION :: res(HUTI_NDIM)

  n = HUTI_NDIM

  IF(ParEnv % PEs>1) THEN
    CALL SParMatrixVector(x,res,ipar)
  ELSE
    CALL CRS_MatrixVectorMultiply(GlobalMatrix,x,res)
  END IF
  res = res - b(1:n)

  err = ParallelReduction(MAXVAL(ABS(res(1:n))),2) / &
      (ParallelReduction(MAXVAL(ABS(x(1:n))),2) + &
      ParallelReduction(MAXVAL(ABS(b(1:n))),2))

!------------------------------------------------------------------------------
END FUNCTION NormwiseBackwardError
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> The normwise relative backward error err = ||b-Ax||/(||A|| ||x|| + ||b||) 
!> where ||.|| is the supremum norm. The matrix norm of A is computed within
!> this subroutine.
!------------------------------------------------------------------------------
FUNCTION NormwiseBackwardErrorGeneralized( x,b,r,ipar,dpar ) RESULT(err)
!------------------------------------------------------------------------------
  USE ParallelUtils

  INTEGER :: ipar(*),n
  DOUBLE PRECISION :: x(HUTI_NDIM),b(HUTI_NDIM),r(HUTI_NDIM),dpar(*),err
  DOUBLE PRECISION :: res(HUTI_NDIM)
  DOUBLE PRECISION :: d(HUTI_NDIM), ANorm

  n = HUTI_NDIM

  d = 1.0d0
  res = 0.0d0
  IF(ParEnv % PEs>1) THEN
    CALL SParABSMatrixVector(d,res,ipar)
  ELSE
    CALL CRS_ABSMatrixVectorMultiply(GlobalMatrix,d,res)
  END IF
  ANorm = ParallelReduction(MAXVAL(ABS(res(1:n))),2)

  res = 0.0d0
  IF(ParEnv % PEs>1) THEN
    CALL SParMatrixVector(x,res,ipar)
  ELSE
    CALL CRS_MatrixVectorMultiply(GlobalMatrix,x,res)
  END IF
  res = res - b(1:n)

  err = ParallelReduction(MAXVAL(ABS(res(1:n))),2) / &
      (ANorm * ParallelReduction(MAXVAL(ABS(x(1:n))),2) + &
      ParallelReduction(MAXVAL(ABS(b(1:n))),2))

!------------------------------------------------------------------------------
END FUNCTION NormwiseBackwardErrorGeneralized
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> The componentwise relative backward error err = max_j{|r|_j /(|A| |x| + |b|)_j}
!> with r the residual r = b-Ax.
!------------------------------------------------------------------------------
FUNCTION ComponentwiseBackwardError( x,b,r,ipar,dpar ) RESULT(err)
!------------------------------------------------------------------------------
  USE ParallelUtils

  INTEGER :: ipar(*),n
  DOUBLE PRECISION :: x(HUTI_NDIM),b(HUTI_NDIM),r(HUTI_NDIM),dpar(*),err
  DOUBLE PRECISION :: d(HUTI_NDIM),res(HUTI_NDIM)    

  n = HUTI_NDIM

  IF(ParEnv % PEs>1) THEN
    CALL SParMatrixVector(x,res,ipar)
  ELSE
    CALL CRS_MatrixVectorMultiply(GlobalMatrix,x,res)
  END IF
  res = res - b(1:n)
     
  IF(ParEnv % PEs>1) THEN
    CALL SParABSMatrixVector(ABS(x),d,ipar)
  ELSE
    CALL CRS_ABSMatrixVectorMultiply(GlobalMatrix,ABS(x),d)
  END IF
     
  d = d + ABS(b(1:n))

  err = 0.0d0
  DO i=1,n
    IF ( ABS(d(i)) < AEPS) THEN
      IF ( ABS(res(i)) > AEPS ) THEN
        err = HUGE(err)
        RETURN
      ELSE
        err = MAX(err,0.0d0)
      END IF
    ELSE
      err = MAX(err,ABS(res(i))/d(i))
    END IF
  END DO

  err = ParallelReduction(err,2)
!------------------------------------------------------------------------------
END FUNCTION ComponentwiseBackwardError
!------------------------------------------------------------------------------
