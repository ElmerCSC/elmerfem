! Flat bed y=0.0
!
! y = bed(x) and y = d bed(x) / dx
!
! --------------------------------------------------------
     FUNCTION fbed(x)
       USE types
       USE DefUtils
     IMPLICIT NONE
     REAL(KIND=dp) :: fbed,x

     fbed = 0.0

     END FUNCTION fbed
! --------------------------------------------------------
     FUNCTION dbed(x)
       USE types
       USE DefUtils
     IMPLICIT NONE
     REAL(KIND=dp) :: dbed,x

     dbed = 0.0

     END FUNCTION dbed
! --------------------------------------------------------
