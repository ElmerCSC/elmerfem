! Sinus bed 
!
! y = bed(x) and y = d bed(x) / dx
!
! --------------------------------------------------------
      FUNCTION fbed(x)
        USE types
        USE DefUtils
      IMPLICIT NONE
      REAL(KIND=dp) :: fbed,x
      REAL(KIND=dp) :: x0,R
      LOGICAL :: FirstTime=.TRUE.
      CHARACTER(len=1) :: Rien
      INTEGER :: i
     
      SAVE FirstTime, x0, R
     
      IF (FirstTime) THEN
        FirstTime = .FALSE.
        write(*,*)'----------------------------------------------------'
        write(*,*)'YOU ARE SOLVING THE CIRCLE PROBLEM (cercle.geo mesh)'
        write(*,*)'----------------------------------------------------'
        Open(10,file='cercle0.geo')
        DO i = 1, 5
          Read(10,*)Rien
        END DO
        Read(10,'(5x,e14.8)')x0
        Read(10,'(5x,e14.8)')R 
        Close(10)
      END IF

      If (x.LE.R+x0) Then
           fbed = Sqrt(R*R - (x-x0)*(x-x0))
      Else If (x.LE.3*R+x0) Then
           fbed = - Sqrt(R*R - (x-2*R-x0)*(x-2*R-x0))
      Else 
           fbed = Sqrt(R*R - (x-4*R-x0)*(x-4*R-x0))
      End If
     

      END FUNCTION fbed
! --------------------------------------------------------
      FUNCTION dbed(x)
        USE types
        USE DefUtils
      IMPLICIT NONE
      REAL(KIND=dp) :: dbed,x,fbed
      REAL(KIND=dp) :: x0,R, DEN, EPS = 1.0e-06_dp
      LOGICAL :: FirstTime=.TRUE.
      CHARACTER(len=1) :: Rien
      INTEGER :: i
     
      SAVE FirstTime, x0, R
     
      IF (FirstTime) THEN
        FirstTime = .FALSE.
        Open(10,file='cercle0.geo')
        DO i = 1, 5
          Read(10,*)Rien
        END DO
        Read(10,'(5x,e14.8)')x0
        Read(10,'(5x,e14.8)')R 
        Close(10)
      END IF

      DEN = fbed(x)
      If (Abs(DEN) < EPS) DEN = SIGN(EPS,DEN) 
      
      If (x.LE.R+x0) Then
           dbed = -(x-x0) / DEN 
      Else If (x.LE.3*R+x0) Then
           dbed = -(x-2*R-x0) / DEN 
      Else 
           dbed = -(x-4*R-x0) / DEN 
      End If


      END FUNCTION dbed
! --------------------------------------------------------
