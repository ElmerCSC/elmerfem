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
      REAL(KIND=dp) :: x0,h,Lf
      LOGICAL :: FirstTime=.TRUE.
      CHARACTER(len=1) :: Rien
      INTEGER :: i
     
      SAVE FirstTime, x0, h, Lf
     
      IF (FirstTime) THEN
        FirstTime = .FALSE.
        write(*,*)'--------------------------------------------------'
        write(*,*)'YOU ARE SOLVING THE SINUS PROBLEM (sinus.geo mesh)'
        write(*,*)'--------------------------------------------------'
        Open(10,file='sinus0.geo')
        DO i = 1, 5
          Read(10,*)Rien
        END DO
        Read(10,'(5x,e14.8)')x0
        Read(10,'(5x,e14.8)')h 
        Read(10,'(5x,e14.8)')Lf
        Close(10)
      END IF

         fbed = 0.5*h*cos(2.*Pi/Lf*(x-x0)) 

      END FUNCTION fbed
! --------------------------------------------------------
      FUNCTION dbed(x)
        USE types
        USE DefUtils
      IMPLICIT NONE
      REAL(KIND=dp) :: dbed,x
      REAL(KIND=dp) :: x0,h,Lf
      LOGICAL :: FirstTime=.TRUE.
      CHARACTER(len=1) :: Rien
      INTEGER :: i
     
      SAVE FirstTime, x0, h, Lf
     
      IF (FirstTime) THEN
        FirstTime = .FALSE.
        Open(10,file='sinus0.geo')
        DO i = 1, 5
          Read(10,*)Rien
        END DO
        Read(10,'(5x,e14.8)')x0
        Read(10,'(5x,e14.8)')h 
        Read(10,'(5x,e14.8)')Lf
        Close(10)
      END IF

      dbed = -h*Pi/Lf*sin(2.*Pi/Lf*(x-x0)) 

      END FUNCTION dbed
! --------------------------------------------------------
