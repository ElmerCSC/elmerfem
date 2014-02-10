      Program plotfbed
      USE Types
      USE DefUtils
      IMPLICIT NONE
      REAL(KIND=dp) :: x, dx, fbed,dbed
      INTEGER :: i,n

      Open(2,file='plotfbed.dat')
 
      x=0.0_dp
      n = 1000
      dx = 1.0_dp/n
      Do i=1,n+1
      write(2,'(3(e14.8,2x))')x,fbed(x),dbed(x)
      x= x+dx
      End Do
      Close(1)
      END PROGRAM plotfbed
