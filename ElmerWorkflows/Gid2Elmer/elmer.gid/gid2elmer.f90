PROGRAM GID2ELMER
!----------------------------------------------------------------------------
! READS GID MESH DATA FROM STDIO AND WRITES IT OUT IN ELMER MESH FILES
!
! Written by : Mikko Lyly 17 May 2005
!
! Modified by:
!----------------------------------------------------------------------------
  IMPLICIT NONE
  CHARACTER(LEN=200) :: lineread
!----------------------------------------------------------------------------
  OPEN(UNIT=10, FILE='mesh.nodes', STATUS='unknown')
  OPEN(UNIT=11, FILE='mesh.boundary', STATUS='unknown')
  OPEN(UNIT=12, FILE='mesh.elements', STATUS='unknown')
  OPEN(UNIT=13, FILE='mesh.header', STATUS='unknown')

1 CONTINUE
!----------------------------------------------------------------------------
  READ(*,'(A200)', ERR=2) lineread
  IF( lineread(1:3)=='eof' ) GOTO 2
  IF( lineread(1:11)=='mesh.nodes:' )    WRITE(10,*) TRIM( lineread(12:200) )
  IF( lineread(1:14)=='mesh.boundary:' ) WRITE(11,*) TRIM( lineread(15:200) )
  IF( lineread(1:14)=='mesh.elements:' ) WRITE(12,*) TRIM( lineread(15:200) )
  IF( lineread(1:12)=='mesh.header:' )   WRITE(13,*) TRIM( lineread(13:200) )
!----------------------------------------------------------------------------
  GOTO 1

2 CONTINUE
  CLOSE(10)
  CLOSE(11)
  CLOSE(12)
  CLOSE(13)

END PROGRAM GID2ELMER
