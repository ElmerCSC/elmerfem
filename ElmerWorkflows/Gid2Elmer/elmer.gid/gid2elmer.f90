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
  INTEGER :: NodesUnit, BoundaryUnit, ElementsUnit, HeaderUnit
!----------------------------------------------------------------------------
  OPEN(NEWUNIT=NodesUnit, FILE='mesh.nodes', STATUS='unknown')
  OPEN(NEWUNIT=BoundaryUnit, FILE='mesh.boundary', STATUS='unknown')
  OPEN(NEWUNIT=ElementsUnit, FILE='mesh.elements', STATUS='unknown')
  OPEN(NEWUNIT=HeaderUnit, FILE='mesh.header', STATUS='unknown')

  DO
!----------------------------------------------------------------------------
  READ(*,'(A200)', ERR=2) lineread
  IF( lineread(1:3)=='eof' ) EXIT
  IF( lineread(1:11)=='mesh.nodes:' )    WRITE(NodesUnit,*) TRIM( lineread(12:200) )
  IF( lineread(1:14)=='mesh.boundary:' ) WRITE(BoundaryUnit,*) TRIM( lineread(15:200) )
  IF( lineread(1:14)=='mesh.elements:' ) WRITE(ElementsUnit,*) TRIM( lineread(15:200) )
  IF( lineread(1:12)=='mesh.header:' )   WRITE(HeaderUnit,*) TRIM( lineread(13:200) )
!----------------------------------------------------------------------------
  END DO

2 CONTINUE
  CLOSE(NodesUnit)
  CLOSE(BoundaryUnit)
  CLOSE(ElementsUnit)
  CLOSE(HeaderUnit)

END PROGRAM GID2ELMER
