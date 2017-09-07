PROGRAM SurfaceMap
!----------------------------------------------------------------
! This program is used to map a given plane mesh to a surface 
! in the threee-dimensional point space. It can also be used to 
! create mesh.director file which defines the director at the
! nodes in a similar way as nodes are defined. 
! To be decided: should mesh.director contain also nodes data
!---------------------------------------------------------------
  IMPLICIT NONE 

  LOGICAL :: CreateDirector

  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12)

  INTEGER, ALLOCATABLE :: Node(:)
  INTEGER :: n, i, j, MappingMode
  
  REAL(KIND=dp), ALLOCATABLE :: x(:), y(:), z(:)
  REAL(KIND=dp) :: R, xi, yi, zi, t, d1, d2, d3
!-------------------------------------------------
  CreateDirector = .FALSE.
  MappingMode = 5 ! 1 = polar mapping, 
                  ! 2 = shear + polar mapping,
                  ! 3 = polar mapping & create mesh.director file
                  ! 4 = shear + polar mapping & create mesh.director file
                  ! 5 = plate body & create mesh.director file
  IF (MappingMode == 3 .OR. MappingMode == 4 .OR. MappingMode == 5) CreateDirector = .TRUE.
  R = 1.0d0

  OPEN(11,file="mesh.header")
  READ(11,*) N
  CLOSE(11)
  ALLOCATE (Node(N), x(N), y(N), z(N))

  OPEN(12,file="mesh.nodes")
  DO i=1,N
    READ(12,*) Node(i),j,x(i),y(i),z(i)
  END DO
  REWIND(12)

  IF (CreateDirector) OPEN(13,file="mesh.director")

  SELECT CASE(MappingMode)
  CASE(1)
     DO i=1,N
        xi = R * SIN(x(i)/R)
        yi = y(i)
        zi = R*(1.0d0-COS(x(i)/R))
        WRITE(12,1200) Node(i),j,xi,yi,zi
     END DO
  CASE(2)
     DO i=1,N
        xi = x(i)
        yi = y(i)
        zi = z(i)
        x(i) = xi + 0.25d0 * yi
        y(i) = yi
        z(i) = zi
        
        WRITE(12,1200) Node(i),j,R*SIN(x(i)/R), y(i), R*(1.0d0-COS(x(i)/R))
     END DO
  CASE(3)
     DO i=1,N
        t = x(i)/R
        xi = R * SIN(t)
        yi = y(i)
        zi = R*(1.0d0-COS(t))
        WRITE(12,1200) Node(i),j,xi,yi,zi
        d1 = R * SIN(t)
        d2 = 0.0d0
        d3 = -R * cos(t)
        WRITE(13,1300) Node(i),d1,d2,d3 
     END DO
  CASE(4)
     DO i=1,N
        xi = x(i)
        yi = y(i)
        zi = z(i)
        x(i) = xi + 0.25d0 * yi
        y(i) = yi
        z(i) = zi

        t = x(i)/R
        WRITE(12,1200) Node(i),j,R*SIN(t), y(i), R*(1.0d0-COS(t))
        d1 = R * SIN(t)
        d2 = 0.0d0
        d3 = -R * cos(t)
        WRITE(13,1300) Node(i),d1,d2,d3 
     END DO
   CASE(5)
     DO i=1,N
       WRITE(12,1200) Node(i),j,x(i),y(i),z(i)
       WRITE(13,1300) Node(i),0.0d0,0.0d0,1.0d0
     END DO
  END SELECT
  CLOSE(12)
  IF (CreateDirector) CLOSE(13)

  DEALLOCATE (Node, x, y, z)
  
1200 FORMAT(i6,2x,i5,3(2x,e22.15)) 
1300 FORMAT(i6,3(2x,e22.15)) 

END PROGRAM SurfaceMap
