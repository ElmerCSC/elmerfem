PROGRAM FindParents
!----------------------------------------------------------------------------
! Determines parents for boundary elements in Elmer mesh files.
!
! Written by: Mikko Lyly 18 May 2004
! Modified by: Mark Smith 7 March 2014
!----------------------------------------------------------------------------
  IMPLICIT NONE

  TYPE HashEntry_t
     INTEGER :: Element
     TYPE(HashEntry_t), POINTER :: Next
  END TYPE HashEntry_t

  TYPE HashTable_t
     TYPE(HashEntry_t), POINTER :: Head
  END TYPE HashTable_t

  TYPE(HashTable_t), ALLOCATABLE :: HashTable(:)
  TYPE(HashEntry_t), POINTER :: HashPtr, HashPtr1

  LOGICAL :: Found
  INTEGER :: i, j, a(30), Istat
  INTEGER :: Candidate, Parent(2)
  INTEGER :: NumberOfParents, Element, CheckCount, NumberOfNodes, Node
  INTEGER :: Nnodes, Nbulk, Nboundary, BoundaryCode, ElementDim, BoundaryDim
  INTEGER, ALLOCATABLE :: ElementCode(:)

!---------------------------------------------------------------------------

! Read header:
! ------------
  OPEN(10, FILE='mesh.header')
  READ(10,*) Nnodes, Nbulk, Nboundary
  CLOSE(10)

! Prepare the element hash table for nodes (inverse connectivity):
! ----------------------------------------------------------------
  OPEN(10, FILE='mesh.elements')

  ISTAT = 0
  ALLOCATE( HashTable( Nnodes ), ElementCode( Nbulk ), STAT = Istat )
  IF( Istat /= 0 ) THEN
     PRINT *,'Memory allocation error. Aborting.'
     STOP
  END IF
  
  DO i = 1,Nnodes
	NULLIFY( HashTable( i ) % Head )
  ENDDO

  DO i = 1, Nbulk
     READ(10,*) A(1:3), A(4:3+MOD(A(3),100))

!    A(1) = element number
!    A(2) = tag
!    A(3) = code
!    A(4:) = node numbers

     ElementCode(i) = A(3)
     NumberOfNodes = MOD(A(3),100)

     DO j = 1,NumberOfNodes
        Node = A(3+j)

        HashPtr => HashTable( Node ) % Head
        Found = .FALSE.

        DO WHILE( ASSOCIATED( HashPtr ) )
           IF( HashPtr % Element == A(1) ) THEN
              Found = .TRUE.
              EXIT
           END IF
           HashPtr => HashPtr % Next
        END DO

        IF( .NOT.Found ) THEN
           ALLOCATE( HashPtr )
           HashPtr % Element = i
           HashPtr % Next => HashTable( Node ) % Head
           HashTable( Node ) % Head => HashPtr
        END IF

     END DO
  END DO
  CLOSE(10)


! Make the mesh.boundary -file with parents:
! ------------------------------------------
  OPEN(10, FILE='mesh.boundary')
  OPEN(11, FILE='mesh.boundary.corrected' )
  
  DO i = 1, Nboundary
     READ(10,*) A(1:5), A(6:5+MOD(A(5),100))
     
!    A(1) = boundaryelement number
!    A(2) = tag
!    A(3) = left parent (assumed unknown)
!    A(4) = right parent (assumed unknown)
!    A(5) = code
!    A(6:) = node numbers

     BoundaryCode = A(5)
     NumberOfNodes = MOD(A(5),100)

     SELECT CASE( INT(BoundaryCode/100) )
     CASE( 1 )
        BoundaryDim = 0
     CASE( 2 )
        BoundaryDim = 1
     CASE( 3, 4 )
        BoundaryDim = 2
     CASE( 5, 6, 7, 8 )
        BoundaryDim = 3
     CASE DEFAULT
        PRINT *,'Cant detect dimension for bounbdary element',A(1)
        BoundaryDim = 0
     END SELECT
     
     NumberOfParents = 0
     HashPtr1 => HashTable( A(6) ) % Head
     Parent = 0

     DO WHILE( ASSOCIATED( HashPtr1 ) )

        Candidate = HashPtr1 % Element
        CheckCount = 0

        Do j = 1,NumberOfNodes
           Node = A(5+j)
           HashPtr => HashTable( Node ) % Head
           DO WHILE( ASSOCIATED( HashPtr ) )
              IF( HashPtr % Element == Candidate ) &
                   CheckCount = CheckCount+1
              HashPtr => HashPtr % Next
           END DO
        END DO

        IF( CheckCount == NumberOfNodes ) THEN
           NumberOfParents = NumberOfParents+1

           IF( NumberOfParents > 2 ) THEN
              PRINT *,'Confused: Found more than 2 parents'
              PRINT *,'for boundary element =',A(1)
           END IF

           SELECT CASE( INT(ElementCode(Candidate)/100) )
           CASE(1)
              ElementDim = 0
           CASE(2)
              ElementDim = 1
           CASE( 3, 4 )
              ElementDim = 2
           CASE( 5, 6, 7, 8 )
              ElementDim = 3
           CASE DEFAULT
              PRINT *,'Cant detect dimension for element',Candidate
              ElementDim = 0
           END SELECT

           IF( ElementDim /= BoundaryDim+1 ) THEN
              PRINT *,'Confused: Dimension of the boundary element',A(1)
              PRINT *,'is incompatible with possible parent',Candidate
           END IF

           IF( NumberOfParents <= 2 ) Parent( NumberOfParents ) = Candidate

        END IF

        HashPtr1 => HashPtr1 % Next
     END DO

     WRITE(11,'(100I8)') A(1:2), Parent(1:2), A(5), A(6:5+MOD(A(5),100))

  END DO

! Close, deallocate, and destroy hash tables:
! -------------------------------------------
  CLOSE(10)
  CLOSE(11)

  DEALLOCATE( ElementCode )

  DO i = 1,Nnodes
     HashPtr => HashTable(i) % Head
     DO WHILE( ASSOCIATED( HashPtr ) )
        HashPtr1 => HashPtr % Next
        DEALLOCATE( HashPtr )
        HashPtr => HashPtr1
     END DO
  END DO

! Done:
! -----

END PROGRAM FindParents
