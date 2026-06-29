!/******************************************************************************
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! *
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! *
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write
! *  to the Free Software Foundation, Inc., 51 Franklin Street,
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{

!------------------------------------------------------------------------------
!>  Mesh boundary condition tagging: intersection BCs, rule-based tagging, body tagging.
!>  Extracted from MeshUtils.
!------------------------------------------------------------------------------

MODULE MeshTagging

    USE MeshBasics
    IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------

  SUBROUTINE CreateIntersectionBCs(Model, Mesh)

    TYPE(Model_t) :: Model
    TYPE(Mesh_t) :: Mesh
    TYPE(Element_t), POINTER :: Element, Element2, Enew, Face, Face2, Parent
    INTEGER, POINTER :: NodeIndexes(:), NodeIndexes2(:), EdgeIndexes(:), EdgeIndexes2(:), ParentBCs(:)
    INTEGER :: i,i2,j,j2,k,k2,e,e2,l,n,n2,m,nbc,nbulk,nold,t,t2,istat,newbcs,newcnt,bc_id
    TYPE(Element_t), POINTER :: NewElements(:)
    TYPE(ValueList_t), POINTER :: BC
    INTEGER, ALLOCATABLE :: BoundaryId(:), IntersectionBCs(:,:)
    LOGICAL, ALLOCATABLE :: EdgeDone(:), NodeDone(:)
    LOGICAL :: Found, Hit, EdgesPresent, NeedEdges
    
    ! Count how many of the BCs are intersection BCs that we need to determine
    j = 0
    DO bc_id=1,Model % NumberOfBCs
      BC => Model % BCs(bc_id) % Values
      IF( ListCheckPresent( BC,'Intersection BC' ) .OR. &
          ListCheckPresent( BC,'Intersection Body') ) j = j+1 
    END DO
    NewBCs = j
    IF(NewBCs==0) RETURN

    CALL Info('CreateIntersectionBCs',&
        'Number of intersection BCs to determine: '//I2S(NewBCs),Level=5)

    ! Create a fast look-up table that define the new BC indexes and the parent BCs
    ALLOCATE(IntersectionBCs(j,5))
    IntersectionBCs = 0
    j = 0
    DO bc_id=1,Model % NumberOfBCs
      BC => Model % BCs(bc_id) % Values
      ParentBCs => ListGetIntegerArray( BC,'Intersection BC',Found )
      k = 0 
      IF(.NOT. Found ) THEN
        ! If the intersection is between two bodies mark it separately!
        ParentBCs => ListGetIntegerArray( BC,'Intersection Body',Found )
        k = 1
      END IF
      IF(.NOT. Found) CYCLE
      j = j + 1
      IF(SIZE(ParentBCs) /= 2 ) CALL Fatal('CreateIntersectionBCs','Only available for two parents!')
      IntersectionBCs(j,1) = Model % BCs(bc_id) % Tag
      IntersectionBCs(j,2:3) = ParentBCs(1:2)
      IntersectionBCs(j,4) = k
    END DO

    nbulk = Mesh % NumberOfBulkElements
    nbc = Mesh % NumberOfBoundaryElements
    nold = nbulk + nbc

    ! If we need to find intersection between boundaries create a helper structure. 
    IF( ANY(IntersectionBCs(:,4) == 0) ) THEN
      ALLOCATE( BoundaryId( nbc ) )
      BoundaryId = 0
     
      DO t=1,nbc
        Element => Mesh % Elements(nbulk+t)

        ! Only treat 2D boundary elements
        IF( Mesh % MeshDim == 3 ) THEN
          IF( Element % TYPE % ElementCode < 300 ) CYCLE
        ELSE
          IF( Element % TYPE % ElementCode < 200 ) CYCLE
        END IF

        DO bc_id=1,Model % NumberOfBCs
          IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc_id) % Tag ) EXIT
        END DO
        IF ( bc_id > Model % NumberOfBCs ) CYCLE

        IF( ANY(IntersectionBCs(:,2)==bc_id .OR. IntersectionBCs(:,3)==bc_id)) THEN
          BoundaryId(t) = bc_id
        END IF
      END DO

      n = COUNT( BoundaryId > 0 )
      CALL Info('CreateIntersectionBCs','Number of candidate intersection parents: '//I2S(n))
    END IF
      
    ! Go the new boundary elements over two times.
    ! On the 1st loop just count the number of new elements.
    ! On the 2nd lopp add the new elements in the element list.
    !-------------------------------------------------------------    
    EdgesPresent = ASSOCIATED( Mesh % Edges )
    NeedEdges = (Mesh % MeshDim == 3 .OR. ANY(IntersectionBCs(:,4)==1) )

    IF(NeedEdges .AND. .NOT. EdgesPresent ) THEN
      CALL Info('CreateInterectionsBCs','Need edges for speedy search!',Level=7)
      CALL FindMeshEdges( Mesh ) 
    END IF

    IF( Mesh % MeshDim == 3 ) THEN
      ALLOCATE( EdgeDone( Mesh % NumberOfEdges ) )       
      CALL CreateIntersection3D(.TRUE.,NewCnt)
      IF(NewCnt==0) THEN
        CALL Info('CreateIntersectionBCs','Could not find any additional interface elements!')        
        GOTO 1
      END IF
      CALL CreateIntersection3D(.FALSE.,NewCnt)
    ELSE
      IF(ANY(IntersectionBCs(:,4) == 0 )) THEN      
        ALLOCATE( NodeDone( Mesh % NumberOfNodes ) )
      END IF
      CALL CreateIntersection2D(.TRUE.,NewCnt)
      IF(NewCnt==0) THEN
        CALL Info('CreateIntersectionBCs','Could not find any additional interface elements!')
        GOTO 1
      END IF
      CALL CreateIntersection2D(.FALSE.,NewCnt)
    END IF
               
    IF( InfoActive(10) ) THEN
      DO i=1,newbcs
        CALL Info('CreateIntersectionBCs','New boundary '//I2S(IntersectionBCs(i,1))//&
            ' with '//I2S(IntersectionBCs(i,5))//' elements')
      END DO
    END IF

1   IF(NeedEdges .AND. .NOT. EdgesPresent ) THEN
      CALL ReleaseMeshEdgeTables( Mesh )
      CALL ReleaseMeshFaceTables( Mesh )
    END IF
      
    CALL Info('CreateIntersectionBCs','All done!',Level=12)

  CONTAINS


    ! Find intersection between 2D boundaries i.e. the result will
    ! be a new 1D boundary, or between 3D bodies and the result will
    ! be a new 2D boundary. 
    !-------------------------------------------------------------
    SUBROUTINE CreateIntersection3D(AllocateOnly,NewCnt)
      LOGICAL :: AllocateOnly
      INTEGER :: NewCnt
      LOGICAL :: BulkMode
            
      EdgeDone = .FALSE.
      NewCnt = 0

      DO l=1,newbcs
        BulkMode = ( IntersectionBCs(l,4) == 1)
        
        IF( BulkMode ) THEN
          DO t=1,Mesh % NumberOfFaces
            Element => Mesh % Faces(t)
            IF(.NOT. ASSOCIATED(Element % BoundaryInfo)) CYCLE

            j = 0
            IF(ASSOCIATED(Element % BoundaryInfo % Left)) THEN
              j = Element % BoundaryInfo % Left % BodyId
            END IF
            j2 = 0
            IF(ASSOCIATED(Element % BoundaryInfo % Right)) THEN
              j2 = Element % BoundaryInfo % Right % BodyId
            END IF

            IF(j == j2) CYCLE
            IF( ALL( IntersectionBCs(l,2:3) /= j)) CYCLE
            IF( ALL( IntersectionBCs(l,2:3) /= j2)) CYCLE

            NodeIndexes => Element % NodeIndexes
            n = Element % TYPE % NumberOfNodes

            NewCnt = NewCnt + 1

            IF(.NOT. AllocateOnly ) THEN
              Enew => Mesh % Elements(nold+NewCnt)        
              ALLOCATE(Enew % BoundaryInfo)         
              Enew % PartIndex = Element % PartIndex
              Enew % ElementIndex = nold + NewCnt
              Enew % TYPE => GetElementType(Element % Type % ElementCode)

              CALL AllocateVector( ENew % NodeIndexes, n)
              Enew % NodeIndexes = NodeIndexes
              Enew % NDOFs = 1
              Enew % BoundaryInfo % Constraint = IntersectionBCs(l,1)
              IF(j==0) THEN
                Enew % BoundaryInfo % Left => NULL()
              ELSE
                Enew % BoundaryInfo % Left => Element % BoundaryInfo % Left 
              END IF
              IF(j2==0) THEN
                Enew % BoundaryInfo % Right => NULL()
              ELSE
                Enew % BoundaryInfo % Right => Element % BoundaryInfo % Right 
              END IF

              Enew % EdgeIndexes => NULL()
              Enew % FaceIndexes => NULL()
              Enew % PDefs => NULL()
              Enew % BubbleIndexes => NULL()

              IntersectionBCs(l,5) = IntersectionBCs(l,5) + 1
            END IF
          END DO
        ELSE
          DO t=1,nbc
            j = BoundaryId(t) 
            IF(j==0) CYCLE
            ! Do we have a suitable pair of indexes for the parents
            IF( ALL( IntersectionBCs(l,2:3) /= j)) CYCLE

            Element => Mesh % Elements(nbulk+t)
            NodeIndexes => Element % NodeIndexes
            n = Element % TYPE % NumberOfNodes

            DO t2=t+1,nbc
              j2 = BoundaryId(t2) 
              IF(j2==0) CYCLE
              IF(j==j2) CYCLE
              IF( ALL( IntersectionBCs(l,2:3) /= j2)) CYCLE

              Element2 => Mesh % Elements(nbulk+t2)
              NodeIndexes2 => Element2 % NodeIndexes
              n2 = Element2 % TYPE % NumberOfNodes

              ! Do we have any common nodes. Some are required...
              k = 0
              DO i=1,n
                IF( ANY(NodeIndexes(i) == NodeIndexes2(1:n2) ) ) k = k+1
              END DO
              IF(k<2) CYCLE

              EdgeIndexes => Element % EdgeIndexes
              IF(ASSOCIATED(EdgeIndexes)) THEN
                Face => Element
              ELSE
                Face => NULL()
                IF( ASSOCIATED( Element % BoundaryInfo % Left ) ) THEN
                  Face => Find_Face(Mesh, Element % BoundaryInfo % Left, Element )                
                END IF
                IF(.NOT. ASSOCIATED(Face) ) THEN
                  IF( ASSOCIATED( Element % BoundaryInfo % Right ) ) THEN
                    Face => Find_Face(Mesh, Element % BoundaryInfo % Right, Element )
                  END IF
                END IF
                IF(ASSOCIATED( Face ) ) THEN
                  EdgeIndexes => Face % EdgeIndexes
                ELSE
                  CALL Fatal('CreateIntersectionBCs','EdgeIndexes not associated!')
                END IF
              END IF

              ! This is a probably candidate as we have two 2D elements of proper type
              ! sharing at least two nodes. Just have to find for which edges the intersection
              ! applies. It could sometimes be a false positive also. 
              DO i=1,Face % TYPE % NumberOfEdges
                e = EdgeIndexes(i)          
                IF( EdgeDone(e) ) CYCLE

                EdgeIndexes2 => Element2 % EdgeIndexes
                IF(ASSOCIATED(EdgeIndexes2)) THEN
                  Face2 => Element2
                ELSE
                  Face2 => NULL()
                  IF( ASSOCIATED( Element2 % BoundaryInfo % Left ) ) THEN
                    Face2 => Find_Face(Mesh, Element2 % BoundaryInfo % Left, Element2 )
                  END IF
                  IF(.NOT. ASSOCIATED(Face2) ) THEN
                    IF( ASSOCIATED( Element2 % BoundaryInfo % Right ) ) THEN
                      Face2 => Find_Face(Mesh, Element2 % BoundaryInfo % Right, Element2 )
                    END IF
                  END IF
                  IF(ASSOCIATED( Face2 ) ) THEN
                    EdgeIndexes2 => Face2 % EdgeIndexes
                  ELSE
                    CALL Fatal('CreateIntersectionBCs','EdgeIndexes2 not associated!')
                  END IF
                END IF

                DO i2=1,Face2 % TYPE % NumberOfEdges
                  e2 = EdgeIndexes2(i2)

                  ! Ok, we have a hit. Same edge appearing in the proper parent
                  ! boundary elements. Create the actual boundary element only if
                  ! we have already allocated for it. 
                  IF(e==e2) THEN
                    EdgeDone(e) = .TRUE.
                    NewCnt = NewCnt + 1

                    IF(.NOT. AllocateOnly ) THEN
                      Enew => Mesh % Elements(nold+NewCnt)        
                      ALLOCATE(Enew % BoundaryInfo)         
                      Enew % PartIndex = Element % PartIndex
                      Enew % ElementIndex = nold + NewCnt

                      Enew % TYPE => Mesh % Edges(e) % TYPE

                      m = Enew % TYPE % NumberOfNodes
                      CALL AllocateVector( ENew % NodeIndexes, m)
                      Enew % NodeIndexes = Mesh % Edges(e) % NodeIndexes
                      Enew % NDOFs = m
                      Enew % BoundaryInfo % Constraint = IntersectionBCs(l,1)
                      Enew % BoundaryInfo % Left => Element
                      Enew % BoundaryInfo % Right => Element2

                      Enew % EdgeIndexes => NULL()
                      Enew % FaceIndexes => NULL()
                      Enew % PDefs => NULL()
                      Enew % BubbleIndexes => NULL()

                      ! Just a simple counter for the new BCs of this type
                      IntersectionBCs(l,4) = IntersectionBCs(l,4) + 1
                    END IF

                    EXIT
                  END IF
                END DO

              END DO
            END DO
          END DO
        END IF
      END DO
        
      ! There is nothing to do since no new elements will be created.
      IF( NewCnt == 0 ) RETURN
              
      IF(AllocateOnly) THEN
        ALLOCATE( NewElements(nold + NewCnt ) )
        CALL Info('CreateIntersectionBCs','Allocated for '//I2S(NewCnt)//' new 1D boundary elements!',Level=6)
        
        NewElements(1:nold) = Mesh % Elements(1:nold)

        DO i=nbulk+1,nold
          Element => Mesh % Elements(i)        
          IF(.NOT. ASSOCIATED(Element % BoundaryInfo)) CYCLE

          Parent => Element % BoundaryInfo % Left
          IF(ASSOCIATED(Parent)) THEN
            NewElements(i) % BoundaryInfo % Left => NewElements(Parent % ElementIndex)
          END IF

          Parent => Element % BoundaryInfo % Right
          IF(ASSOCIATED(Parent)) THEN
            NewElements(i) % BoundaryInfo % Right => NewElements(Parent % ElementIndex)
          END IF
        END DO

        DO t=1,Mesh % NumberOfFaces
          Element => Mesh % Faces(t)
          IF(ASSOCIATED(Element % BoundaryInfo % Left)) THEN
            Element % BoundaryInfo % Left => &
                NewElements(Element % BoundaryInfo % Left % ElementIndex)
          END IF
          IF(ASSOCIATED(Element % BoundaryInfo % Right)) THEN
            Element % BoundaryInfo % Right => &
                NewElements(Element % BoundaryInfo % Right % ElementIndex)
          END IF
        END DO
                
        DEALLOCATE(Mesh % Elements)
        Mesh % Elements => NewElements
        Mesh % NumberOfBoundaryElements = nbc + NewCnt
      END IF

    END SUBROUTINE CreateIntersection3D

    
    ! Find intersection between 1D boundaries i.e. the result will
    ! be a new 0D boundary (=node). 
    !-------------------------------------------------------------
    SUBROUTINE CreateIntersection2D(AllocateOnly,NewCnt)
      LOGICAL :: AllocateOnly
      INTEGER :: NewCnt
      LOGICAL :: BulkMode
      
      NewCnt = 0

      DO l=1,newbcs
        BulkMode = ( IntersectionBCs(l,4) == 1)
        
        IF( BulkMode ) THEN
          DO t=1,Mesh % NumberOfEdges
            Element => Mesh % Edges(t)
            IF(.NOT. ASSOCIATED(Element % BoundaryInfo)) CYCLE

            j = 0
            IF(ASSOCIATED(Element % BoundaryInfo % Left)) THEN
              j = Element % BoundaryInfo % Left % BodyId
            END IF

            j2 = 0
            IF(ASSOCIATED(Element % BoundaryInfo % Right)) THEN
              j2 = Element % BoundaryInfo % Right % BodyId
            END IF

            IF(j == j2) CYCLE
            IF( ALL( IntersectionBCs(l,2:3) /= j)) CYCLE
            IF( ALL( IntersectionBCs(l,2:3) /= j2)) CYCLE

            NodeIndexes => Element % NodeIndexes
            n = Element % TYPE % NumberOfNodes

            NewCnt = NewCnt + 1

            IF(.NOT. AllocateOnly ) THEN
              Enew => Mesh % Elements(nold+NewCnt)        
              ALLOCATE(Enew % BoundaryInfo)         
              !Enew % PartIndex = Element % PartIndex
              Enew % ElementIndex = nold + NewCnt
              Enew % TYPE => GetElementType(Element % Type % ElementCode)

              CALL AllocateVector( ENew % NodeIndexes, n)
              Enew % NodeIndexes = NodeIndexes
              Enew % NDOFs = 1
              Enew % BoundaryInfo % Constraint = IntersectionBCs(l,1)
              IF(j==0) THEN
                Enew % BoundaryInfo % Left => NULL()
              ELSE
                Enew % BoundaryInfo % Left => Element % BoundaryInfo % Left 
              END IF
              IF(j2==0) THEN
                Enew % BoundaryInfo % Right => NULL()
              ELSE
                Enew % BoundaryInfo % Right => Element % BoundaryInfo % Right 
              END IF

              Enew % EdgeIndexes => NULL()
              Enew % FaceIndexes => NULL()
              Enew % PDefs => NULL()
              Enew % BubbleIndexes => NULL()

              IntersectionBCs(l,5) = IntersectionBCs(l,5) + 1
            END IF
          END DO
        ELSE
          DO t=1,nbc
            NodeDone = .FALSE.
            j = BoundaryId(t) 
            IF(j==0) CYCLE
            IF( ALL( IntersectionBCs(l,2:3) /= j)) CYCLE

            Element => Mesh % Elements(nbulk+t)
            NodeIndexes => Element % NodeIndexes
            n = Element % TYPE % NumberOfNodes

            DO t2=t+1,nbc
              j2 = BoundaryId(t2)
              IF(j2==0) CYCLE

              IF(j==j2) CYCLE
              IF( ALL( IntersectionBCs(l,2:3) /= j2)) CYCLE

              Element2 => Mesh % Elements(nbulk+t2)
              NodeIndexes2 => Element2 % NodeIndexes
              n2 = Element2 % TYPE % NumberOfNodes

              ! Ok, so BC elements t and t2 have suitable indeces. 
              ! Do we have any common nodes. Some are required...
              k = 0
              e = 0
              DO i=1,n
                IF( ANY(NodeIndexes(i) == NodeIndexes2(1:n2) ) ) THEN
                  e = NodeIndexes(i)
                  k = k+1
                END IF
              END DO
              IF(k/=1) CYCLE

              IF(NodeDone(e)) CYCLE
              NodeDone(e) = .TRUE.
              NewCnt = NewCnt + 1
              
              PRINT *,'BC node:',l,t,t2,e,n
              
              IF(.NOT. AllocateOnly ) THEN
                Enew => Mesh % Elements(nold+NewCnt)        
                ALLOCATE(Enew % BoundaryInfo)         
                Enew % PartIndex = Element % PartIndex
                Enew % ElementIndex = nold + NewCnt
                !Enew % TYPE => Element % Type 
                Enew % TYPE => GetElementType(101)

                CALL AllocateVector( ENew % NodeIndexes, 1)
                Enew % NodeIndexes = e
                Enew % NDOFs = 1
                Enew % BoundaryInfo % Constraint = IntersectionBCs(l,1)
                Enew % BoundaryInfo % Left => Element
                Enew % BoundaryInfo % Right => Element2

                Enew % EdgeIndexes => NULL()
                Enew % FaceIndexes => NULL()
                Enew % PDefs => NULL()
                Enew % BubbleIndexes => NULL()
                
                IntersectionBCs(l,5) = IntersectionBCs(l,5) + 1
              END IF
            END DO
          END DO
        END IF
      END DO

      IF( NewCnt == 0 ) RETURN
              
      IF(AllocateOnly) THEN
        ALLOCATE( NewElements(nold + NewCnt ) )
        CALL Info('CreateIntersectionBCs','Allocated for '//I2S(NewCnt)//' new boundary elements in 2D!',Level=6)
        
        NewElements(1:nold) = Mesh % Elements(1:nold)
        
        DO i=nbulk+1,nold
          Element => Mesh % Elements(i)        
          IF(.NOT. ASSOCIATED(Element % BoundaryInfo)) CYCLE

          Parent => Element % BoundaryInfo % Left
          IF(ASSOCIATED(Parent)) THEN
            NewElements(i) % BoundaryInfo % Left => NewElements(Parent % ElementIndex)
          END IF

          Parent => Element % BoundaryInfo % Right
          IF(ASSOCIATED(Parent)) THEN
            NewElements(i) % BoundaryInfo % Right => NewElements(Parent % ElementIndex)
          END IF
        END DO

        DO t=1,Mesh % NumberOfEdges
          Element => Mesh % Edges(t)
          IF(ASSOCIATED(Element % BoundaryInfo % Left)) THEN
            Element % BoundaryInfo % Left => &
                NewElements(Element % BoundaryInfo % Left % ElementIndex)
          END IF
          IF(ASSOCIATED(Element % BoundaryInfo % Right)) THEN
            Element % BoundaryInfo % Right => &
                NewElements(Element % BoundaryInfo % Right % ElementIndex)
          END IF
        END DO
        
        DEALLOCATE(Mesh % Elements)
        Mesh % Elements => NewElements
        Mesh % NumberOfBoundaryElements = nbc + NewCnt
      END IF
      
    END SUBROUTINE CreateIntersection2D
    
    
  END SUBROUTINE CreateIntersectionBCs

  SUBROUTINE TagBCsUsingRule(Model, Mesh)

    TYPE(Model_t) :: Model
    TYPE(Mesh_t) :: Mesh

    TYPE(Element_t), POINTER :: Element, Parent
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER :: i,j,k,n,m,t,t0
    INTEGER :: bc_ind, pSign, rSign, dim, BCsTagged, RuleInd
    TYPE(ValueList_t), POINTER :: BC
    REAL(KIND=dp) :: Coord(3), eps, val, r, rad, phi, phimin, phimax, RuleC
    LOGICAL :: Found, Hit, Parallel, CreateBCs, SplitBC, DoIt
    INTEGER, ALLOCATABLE :: EdgeConstraint(:)
    CHARACTER(:), ALLOCATABLE :: RuleStr
    CHARACTER(*), PARAMETER :: Caller = 'TagBCsUsingRule'
     

    Parallel = ( ParEnv % PEs > 1 )
    IF( Parallel ) THEN
      IF( ListGetLogical( Model % Simulation,'Single Mesh',Found ) ) THEN
        Parallel = .FALSE.
        CALL Info(Caller,'Working on single mesh, so reverting parallel mode to serial!',Level=8)
      END IF
    END IF
    
    ! We may need the rotor radius in defining certain BCs.
    DoIt = ListGetLogical( Model % Simulation,'Rotor Mode',Found) .AND. &
        .NOT. ListCheckPresent( Model % Simulation,'Rotor Radius')
    DoIt = DoIt .OR. ListGetLogical( Model % Simulation,'Determine Rotor Radius',Found)
    IF(DoIt) THEN
      IF( Parallel ) THEN
        CALL Fatal(Caller,'Cannot determine "Rotor Radius" yet in parallel!')
      ELSE        
        Rad = DetermineRotorRadius(Mesh)
        IF(Rad>0) THEN
          CALL ListAddConstReal(Model % Simulation,'Rotor Radius',Rad)
          WRITE(Message,'(A,ES14.6)') '"Rotor Radius" is found to be: ',Rad
          CALL Info(Caller,Message)
        ELSE
          CALL Fatal(Caller,'Could not determine "Rotor Radius", maybe there are not two pieces!?')
        END IF
      END IF
    END IF    
    
    ! Nothing to do with any boundary   
    CreateBCs = ListCheckPresentAnyBC(Model,'Boundary Create')
    IF(.NOT. ( ListCheckPresentAnyBC( Model,'Boundary Levelset') .OR. &
        ListCheckPresentAnyBC(Model,'Boundary Detect') .OR. CreateBCs ) ) RETURN 

    CALL Info(Caller,'Tagging BCs using simple geometric detection rules')

    IF(CreateBCs) THEN
      !IF(.NOT. ASSOCIATED(Mesh % Edges) ) THEN
        CALL FindMeshEdges2D(Mesh)
      !END IF
      ALLOCATE(EdgeConstraint(Mesh % NumberofEdges) )
      EdgeConstraint = 0
    END IF

    SplitBC = .FALSE.
    dim = Mesh % MeshDim
    BCsTagged = 0
    t0 = Mesh % NumberOfBulkElements

    n = 0
    DO t=1, Mesh % NumberOfBoundaryElements
      Element => Mesh % Elements(t0+t)      
      IF(Element % BoundaryInfo % Constraint == 0) n=n+1
    END DO
    CALL Info(Caller,'Number of unconstrained boundary elements: '//I2S(n))

    IF(.NOT. CreateBCs ) THEN
      m = n
      IF(Parallel) m = ParallelReduction(m)
      IF(m == 0) THEN
        CALL Warn(Caller,'Boundary detection requested but all boundaries already set!')
        RETURN
      END IF
    END IF
    
    n = Mesh % NumberOfNodes 

    DO bc_ind = 1, Model % NumberOfBCs
      BC => Model % BCs(bc_ind) % Values
      
      eps = ListGetCReal(BC,'Boundary Detect Epsilon',Found )
      IF(.NOT. Found) eps = 1.0e-6
    
      pSign = 0
      RuleInd = 0
      !PRINT *,'bc tag:',bc_ind, Model % BCs(bc_ind) % Tag
      
      ! Do we have a simple rule? 
      ! These rules have been designed such that particularly electrical machines
      ! that have rather constant set of BCs can be treated easily. 
      RuleStr = ListGetString(BC,'Boundary Detect',Found )
      IF(.NOT. Found) RuleStr = ListGetString(BC,'Boundary Create',Found )

      ! When we have rules "phimax" and "phimin" they may be augmented with inner and outer
      rSign = 0
      IF( LEN_TRIM(RuleStr) == 12 ) THEN
        IF( RuleStr(8:12) == 'inner' ) THEN
          Rad = ListGetCReal( CurrentModel % Simulation,'Rotor Radius',Found)
          rSign = -1
        ELSE IF(RuleStr(8:12) == 'outer') THEN
          Rad = ListGetCReal( CurrentModel % Simulation,'Rotor Radius',Found)
          rSign = 1
        END IF
      END IF

      IF( Found ) THEN
        SELECT CASE( RuleStr )
        CASE('xmin')
          RuleInd = 1
          RuleC = MINVAL( Mesh % Nodes % x )
          IF(Parallel) RuleC = ParallelReduction(RuleC,1)
        CASE('xmax')
          RuleInd = 1
          RuleC = MAXVAL( Mesh % Nodes % x )
          IF(Parallel) RuleC = ParallelReduction(RuleC,2)
        CASE('ymin')
          RuleInd = 2
          RuleC = MINVAL( Mesh % Nodes % y )
          IF(Parallel) RuleC = ParallelReduction(RuleC,1)
        CASE('ymax')
          RuleInd = 2
          RuleC = MAXVAL( Mesh % Nodes % y )
          IF(Parallel) RuleC = ParallelReduction(RuleC,2)
        CASE('zmin')
          RuleInd = 3
          RuleC = MINVAL( Mesh % Nodes % z )
          IF(Parallel) RuleC = ParallelReduction(RuleC,1)
        CASE('zmax')
          RuleInd = 3
          RuleC = MAXVAL( Mesh % Nodes % z )
          IF(Parallel) RuleC = ParallelReduction(RuleC,2)
        CASE('r inner')
          RuleInd = 4 
          RuleC = ListGetCReal( CurrentModel % Simulation,'Rotor Radius',Found)
          IF(.NOT. Found) CALL Fatal(Caller,'boundary detect "r inner" requires "Rotor Radius" be given!')
          pSign = -1
        CASE('r outer')
          RuleInd = 4
          RuleC = ListGetCReal( CurrentModel % Simulation,'Rotor Radius',Found)
          IF(.NOT. Found) CALL Fatal(Caller,'boundary detect "r outer" requires "Rotor Radius" be given!')
          pSign = 1
        CASE('phimax','phimax inner','phimax outer')
          ! The "inner" and "outer" rules allow to separate rotor and stator
          phimax = -2*PI
          phimin = 2*PI
          m = 0
          DO i=1,n
            ! Skip the nodes at exact origin
            r = SQRT(Mesh % Nodes % x(i)**2 + Mesh % Nodes % y(i)**2)
            IF(r < eps) CYCLE
            IF(rSign == -1 .AND. r > rad-eps ) CYCLE
            IF(rSign == 1 .AND. r < rad+eps ) CYCLE            
            phi = ATAN2(Mesh % Nodes % y(i), Mesh % Nodes % x(i) )
            phimax = MAX(phimax,phi)
            phimin = MIN(phimin,phi)
            m = m+1
          END DO
          IF(Parallel) THEN
            phimin = ParallelReduction(phimin,1)
            phimax = ParallelReduction(phimax,2)
          END IF
          ! There is a discontinuity at PI. Warn about it so we can code some more. 
          IF(phimax - phimin > PI ) CALL Fatal(Caller,'dPhi bigger than PI?')
          RuleC = phimax
          RuleInd = 5           
        CASE('phimin','phimin inner','phimin outer')
          phimax = -2*PI
          phimin = 2*PI
          m  = 0
          DO i=1,n
            r = SQRT(Mesh % Nodes % x(i)**2 + Mesh % Nodes % y(i)**2)
            IF(r < eps) CYCLE
            IF(rSign == -1 .AND. r > rad-eps ) CYCLE
            IF(rSign == 1 .AND. r < rad+eps ) CYCLE            
            phi = ATAN2(Mesh % Nodes % y(i), Mesh % Nodes % x(i) )
            phimax = MAX(phimax,phi)
            phimin = MIN(phimin,phi)
            m = m+1
          END DO
          IF(Parallel) THEN
            phimin = ParallelReduction(phimin,1)
            phimax = ParallelReduction(phimax,2)
          END IF
          IF(phimax - phimin > PI ) CALL Fatal(Caller,'dPhi bigger than PI?')
          RuleC = phimin
          RuleInd = 5          
        CASE('rmin')
          RuleInd = 6
          RuleC = HUGE(RuleC)
          DO i=1,n 
            RuleC = MIN(RuleC, Mesh % Nodes % x(i)**2+Mesh % Nodes % y(i)**2)
          END DO
          RuleC = SQRT(RuleC)
          IF(Parallel) RuleC = ParallelReduction(RuleC,1)
        CASE('rmax')
          RuleInd = 6
          RuleC = 0.0_dp
          DO i=1,n 
            RuleC = MAX(RuleC, Mesh % Nodes % x(i)**2+Mesh % Nodes % y(i)**2)
          END DO
          RuleC = SQRT(RuleC)
          IF(Parallel) RuleC = ParallelReduction(RuleC,2)
        CASE DEFAULT
          CALL Info(Caller,"Available rules: xmin, xmax, ymin, ymax, rmin, rmax, r inner, r outer'&
              //', phimin (inner/outer), phimax (inner/outer)",Level=3)
          CALL Fatal(Caller,'Uknown "Boundary Detect" method: '//TRIM(RuleStr))
        END SELECT
      ELSE          
        IF( .NOT. ListCheckPresent(BC,'Boundary Levelset') ) CYCLE
        ! Should we check the sign of the parent too?
        
        pSign = ListGetInteger(BC,'Boundary Levelset Parent Sign',Found )
        IF(Found) THEN
          IF(ABS(pSign) /= 1 ) THEN
            CALL Fatal(Caller,'"Boundary Levelset Parent Sign" should be either 1 or -1')
          END IF
        END IF        
      END IF

      IF( RuleInd == 3 .AND. Mesh % MeshDim < 3 ) THEN
        CALL Fatal(Caller,'Cannot use z-rules for 2D mesh!')
      END IF
      
      CALL Info(Caller,'Trying to tag elements to boundary: '//I2S(bc_ind),Level=20)
      BCsTagged = 0
      
      CALL TagElements()

      IF(Parallel) BCsTagged = ParallelReduction(BCsTagged)
      
      CALL Info(Caller,'Number of boundary elements "'//TRIM(RuleStr)//'" tagged to '&
          //I2S(bc_ind)//' is: '//I2S(BCsTagged),Level=7)
      IF( BCsTagged == 0 ) THEN
        CALL Fatal(Caller,'Could not find any boundary elements with rule!')
      END IF

    END DO

    IF( CreateBCs ) THEN
      CALL EdgesToBoundaryElements()
      IF(SplitBC) CALL Info(Caller,'Some of the boundaries was an internal one!',Level=10)
    END IF

    CALL Info(Caller,'Done creating additional BCs',Level=10)
    
    
  CONTAINS

    FUNCTION DiffPhi(Coord,Phi0) RESULT( val )
      REAL(KIND=dp) :: Coord(:)
      REAL(KIND=dp) :: Phi0, val
      INTEGER :: i

      IF(SQRT(SUM(Coord(1:2)**2)) < eps) THEN
        val = 0.0_dp
        RETURN
      END IF
      
      val = ATAN2(Coord(2),Coord(1)) - Phi0
      IF( val > PI ) THEN
        val = val - 2*PI
      ELSE IF( val < -PI ) THEN
        val = val + 2*PI
      END IF                 
    END FUNCTION DiffPhi

    
    SUBROUTINE TagElements()      
      INTEGER :: nc, np
      LOGICAL :: Hit
      REAL(KIND=dp) :: val
      INTEGER :: CandElems

      IF(CreateBCs) THEN
        CandElems = Mesh % NumberOfEdges
      ELSE
        CandElems = Mesh % NumberOfBoundaryElements
      END IF
      
      DO t=1, CandElems
        IF(CreateBCs) THEN
          Element => Mesh % Edges(t)
        ELSE
          Element => Mesh % Elements(t0+t)
          IF(Element % BoundaryInfo % Constraint > 0) CYCLE
        END IF
          
        ! Number of corners
        nc = Element % TYPE % ElementCode / 100
        Hit = .TRUE.
        
        DO i=1,nc
          j = Element % NodeIndexes(i)
          Coord(1) = Mesh % Nodes % x(j)
          Coord(2) = Mesh % Nodes % y(j)
          Coord(3) = Mesh % Nodes % z(j)

          SELECT CASE( RuleInd )

          CASE(1,2,3)
            val = Coord(RuleInd) - RuleC
          CASE(4,6)
            val = RuleC - SQRT(SUM(Coord(1:2)**2))                 
          CASE(5)
            val = DiffPhi(Coord,RuleC)
          CASE DEFAULT
            val = ListGetFunVec(BC,'Target Levelset',Coord(1:dim),dim)
          END SELECT
            
          IF(ABS(val) > eps) THEN
            Hit = .FALSE.            
            EXIT
          END IF
        END DO
        IF(.NOT. Hit) CYCLE

        ! We may additionally test inner/outer rule for the radius
        IF(rSign /= 0) THEN
          DO i=1,nc
            j = Element % NodeIndexes(i)
            Coord(1) = Mesh % Nodes % x(j)
            Coord(2) = Mesh % Nodes % y(j)
            Coord(3) = Mesh % Nodes % z(j)
            
            IF( SQRT(SUM(Coord(1:2)**2)) < eps ) CYCLE

            ! This should be negative for rotor, positive for stator
            val = SQRT(SUM(Coord(1:2)**2)) - Rad                
            IF(ABS(val) < eps) CYCLE

            IF( val*rSign < 0 ) THEN
              Hit = .FALSE.            
              EXIT
            END IF
          END DO
          IF(.NOT. Hit) CYCLE
        END IF
          
        
        IF(pSign /= 0) THEN
          IF( ASSOCIATED( Element % BoundaryInfo % Left ) .AND. &
              ASSOCIATED( Element % BoundaryInfo % Right ) ) SplitBC = .TRUE.

          DO k=1,2
            IF(k==1) THEN
              Parent => Element % BoundaryInfo % Left
            ELSE
              Parent => Element % BoundaryInfo % Right
            END IF

            IF(.NOT. ASSOCIATED(Parent)) CYCLE
            
            ! Number of corners, in 3D we must treat tets, pyramids, and wedges.
            np = Parent % TYPE % ElementCode / 100
            IF(np >= 5 .AND. np <= 7) np = np-1
            
            Hit = .TRUE.
            DO i=1,np
              j = Parent % NodeIndexes(i)
              IF(ANY(Element % NodeIndexes(1:nc) == j)) CYCLE
              
              ! Use the 1st non-boundary corner node to test the condition
              Coord(1) = Mesh % Nodes % x(j)
              Coord(2) = Mesh % Nodes % y(j)
              Coord(3) = Mesh % Nodes % z(j)
              
              SELECT CASE( RuleInd )
                
              CASE(4)
                val = RuleC**2 - SUM(Coord(1:2)**2)                 
              CASE DEFAULT
                val = ListGetFunVec(BC,'Target Levelset',Coord(1:dim),dim)
              END SELECT
              
              IF( val*pSign < 0.0_dp ) Hit = .FALSE.
              ! One node is representative
              EXIT
            END DO

            IF(Hit) EXIT
          END DO
          IF(.NOT. Hit) CYCLE
        END IF
        
        IF(CreateBCs) THEN
          EdgeConstraint(t) = Model % BCs(bc_ind) % Tag
        ELSE
          Element % BoundaryInfo % Constraint = Model % BCs(bc_ind) % Tag
        END IF
        BCsTagged = BCsTagged + 1
      END DO
    
    END SUBROUTINE TagElements   

    
    SUBROUTINE EdgesToBoundaryElements()
      
      INTEGER :: nbulk,nbc,nold,nadd,npar,parentcnt(0:2)
      TYPE(Element_t), POINTER :: NewElements(:),Element,Parent,Edge
      
      nbulk = Mesh % NumberOfBulkElements
      nbc = Mesh % NumberOfBoundaryElements
      nold = nbulk + nbc
      nadd = COUNT(EdgeConstraint > 0)

      IF(nadd == 0) THEN
        CALL Info(Caller,'No new boundary elements to add!',Level=6)
        RETURN
      END IF
        
      ALLOCATE( NewElements(nold + nadd ) )
      CALL Info(Caller,'Allocated for '//I2S(nadd)//' new boundary elements!',Level=6)
      
      NewElements(1:nold) = Mesh % Elements(1:nold)
      
      DO i=nbulk+1,nold
        Element => Mesh % Elements(i)        
        IF(.NOT. ASSOCIATED(Element % BoundaryInfo)) CYCLE
        
        Parent => Element % BoundaryInfo % Left
        IF(ASSOCIATED(Parent)) THEN
          NewElements(i) % BoundaryInfo % Left => NewElements(Parent % ElementIndex)
        END IF
        
        Parent => Element % BoundaryInfo % Right
        IF(ASSOCIATED(Parent)) THEN
          NewElements(i) % BoundaryInfo % Right => NewElements(Parent % ElementIndex)
        END IF
      END DO
      
      DEALLOCATE(Mesh % Elements)
      Mesh % Elements => NewElements
      Mesh % NumberOfBoundaryElements = nbc + nadd

      k = nold
      parentCnt = 0

      DO i=1,Mesh % NumberOfEdges
        j = EdgeConstraint(i)
        IF(j==0) CYCLE
        k = k+1

        Element => Mesh % Elements(k)
        Edge => Mesh % Edges(i)

        IF(.NOT. ASSOCIATED(Element)) THEN
          CALL Fatal(Caller,'Element not associated!?')
        END IF
        IF(.NOT. ASSOCIATED(Edge)) THEN
          CALL Fatal(Caller,'Edge not associated!?')
        END IF
          
        Element % TYPE => Edge % TYPE

        IF(.NOT. ASSOCIATED(Element % BoundaryInfo)) THEN
          ALLOCATE(Element % BoundaryInfo)
        END IF

        npar = 0
        IF(ASSOCIATED(Edge % BoundaryInfo) ) THEN
          Parent => Edge % BoundaryInfo % Left        
          IF(ASSOCIATED(Parent)) THEN
            Element % BoundaryInfo % Left => Mesh % Elements(Parent % ElementIndex)
            npar = npar + 1
          END IF
          Parent => Edge % BoundaryInfo % Right
          IF(ASSOCIATED(Parent)) THEN
            Element % BoundaryInfo % Right => Mesh % Elements(Parent % ElementIndex)
            npar = npar + 1
          END IF
        END IF
        Element % BoundaryInfo % Constraint = j

        parentCnt(npar) = parentCnt(npar) + 1
                
        n = Element % TYPE % NumberOfNodes
        ALLOCATE(Element % NodeIndexes(n))
        Element % NodeIndexes = Edge % NodeIndexes
      END DO
      
      DO i=0,2
        j = parentCnt(i)
        IF(j>0) CALL Info(Caller,'New boundary elements with '//I2S(i)//' parents: '//I2S(j),Level=6)
      END DO
      
    END SUBROUTINE EdgesToBoundaryElements
          
  END SUBROUTINE TagBCsUsingRule

  SUBROUTINE TagBodiesUsingCondition(Model, Mesh)

    TYPE(Model_t) :: Model
    TYPE(Mesh_t), TARGET :: Mesh

    TYPE(Element_t), POINTER :: Element
    INTEGER :: i,j,k,n,m,t
    INTEGER :: body_id
    REAL(KIND=dp) :: bodyCond(27)
    LOGICAL :: Found, Parallel, Sloppy, Conservative, Switch
    CHARACTER(:), ALLOCATABLE :: str
    CHARACTER(*), PARAMETER :: Caller = 'TagBodiesUsingCondition'

    
    ! Check that there is something to do, else exit
    IF(.NOT. ListCheckPrefix(Model % Simulation,'Body Define Condition') ) RETURN
    
    Parallel = ( ParEnv % PEs > 1 )
    IF( Parallel ) THEN
      IF( ListGetLogical( Model % Simulation,'Single Mesh',Found ) ) THEN
        Parallel = .FALSE.
        CALL Info(Caller,'Working on single mesh, so reverting parallel mode to serial!',Level=8)
      END IF
    END IF
    
    CALL Info(Caller,'Redefining bodies using geometric detection rules')
            
    Sloppy = ListGetLogical( Model % Simulation,'Body Define Sloppy', Found ) 
    Conservative = ListGetLogical( Model % Simulation,'Body Define Conservative', Found ) 

    ! If these are not set we cannot use ListGetReal later on...
    Model % Mesh => Mesh
    Model % Variables => Mesh % Variables
    
    DO i=1,100       
      str = 'Body Define Condition '//I2S(i)
      IF(.NOT. ListCheckPresent(Model % Simulation,str)) EXIT
      
      body_id = ListGetInteger(Model % Simulation,'Body Define Index '//I2S(i),UnfoundFatal=.TRUE.)
      
      k = 0
      DO t=1,Mesh % NumberOfBulkElements
        Element => Mesh % Elements(t)
        Model % CurrentElement => Element
        n = Element % Type % NumberOfNodes

        IF(Element % BodyId == body_id) CYCLE
        
        bodyCond(1:n) = ListGetReal( Model % Simulation, str, n, Element % NodeIndexes )
        m = COUNT(bodyCond(1:n) > 0)
        
        Switch = .FALSE.
        IF( Conservative ) THEN
          Switch = (m==n)
        ELSE IF( Sloppy ) THEN
          Switch = (m>0)
        ELSE
          Switch = (2*m>n)
        END IF

        IF(switch) THEN
          k = k+1
          Element % BodyId = body_id
        END IF
      END DO

      CALL Info(Caller,'Defining body index '//I2S(body_id)//' in '//I2S(k)//' elements')
      
    END DO
      
  END SUBROUTINE TagBodiesUsingCondition


END MODULE MeshTagging
