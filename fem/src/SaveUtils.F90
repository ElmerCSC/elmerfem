!*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation; either
! * version 2.1 of the License, or (at your option) any later version.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! * 
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library (in file ../LGPL-2.1); if not, write 
! * to the Free Software Foundation, Inc., 51 Franklin Street, 
! * Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen, Peter Råback
! *  Email:   elmeradm@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *****************************************************************************/
  
!> \ingroup ElmerLib
!> \{

!------------------------------------------------------------------------------
!>  Modules for generic stuff when saving results etc.
!------------------------------------------------------------------------------

MODULE SaveUtils

  USE Types
  USE SParIterGlobals
  USE MatrixAssembly
  USE Lists
  USE Messages
  USE MeshUtils, ONLY: GetLagrangeIndexes, CopyElementNodesFromMesh
  USE ElementDescription
  
  IMPLICIT NONE

CONTAINS


  ! Map element code of Elmer to the code used by VTK.
  !-----------------------------------------------------------------------------------
  FUNCTION Elmer2VtkElement( ElmerCode, SaveLinear ) RESULT ( VTKCode )
    INTEGER :: ElmerCode
    LOGICAL :: SaveLinear
    INTEGER :: VTKCode
    
    SELECT CASE (ElmerCode)
    CASE( 101 )
      VTKCode = 1
    CASE( 202 )
      VTKCode = 3
    CASE( 203 )
      VTKCode = 21
    CASE( 204 )
      VTKCode = 35  ! VTK_CUBIC_LINE but 68, VTK_LAGRANGE_CURVE, tested to work as well 
    CASE( 205, 206, 207, 208, 209 )
      VTKCode = 68
    CASE( 303 )
      VTKCode = 5
    CASE( 306 )
      VTKCode = 22
    CASE( 310, 315, 321, 328, 336, 345 )
      VTKCode = 69  ! VTK_LAGRANGE_TRIANGLE
    CASE( 404 )
      VTKCode = 9
    CASE( 408 )
      VTKCode = 23
    CASE( 409 )
      VTKCode = 28
    CASE( 416, 425, 436, 449, 464, 481 )
      VTKCode = 70  ! VTK_LAGRANGE_QUADRILATERAL
    CASE( 504 )
      VTKCode = 10
    CASE( 510 )
      VTKCode = 24
    CASE( 520, 535 )
      VTKCode = 71  ! VTK_LAGRANGE_TETRAHEDRON
    CASE( 605 )
      VTKCode = 14
    CASE( 613 )
      VTKCode = 27
    CASE( 706 )
      VTKCode = 13
    CASE( 715 ) 
      VTKCode = 26
    CASE( 718, 740, 775 ) 
      VTKCode = 73  ! VTK_LAGRANGE_WEDGE
    CASE( 808 )
      VTKCode = 12
    CASE( 820 )
      VTKCode = 25
    CASE( 827 )
      VTKCode = 29
    CASE( 864, 8125, 8216, 8343, 8512, 8729)
      VTKCode = 72  ! VTK_LAGRANGE_HEXAHEDRON
    CASE DEFAULT
      WRITE(Message,'(A,I0)') 'Not implemented for elementtype: ',ElmerCode
      CALL Fatal('Elmer2VtkElement',Message)

    END SELECT

    ! If requested, return the 1st order element corresponding to the higher order elements
    IF( SaveLinear ) THEN
      SELECT CASE (VTKCode)
      CASE( 21, 35 )
        VTKCode = 3
      CASE( 22, 69 )
        VTKCode = 5
      CASE( 23, 28, 70 )
        VTKCode = 9
      CASE( 24 )
        VTKCode = 10
      CASE( 27 )
        VTKCode = 14
      CASE( 26 )
        VTKCode = 13
      CASE( 25, 29 )
        VTKCode = 12
      END SELECT
    END IF

  END FUNCTION Elmer2VtkElement


  ! Map elemental node indexes of Elmer to the order used by VTK.
  !-----------------------------------------------------------------------------------
  SUBROUTINE Elmer2VtkIndexes( Element, DgElem, SaveLinear, NodeIndexes )
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: DgElem
    LOGICAL :: SaveLinear
    INTEGER :: NodeIndexes(:)

    TYPE(Element_t), POINTER :: Parent
    INTEGER, POINTER :: UseIndexes(:)
    INTEGER, TARGET :: BCIndexes(27)
    INTEGER :: ElmerCode, i,j,k,n,hits
    INTEGER, POINTER :: Order(:)
    INTEGER, TARGET, DIMENSION(16) :: &
        Order416 = (/1,2,3,4,5,6,7,8,10,9,12,11,13,14,16,15/)
    INTEGER, TARGET, DIMENSION(20) :: &
        Order820 = (/1,2,3,4,5,6,7,8,9,10,11,12,17,18,19,20,13,14,15,16/)
    INTEGER, TARGET, DIMENSION(27) :: &
        Order827 = (/1,2,3,4,5,6,7,8,9,10,11,12,17,18,19,20,13,14,15,16,24,22,21,23,25,26,27/)
    LOGICAL :: DoReorder


    ElmerCode = Element % Type % ElementCode


    IF( DGElem ) THEN
      UseIndexes => NULL()
      IF( ASSOCIATED( Element % DGIndexes ) ) THEN
        UseIndexes => Element % DGIndexes
      ELSE IF ( ASSOCIATED(Element % BoundaryInfo) ) THEN
        Parent => Element % BoundaryInfo % Left
        IF (.NOT.ASSOCIATED(Parent) ) THEN
          Parent => Element % BoundaryInfo % Right        
        END IF
        IF ( ASSOCIATED(Parent) ) THEN
          IF (.NOT. ASSOCIATED(Parent % DGIndexes) ) THEN
            ! This could happen if we have parents of parents i.e. the original element
            ! is a line element, has parents that are face elements, having parents being volume elements. 
            IF( ASSOCIATED( Parent % BoundaryInfo ) ) THEN
              IF( ASSOCIATED( Parent % BoundaryInfo % Left ) ) THEN
                Parent => Parent % BoundaryInfo % Left
              ELSE IF( ASSOCIATED( Parent % BoundaryInfo % Right ) ) THEN
                Parent => Parent % BoundaryInfo % Right
              END IF
            END IF
          END IF
          IF( ASSOCIATED( Parent % DGIndexes ) ) THEN
            n = Element % TYPE % NumberOfNodes 
            hits = 0
            DO j=1,n
              DO k=1,Parent % TYPE % NumberOfNodes
                IF(Element % NodeIndexes(j) == Parent % NodeIndexes(k)) THEN
                  BCIndexes(j) = Parent % DGIndexes(k) 
                  hits = hits + 1
                  EXIT
                END IF
              END DO
            END DO
            UseIndexes => BCIndexes
            IF( Hits < n ) THEN
              CALL Fatal('Elmer2VtkIndexes','Could not determine DG boundary indexes')
            END IF
          END IF
        END IF
      ENDIF

      IF(.NOT. ASSOCIATED( UseIndexes ) ) THEN
        CALL Warn('Elmer2VtkIndexes','Could not set DG indexes for boundary element!')        
        UseIndexes => Element % NodeIndexes
      END IF
    ELSE
      UseIndexes => Element % NodeIndexes
    END IF

    n = Element % TYPE % NumberOfNodes 


    ! Linear elements never require reordering 
    IF( .NOT. SaveLinear ) THEN
      SELECT CASE (ElmerCode)

      CASE( 416 )
        Order => Order416
        DoReOrder = .TRUE.

      CASE( 820 )
        Order => Order820
        DoReOrder = .TRUE.

      CASE( 827 ) 
        Order => Order827
        DoReOrder = .TRUE.

      CASE DEFAULT
        DoReorder = .FALSE.

      END SELECT
    ELSE
      DoReOrder = .FALSE.
    END IF

    IF( DoReorder ) THEN
      NodeIndexes(1:n) = UseIndexes( Order(1:n) )
    ELSE
      NodeIndexes(1:n) = UseIndexes(1:n)
    END IF

  END SUBROUTINE Elmer2VtkIndexes


  ! Map elemental node indexes of Elmer to the order used by Gmsh.
  !-----------------------------------------------------------------------------------
  SUBROUTINE ElmerToGmshIndex(Code,ElmerIndexes,GmshIndexes)

    INTEGER :: Code
    INTEGER :: ElmerIndexes(:),GmshIndexes(:)
    INTEGER :: i,n
    LOGICAL :: reorder, Visited = .FALSE.

    INTEGER, TARGET, SAVE :: order510(10),order613(13),order715(15),order820(20)
    INTEGER, POINTER :: order(:)

    SAVE Visited

    IF(.NOT. Visited ) THEN
      order510(:) = (/ 0,1,2,3,4,5,6,7,9,8 /)
      order613(:) = (/ 0,1,2,3,4,5,8,10,6,7,9,11,12 /)
      order715(:) = (/ 0,1,2,3,4,5,6,9,7,8,10,11,12,14,13 /)
      order820(:) = (/ 0,1,2,3,4,5,6,7,8,11,13,9,10,12,14,15,16,18,19,17 /)
      Visited = .TRUE.
    END IF

    reorder = .FALSE.

    SELECT CASE( Code )
      
    CASE (510)
      reorder = .TRUE.
      order => order510
      
    CASE (613)
      reorder = .TRUE.
      order => order613
      
    CASE (715)
      reorder = .TRUE.
      order => order715
      
    CASE (820)
      reorder = .TRUE.
      order => order820
     
    CASE DEFAULT
      
    END SELECT

    n = MOD(Code,100) 
    IF( reorder ) THEN
      DO i=1,n 
        GmshIndexes(order(i)+1) = ElmerIndexes(i)
      END DO
    ELSE
      GmshIndexes(1:n) = ElmerIndexes(1:n)      
    END IF


  END SUBROUTINE ElmerToGmshIndex

  
  
  ! Given different criteria for saving create a geometrical mask for elements
  ! and continuous numbering for the associated nodes.
  !------------------------------------------------------------------------------  
  SUBROUTINE GenerateSaveMask(Mesh,Params,Parallel,GroupId,SaveLinear,&
      NodePerm,ActiveElem,NumberOfGeomNodes,NumberOfElements, &
      ElemFirst,ElemLast)
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL :: Parallel
    INTEGER :: GroupId
    LOGICAL :: SaveLinear
    INTEGER, ALLOCATABLE :: NodePerm(:)
    LOGICAL, ALLOCATABLE :: ActiveElem(:)
    INTEGER :: NumberOfGeomNodes
    INTEGER :: NumberOfElements

    LOGICAL GotIt, SkipHalo, SaveOnlyHalo, GotMaskName, GotMaskCond, MaskExists, &
        SaveBulkOnly, SaveBoundariesOnly, GroupCollection, IsBoundaryElement, &
        IsHalo, Hit
    CHARACTER(:), ALLOCATABLE :: Str, MaskName
    REAL(KIND=dp), ALLOCATABLE :: MaskCond(:)
    INTEGER :: LeftIndex, RightIndex, ElemFirst, ElemLast, n, m, i, k, l
    TYPE(Variable_t), POINTER :: MaskVar
    INTEGER, POINTER :: MaskPerm(:), Indexes(:)
    TYPE(Element_t), POINTER :: Element, LeftElem, RightElem
    TYPE(Model_t), POINTER :: Model
    CHARACTER(*), PARAMETER :: Caller = 'GenerateSaveMask'
    
    Model => CurrentModel

    GroupCollection = ( GroupId > 0 ) 
    
    IF(.NOT. ALLOCATED( NodePerm ) ) THEN
      n = Mesh % NumberOfNodes
      CALL Info(Caller,'Allocating NodePerm of size: '//I2S(n),Level=15)
      ALLOCATE(NodePerm(n))
    END IF
    NodePerm = 0
    
    IF(.NOT. ALLOCATED(ActiveElem) ) THEN
      n = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      CALL Info(Caller,'Allocating ActiveElem of size: '//I2S(n),Level=15)
      ALLOCATE(ActiveElem(n) )
    END IF
    ActiveElem = .FALSE.
    
    IF( Parallel ) THEN
      SkipHalo = ListGetLogical( Params,'Skip Halo Elements', GotIt )
      IF(.NOT. GotIt) SkipHalo = .TRUE.
      SaveOnlyHalo = ListGetLogical( Params,'Save Halo Elements Only', GotIt )
    ELSE
      SkipHalo = .FALSE.
      SaveOnlyHalo = .FALSE.
    END IF

    GotMaskName = .FALSE.
    Str = ListGetString( Params,'Mask Variable',MaskExists)
    IF( MaskExists ) THEN
      MaskVar => VariableGet(Mesh % Variables,TRIM(Str),ThisOnly=.TRUE.)
      IF( ASSOCIATED(MaskVar)) MaskPerm => MaskVar % Perm
      MaskExists = ASSOCIATED(MaskPerm)
      IF( MaskExists ) THEN
        CALL Info(Caller,'Using > '// TRIM(Str) // ' < as mask variable')
      END IF
    ELSE
      ! Check if there is an additional mask name given
      GotIt = .FALSE.
      IF( Mesh % MeshDim == 2 ) THEN
        MaskName = ListGetString( Params,'2D Mask Name',GotIt)    
      ELSE IF( Mesh % MeshDim == 3 ) THEN  
        MaskName = ListGetString( Params,'3D Mask Name',GotIt)    
      END IF
      IF(.NOT. GotIt) MaskName = ListGetString( Params,'Mask Name',GotIt) 
      GotMaskName = GotIt
    END IF

    GotMaskCond = .FALSE.
    IF( .NOT. GotMaskName ) THEN
      MaskName = ListGetString( Params,'Mask Condition',GotMaskCond)
      IF( GotMaskCond ) THEN
        CALL Info(Caller,'Using mask condition: '//TRIM(MaskName),Level=8)
        n = Mesh % MaxElementNodes
        ALLOCATE( MaskCond(n) )
        MaskCond = 0.0_dp
      END IF
    END IF
    
    SaveBoundariesOnly = ListGetLogical( Params,'Save Boundaries Only',GotIt ) 
    IF( SaveBoundariesOnly ) CALL Info(Caller,'Saving only boundary elements!',Level=8)
    
    SaveBulkOnly = ListGetLogical( Params,'Save Bulk Only',GotIt ) 
    IF( SaveBulkOnly ) CALL Info(Caller,'Saving only bulk elements!',Level=8)
    
    NumberOfGeomNodes = Mesh % NumberOfNodes
    IF( MaskExists ) THEN
      NumberOfGeomNodes = COUNT( MaskPerm(1:NumberOfGeomNodes) > 0 ) 
      CALL Info(Caller,'Mask is positive for nodes: '//I2S(NumberOfGeomNodes),Level=15)
      IF( NumberOfGeomNodes == 0 ) THEN
        CALL Info(Caller,'Leaving early since mask is negative everywhere')
        RETURN
      END IF
    END IF

    NumberOfElements = 0
    ElemFirst = HUGE( ElemFirst )
    ElemLast = 0 
    
    ! Count the true number of elements and mark the 1st and last element
    !-----------------------------------------------------------------------
    DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

      IsBoundaryElement = ( i > Mesh % NumberOfBulkElements )
      
      IF( IsBoundaryElement ) THEN
        IF( SaveBulkOnly ) CYCLE
      ELSE
        IF( SaveBoundariesOnly ) CYCLE
      END IF

      Element => Mesh % Elements(i)
      CurrentModel % CurrentElement => Element

      IF( GroupCollection ) THEN
        IF( .NOT. IsBoundaryElement ) THEN
          IF( Element % BodyId /= GroupId ) CYCLE
        ELSE
          IF( .NOT. ASSOCIATED( Element % BoundaryInfo ) ) CYCLE
          IF( Element % BoundaryInfo % Constraint /= &
              GroupId - CurrentModel % NumberOfBodies ) CYCLE
        END IF
      END IF

      IF( Element % Type % ElementCode < 200 ) CYCLE          
      !IF (.NOT. IsBoundaryElement .AND. Element % BodyId < 1) CYCLE

      IF( SkipHalo .OR. SaveOnlyHalo ) THEN
        IF( IsBoundaryElement ) THEN
          IF( ASSOCIATED( Element % BoundaryInfo ) ) THEN
            LeftElem => Element % BoundaryInfo % Left
            IF( ASSOCIATED( LeftElem ) ) THEN
              LeftIndex = LeftElem % ElementIndex
              IF( LeftIndex > 0 ) THEN
                IF( Mesh % Elements(LeftIndex) % PartIndex /= ParEnv % MyPe ) LeftIndex = 0
              END IF
            ELSE
              LeftIndex = 0
            END IF
            RightElem => Element % BoundaryInfo % Right
            IF( ASSOCIATED( RightElem ) ) THEN
              RightIndex = RightElem % ElementIndex
              IF( RightIndex > 0 ) THEN
                IF( Mesh % Elements(RightIndex) % PartIndex /= ParEnv % MyPe ) RightIndex = 0
              END IF
            ELSE
              RightIndex = 0
            END IF
            IsHalo = ( LeftIndex == 0 .AND. RightIndex == 0 )
          ELSE
            IsHalo = .FALSE.
          END IF
        ELSE
          IsHalo = ( Element % PartIndex /= ParEnv % MyPe )
        END IF

        IF( IsHalo ) THEN
          IF( SkipHalo ) CYCLE
        ELSE
          IF( SaveOnlyHalo ) CYCLE
        END IF
      END IF


      IF( MaskExists ) THEN
        IF( ANY(MaskPerm(Element % NodeIndexes) <= 0) ) CYCLE
      END IF

      IF( GotMaskName ) THEN
        Hit = .FALSE.
        IF( .NOT. IsBoundaryElement ) THEN
          l = Element % BodyId
          k = ListGetInteger( CurrentModel % Bodies(l) % Values,'Body Force',GotIt)
          IF( GotIt ) THEN
            Hit = ListGetLogical( CurrentModel % BodyForces(k) % Values, TRIM(MaskName), GotIt)
          END  IF
          IF( .NOT. Hit ) THEN
            k = ListGetInteger( Model % Bodies(l) % Values,'Equation',GotIt)
            IF( GotIt ) THEN
              Hit = ListGetLogical( Model % Equations(k) % Values, TRIM(MaskName), GotIt)
            END IF
          END IF
        ELSE
          DO l=1, Model % NumberOfBCs
            IF ( Model % BCs(l) % Tag /= Element % BoundaryInfo % Constraint ) CYCLE
            Hit = ListGetLogical(Model % BCs(l) % Values, MaskName, GotIt ) 
            EXIT
          END DO
        END IF
        IF(.NOT. Hit ) CYCLE
      END IF ! GotMaskName

      IF( GotMaskCond ) THEN
        n = Element % TYPE % NumberOfNodes
        Indexes => Element % NodeIndexes
        GotIt = .FALSE.
        
        IF( .NOT. IsBoundaryElement ) THEN
          l = Element % BodyId
          k = ListGetInteger( Model % Bodies(l) % Values,'Body Force',GotIt)
          IF( GotIt ) THEN
            MaskCond(1:n) = ListGetReal( Model % BodyForces(k) % Values, TRIM(MaskName), &
                n, Indexes, GotIt)
          END  IF

          IF( .NOT. GotIt ) THEN
            k = ListGetInteger( Model % Bodies(l) % Values,'Equation',GotIt)
            IF( GotIt ) THEN
              MaskCond(1:n) = ListGetReal( Model % Equations(k) % Values, TRIM(MaskName), &
                  n, Indexes, GotIt)
            END IF
          END IF
        ELSE
          IF( ASSOCIATED( Element % BoundaryInfo ) ) THEN
            DO l=1, Model % NumberOfBCs
              IF ( Model % BCs(l) % Tag /= Element % BoundaryInfo % Constraint ) CYCLE
              MaskCond(1:n) = ListGetReal(Model % BCs(l) % Values, MaskName, &
                  n, Indexes, GotIt ) 
              EXIT
            END DO
          END IF
        END IF
        IF( .NOT. GotIt ) CYCLE
        IF( ANY(MaskCond(1:n) < 0.0_dp ) ) THEN
          CYCLE
        END IF
      END IF


      ActiveElem(i) = .TRUE.
      NumberOfElements = NumberOfElements + 1
      ElemFirst = MIN( ElemFirst, i )
      ElemLast = MAX( ElemLast, i )

      IF( SaveLinear ) THEN
        m = Element % TYPE % ElementCode / 100
        IF( m >= 5 .AND. m <= 7 ) m = m-1
        NodePerm( Element % NodeIndexes(1:m) ) = 1
      ELSE
        IF( MAXVAL( Element % NodeIndexes ) > SIZE( NodePerm ) ) THEN
          PRINT *,'too big:',SIZE(NodePerm), Element % NodeIndexes
        END IF
        NodePerm( Element % NodeIndexes ) = 1
      END IF

    END DO

    NumberOfGeomNodes = COUNT( NodePerm > 0 ) 
    IF( NumberOfElements == 0 ) THEN
      CALL Info(Caller,'No active elements forthis mask',Level=12)
    ELSE
      CALL Info(Caller,'Number of active elements '//I2S(NumberOfElements)//&
          ' out of '//I2S(Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements),Level=10)      
      CALL Info(Caller,'Number of geometry nodes '//I2S(NumberOfGeomNodes)//&
          ' out of '//I2S(Mesh % NumberOfNodes),Level=10)
    END IF
      
  END SUBROUTINE GenerateSaveMask
    
  
  ! Given the geometric permutation, create the dof permutation used in saving
  ! the different parts.
  !-----------------------------------------------------------------------------  
  SUBROUTINE GenerateSavePermutation(Mesh,DG,DN,LagN,SaveLinear,ActiveElem,NumberOfGeomNodes,&
      NoPermutation,NumberOfDofNodes,DgPerm,InvDgPerm,NodePerm,InvNodePerm)
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL :: DG, DN, SaveLinear
    INTEGER :: LagN
    LOGICAL, ALLOCATABLE :: ActiveElem(:)
    INTEGER :: NumberOfGeomNodes,NumberOfDofNodes
    LOGICAL :: NoPermutation
    INTEGER, ALLOCATABLE :: DgPerm(:),InvDgPerm(:),NodePerm(:),InvNodePerm(:)

    INTEGER, ALLOCATABLE :: BodyVisited(:)
    INTEGER :: i,j,k,l,m,n
    INTEGER :: Sweep
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER, ALLOCATABLE :: pIndexes(:)
    TYPE(Element_t), POINTER :: Element 
    TYPE(Model_t), POINTER :: Model
    CHARACTER(*), PARAMETER :: Caller = 'GenerateSavePermutation'

    
    Model => CurrentModel
            
    NumberOfDofNodes = 0
    IF( DG .OR. DN ) THEN
      NoPermutation = .FALSE.

      IF( LagN > 0 ) THEN
        CALL Fatal(Caller,'Cannot combine DG and higher order Lagrange elements!')
      END IF
      
      IF( DN ) THEN      
        CALL Info(Caller,'Saving results as discontinuous over bodies',Level=15)
        ALLOCATE( BodyVisited( Mesh % NumberOfNodes ) )
      ELSE
        CALL Info(Caller,'Saving results as discontinuous DG fields',Level=15)
      END IF

      IF( .NOT. ALLOCATED( DgPerm ) ) THEN
        k = 0
        DO i=1,Mesh % NumberOfBulkElements         
          Element => Mesh % Elements(i)
          k = k + Element % TYPE % NumberOfNodes
        END DO
        CALL Info(Caller,'Maximum number of dofs in DG: '//I2S(k),Level=12)
        ALLOCATE( DgPerm(k) )
      END IF
      DgPerm = 0
        
      DO Sweep=1,2
        l = 0
        IF( DG ) THEN
          DO i=1,Mesh % NumberOfBulkElements         
            IF( .NOT. ActiveElem(i) ) CYCLE
            Element => Mesh % Elements(i)
            NodeIndexes => Element % NodeIndexes

            IF( SaveLinear ) THEN
              m = Element % TYPE % ElementCode / 100
              IF( m >= 5 .AND. m <= 7 ) m = m-1
            ELSE
              m = Element % Type % NumberOfNodes
            END IF

            DO k=1,m
              IF( NodePerm( NodeIndexes(k) ) == 0 ) CYCLE
              l = l + 1
              IF( Sweep == 2 ) THEN
                InvNodePerm(l) = NodeIndexes(k)
                DgPerm( Element % DGIndexes(k) ) = l
                InvDgPerm(l) = Element % DGIndexes(k)
              END IF
            END DO
          END DO
        ELSE      
          DO i=1,Model % NumberOfBodies
            BodyVisited = 0
            DO j=1,Mesh % NumberOfBulkElements         
              IF(.NOT. ActiveElem(j) ) CYCLE
              Element => Mesh % Elements(j)
              IF( Element % BodyId /= i ) CYCLE
              NodeIndexes => Element % NodeIndexes

              IF( SaveLinear ) THEN
                m = Element % TYPE % ElementCode / 100
                IF( m >= 5 .AND. m <= 7 ) m = m-1
              ELSE
                m = Element % Type % NumberOfNodes
              END IF

              DO k=1,m
                IF( NodePerm( NodeIndexes(k) ) == 0 ) CYCLE
                IF( BodyVisited( NodeIndexes(k) ) > 0 ) THEN
                  DgPerm( Element % DGIndexes(k) ) = BodyVisited( NodeIndexes(k) )
                  CYCLE
                END IF
                l = l + 1
                BodyVisited(NodeIndexes(k)) = l
                IF( Sweep == 2 ) THEN
                  InvNodePerm(l) = NodeIndexes(k)
                  DgPerm( Element % DGIndexes(k) ) = l
                  InvDgPerm(l) = Element % DGIndexes(k)
                END IF
              END DO
            END DO
          END DO
        END IF

        IF( Sweep == 1 ) THEN
          CALL Info(Caller,'Independent dofs in discontinuous mesh: '//I2S(l),Level=10)
          NumberOfDofNodes = l
          IF(ALLOCATED(InvNodePerm)) DEALLOCATE( InvNodePerm )
          IF(ALLOCATED(InvDgPerm)) DEALLOCATE( InvDgPerm ) 
          ALLOCATE( InvNodePerm(l), InvDgPerm(l) ) 
          InvNodePerm = 0
          InvDgPerm = 0
        END IF
      END DO

      IF( DN ) DEALLOCATE( BodyVisited ) 

    ELSE IF( LagN > 0 ) THEN
      CALL Info(Caller,'Creating permutation for order '//I2S(LagN)//' Lagrange nodes!', Level=12)

      ! Calling without Element as argument returns the max. index value
      n = GetLagrangeIndexes( Mesh, LagN )
      ALLOCATE(DgPerm(n), pIndexes(n))
      DgPerm = 0        
      pIndexes = 0 
      
      ! Now call and then number the indexes!
      ! We use the same elemental subroutine to get the indexes as is done in the interpolation
      ! to avoid problems related to code inconsistency. There could be faster global ways too...
      DO i=1,Mesh % NumberOfBulkElements
        Element => Mesh % Elements(i)
        m = GetLagrangeIndexes( Mesh, LagN, Element, pIndexes )
        DgPerm(pIndexes(1:m)) = 1
      END DO

      m = 0
      DO i=1,n
        IF(DgPerm(i) > 0) THEN
          m = m+1
          DgPerm(i) = m
        END IF
      END DO

      ! Both the number of nodes and number of dofs will now follow the new higher order L-elements
      ! We will use no permutation for dofs or coordinates since we create a permutation-free temporal
      ! solution vectors and coordinates. 
      NumberOfDofNodes = m
      NumberOfGeomNodes = m
      NoPermutation = .TRUE.
      
      CALL Info(Caller,'Number of dofs for higher order Lagrange elements: '//I2S(m),Level=12)
    ELSE
      NoPermutation = ( NumberOfGeomNodes == Mesh % NumberOfNodes )    
      IF( NoPermutation ) THEN
        DEALLOCATE( NodePerm ) 
      ELSE
        CALL Info(Caller,'Not saving all nodes, creating permutation!',Level=12)
        IF( ALLOCATED( InvNodePerm ) ) DEALLOCATE( InvNodePerm ) 
        
        ALLOCATE( InvNodePerm( NumberOfGeomNodes ) ) 
        CALL Info(Caller,'Allocating InvNodePerm of size: '//I2S(NumberOfGeomNodes),Level=15)
        
        InvNodePerm = 0
        j = 0
        DO i=1,Mesh % NumberOfNodes
          IF( NodePerm(i) > 0 ) THEN
            j = j + 1       
            NodePerm(i) = j
            InvNodePerm(j) = i
          END IF
        END DO
        
      END IF
      NumberOfDofNodes = NumberOfGeomNodes 
    END IF
    
  END SUBROUTINE GenerateSavePermutation


  SUBROUTINE EvaluateVariableAtGivenPoint(No,Values,Mesh,Var,Var2,Var3,Element,LocalCoord,&
      LocalBasis,LocalNode,LocalDGNode,DoGrad,DoDiv,GotEigen,GotEdge)

    INTEGER :: No
    REAL(KIND=dp) :: Values(:)
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Variable_t), POINTER :: Var
    TYPE(Variable_t), POINTER, OPTIONAL :: Var2, Var3
    TYPE(Element_t), POINTER, OPTIONAL :: Element
    REAL(KIND=dp), OPTIONAL :: LocalCoord(3)
    INTEGER, OPTIONAL :: LocalNode
    INTEGER, OpTIONAL :: LocalDGNode
    REAL(KIND=dp), OPTIONAL, TARGET :: LocalBasis(:)
    LOGICAL, OPTIONAL :: DoGrad, DoDiv
    LOGICAL, OPTIONAL :: GotEigen, GotEdge
    
    LOGICAL :: Found, EdgeBasis, AVBasis, DgVar, IpVar, ElemVar, DoEigen, &
        PiolaVersion, PElem, NeedDerBasis, Stat, UseGivenNode, IsGrad, IsDiv, &
        IsEigen
    INTEGER :: i1,i2,ii,i,j,k,l,comps, n, nd, np, NoEigenValues, iMode
    INTEGER, TARGET :: DGIndexes(27), Indexes(100), DofIndexes(100), NodeIndex(1)
    INTEGER, POINTER :: pToIndexes(:)
    REAL(KIND=dp) :: u,v,w,detJ
    TYPE(Nodes_t), SAVE :: Nodes
    INTEGER :: lr
    REAL(KIND=dp), TARGET :: Basis(100),dBasisdx(100,3),NodeBasis(1)
    REAL(KIND=dp) :: WBasis(54,3),RotWBasis(54,3)
    REAL(KIND=dp), ALLOCATABLE, SAVE :: fdg(:), fip(:)
    REAL(KIND=dp), POINTER :: pToBasis(:)
    TYPE(Variable_t), POINTER :: pVar
    COMPLEX(KIND=dp), POINTER :: cValues(:)

    INTERFACE 
      SUBROUTINE Ip2DgFieldInElement( Mesh, Parent, nip, fip, np, fdg )
        USE Types
        TYPE(Mesh_t), POINTER :: Mesh
        TYPE(Element_t), POINTER :: Parent
        INTEGER :: nip, np
        REAL(KIND=dp) :: fip(:), fdg(:)
      END SUBROUTINE Ip2DgFieldInElement
    END INTERFACE

    IF(PRESENT(GotEdge)) GotEdge = .FALSE.
    IF(PRESENT(GotEigen)) GotEIgen = .FALSE.
        
    IF(.NOT. ASSOCIATED(Var)) RETURN
    IF(.NOT. ASSOCIATED(Var % Values)) RETURN

    ! Check that a vector field was not given by components 
    comps = 1
    IF(PRESENT(Var2)) THEN
      IF(ASSOCIATED(Var2)) THEN
        comps = 2
        IF( PRESENT(Var3)) THEN
          IF(ASSOCIATED(Var3)) THEN
            comps = 3
          END IF
        END IF
      END IF
    END IF

    IF(.NOT. PRESENT(Element) ) THEN
      CALL Fatal('EvaluteVariableAtGivenPoint','This routine most often really likes to have Element provided!')
    END IF
    
    EdgeBasis = ( Var % TYPE == variable_on_edges )

    DGVar = ( Var % TYPE == variable_on_nodes_on_elements ) 
    IpVar = ( Var % TYPE == variable_on_gauss_points )
    ElemVar = ( Var % TYPE == Variable_on_elements )       
    PiolaVersion = .FALSE.
    pElem = .FALSE.
    n = Element % TYPE % NumberOfNodes                  

    ! Do we have p-elements
    IF( ASSOCIATED( Var % Solver ) ) THEN
      nd = mGetElementDOFs( Indexes, Element, Var % Solver)
    ELSE
      nd = mGetElementDOFs( Indexes, Element )
    END IF
    IF(nd > n) pElem = isActivePElement(Element,Var % Solver)                              

    ! Elemental variable is a special case which is constant on the element
    ! This does not depend on any other things. 
    IF( ElemVar ) THEN
      l = Element % ElementIndex 
      IF( SIZE( Var % Perm ) >= l ) THEN
        l = Var % Perm(l)
      END IF
      IF( l > 0 ) THEN
        IF( Var % Dofs > 1 ) THEN
          DO ii=1,Var % Dofs
            Values(No+ii) = Var % Values(Var%Dofs*(l-1)+ii)
          END DO
        ELSE
          Values(No+1) = Var % Values(l)              
          IF( comps >= 2 ) Values(No+2) = Var % Values(l)
          IF( comps >= 3 ) Values(No+3) = Var % Values(l)
        END IF
      END IF
      No = No + MAX( Var % Dofs, comps )       
      RETURN
    END IF
    
    ! We may need to compute dBasisdx if we want derivative or divergence at the point. 
    NeedDerBasis = .FALSE.
    IsGrad = .FALSE.
    IsDiv = .FALSE.
    IF(PRESENT(DoGrad)) IsGrad = DoGrad
    IF(PRESENT(DoDiv)) IsDiv = DoDiv
    NeedDerBasis = IsGrad .OR. IsDiv
    
    pToBasis => NULL()
    pToIndexes => NULL()
    
    IsEigen = .FALSE.
    IF( PRESENT( GotEigen ) ) THEN
      IsEigen = ASSOCIATED( Var % EigenValues )
      IF(IsEigen) NoEigenValues = SIZE( Var % EigenValues )
      GotEigen = IsEigen
      IF( comps > 1 ) THEN
        CALL Warn('EvaluetVariableAtGivenPoint','Eigenmode cannot be given in components!')
        IsEigen = .FALSE.
      END IF
    END IF
    
    ! Given node is the quickest way to estimate the values at nodes.
    ! However, it only applies to nodal fields. 
    NodeIndex(1) = 0
    IF(.NOT. (EdgeBasis .OR. IpVar .OR. ElemVar .OR. NeedDerBasis .OR. pElem) ) THEN
      
      IF( DgVar ) THEN        
        IF(PRESENT(LocalDGNode)) NodeIndex(1) = LocalDGNode
        IF(NodeIndex(1) == 0 ) THEN
          PRINT *,'Find from base',PRESENT(Element), LocalNode
          IF( PRESENT(LocalNode)) THEN
            PToIndexes => PickDGIndexes(Element)
            DO i=1, n
              IF( Element % NodeIndexes(i) == LocalNode ) THEN
                NodeIndex(1) = pToIndexes(i)
                EXIT
              END IF
            END DO
          END IF
        END IF
      ELSE IF( PRESENT(LocalNode) ) THEN
        NodeIndex(1) = LocalNode
      END IF

      IF( NodeIndex(1) > 0 ) THEN
        pToIndexes => NodeIndex        
        NodeBasis(1) = 1.0_dp
        pToBasis => NodeBasis
        n = 1
        nd = 1        
      ELSE IF( PRESENT( LocalBasis ) ) THEN
        ! The 2nd quickest is to use existing local nodal basis functions
        PtoBasis => LocalBasis
        nd = n        
        IF( DgVar ) THEN
          PtoIndexes => PickDGIndexes(Element)
        ELSE
          PtoIndexes => Element % NodeIndexes 
        END IF
      END IF
    END IF
    

    IF(.NOT. ASSOCIATED(PtoIndexes) ) THEN
      IF(.NOT. PRESENT(LocalCoord) ) THEN
        CALL Fatal('EvaluteVariableAtGivenPoint',&
            'No recipe to evaluate variable without local coordinates: '//TRIM(Var % Name))
      END IF

      IF( EdgeBasis ) THEN
        IF( ListGetLogical(Var % Solver % Values, 'Quadratic Approximation', Found) ) THEN
          PiolaVersion = .TRUE.
        ELSE
          PiolaVersion = ListGetLogical(Var % Solver % Values,'Use Piola Transform', Found )   
        END IF
        np = n * Var % Solver % Def_Dofs(Element % type % ElementCode/100,Element % BodyId,1)
        AVBasis = (np > 0 )
        IF( PRESENT( GotEdge ) ) GotEdge = .TRUE.
      END IF
            
      u = LocalCoord(1)
      v = LocalCoord(2)
      w = LocalCoord(3)

      IF( EdgeBasis ) THEN
        CALL CopyElementNodesFromMesh(Nodes,Mesh,n,Element % NodeIndexes)
        stat = ElementInfo( Element, Nodes, u, v, w, &
            detJ, Basis, dBasisdx,  EdgeBasis = WBasis, &
            RotBasis = RotWBasis, USolver = Var % Solver)
      ELSE          
        IF( pElem ) THEN
          ! The standard element of SaveLine, SaveScalars etc. is most likely standard nodal element.
          ! If the user gives something else we may be trouble...
          CALL CopyElementNodesFromMesh(Nodes,Mesh,n,Element % NodeIndexes)
          CALL ConvertToPReference(Element % TYPE % ElementCode,u,v,w)            
          IF( NeedDerBasis ) THEN
            stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis, dBasisdx, &
                USolver = Var % Solver )
          ELSE
            stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis, &
                USolver = Var % Solver )
          END IF
        ELSE
          IF( NeedDerBasis ) THEN
            CALL CopyElementNodesFromMesh(Nodes,Mesh,n,Element % NodeIndexes)
            stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis, dBasisdx )
          ELSE
            CALL CopyElementNodesFromMesh(Nodes,Mesh,n,Element % NodeIndexes)          
            stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis )
          END IF
        END IF

        IF( DgVar ) THEN
          PtoIndexes => PickDGIndexes(Element)
        END IF
      END IF
      IF(.NOT. ASSOCIATED(pToIndexes)) PtoIndexes => Indexes
      IF(.NOT. ASSOCIATED(pToBasis)) PtoBasis => Basis
    END IF
    
    IF( EdgeBasis ) THEN
      DofIndexes(1:nd) = Var % Perm(PToIndexes(1:nd))
      IF( IsEigen ) THEN
        DO iMode = 1, NoEigenValues
          DO j=1,3          
            No = No + 1
            IF( ALL(DofIndexes(np+1:nd) > 0 ) ) THEN            
              cValues => Var % EigenVectors(iMode,:)
              Values(No) = SUM( WBasis(1:nd-np,j) * cValues(DofIndexes(np+1:nd)))
            ELSE
              Values(No) = 0.0_dp
            END IF
          END DO
        END DO
      ELSE
        DO j=1,3
          No = No + 1
          IF( ALL(DofIndexes(np+1:nd) > 0 ) ) THEN
            Values(No) = SUM( WBasis(1:nd-np,j) * Var % Values(DofIndexes(np+1:nd)))
          ELSE
            Values(No) = 0.0_dp
          END IF
        END DO
      END IF
      IF( AVBasis ) THEN
        No = No + 1
        Values(No) = 0.0_dp
        IF( ALL(DofIndexes(1:np) > 0 ) ) THEN
          Values(No) = SUM( Basis(1:np) * Var % Values(DofIndexes(1:np)))
        END IF
      END IF

    ELSE IF ( IpVar ) THEN
      l = 0 
      IF( SIZE( Var % Perm ) > Element % ElementIndex ) THEN 
        i1 = Var % Perm(Element % ElementIndex)
        i2 = Var % Perm(Element % ElementIndex+1)
        l = i2-i1
      END IF

      IF(l>0) THEN
        IF( .NOT. ALLOCATED(fip) .OR. SIZE(fip) < l ) THEN
          IF( ALLOCATED( fip ) ) DEALLOCATE( fip )
          ALLOCATE( fip(l) )
        END IF

        IF( .NOT. ALLOCATED(fdg) .OR. SIZE(fdg) < n ) THEN
          IF( ALLOCATED( fdg ) ) DEALLOCATE( fdg )
          ALLOCATE( fdg(n) )
        END IF

        DO ii=1,MAX(Var % Dofs,comps)
          IF( Var % Dofs > 1 ) THEN
            CONTINUE
          ELSE          
            IF( ii == 1 ) THEN
              pVar => Var
            ELSE IF( ii == 2 ) THEN
              pVar => Var2
            ELSE IF( ii == 3 ) THEN
              pVar => Var3
            END IF
            fip(1:l) = pVar % Values(i1+1:i2)
          END IF

          CALL Ip2DgFieldInElement( Mesh, Element, l, fip, n, fdg )              
          Values(No+ii) = SUM( PtoBasis(1:n) * fdg(1:n) )
        END DO
      END IF
      
      No = No + MAX(Var % Dofs, Comps )      
    ELSE
      IF(.NOT. ASSOCIATED(pToBasis)) THEN
        CALL Fatal('EvaluteVariableAtGivenPoint',&
            'pToBasis not associated for variable: '//TRIM(Var % Name))
      END IF
      IF(.NOT. ASSOCIATED(pToIndexes)) THEN
        CALL Fatal('EvaluteVariableAtGivenPoint',&
            'pToIndexes not associated for variable: '//TRIM(Var % Name))
      END IF
      
      IF( ASSOCIATED(Var % Perm) ) THEN
        DofIndexes(1:nd) = Var % Perm(PToIndexes(1:nd))
      ELSE
        DofIndexes(1:nd) = PtoIndexes(1:nd)
      END IF
      
      IF( ALL(Dofindexes(1:nd) > 0 ) ) THEN
        IF( IsGrad ) THEN
          IF( Var % Dofs /= 1 ) THEN
            CALL Fatal('EvaluteVariableAtGivenPoint','Gradient only possible for one dof!')
          END IF
          DO ii=1,3
            Values(No+ii) = SUM( dBasisdx(1:nd,ii) * Var % Values(DofIndexes(1:nd)))
          END DO
          No = No + 3
        ELSE IF( IsEigen ) THEN          
          DO iMode = 1, NoEigenValues
            cValues => Var % EIgenVectors(iMode,:)            
            DO ii=1,Var % DOfs              
              Values(No+ii) = Values(No+ii) + SUM( PtoBasis(1:nd) * &
                  cValues(Var%Dofs*(DofIndexes(1:nd)-1)+ii))
            END DO
          END DO
          No = No + Var % Dofs
        ELSE
          IF( Var % Dofs > 1 ) THEN
            DO ii=1,Var % Dofs
              Values(No+ii) = SUM( PtoBasis(1:nd) * &
                  Var % Values(Var%Dofs*(DofIndexes(1:nd)-1)+ii))
            END DO
          ELSE
            Values(No+1) = SUM(PToBasis(1:nd) * Var % Values(DofIndexes(1:nd)))            
            IF( comps >= 2 ) Values(No+2) = SUM(PToBasis(1:nd) * Var2 % Values(DofIndexes(1:nd)))
            IF( comps >= 3 ) Values(No+3) = SUM(PToBasis(1:nd) * Var3 % Values(DofIndexes(1:nd)))            
          END IF
          No = No + MAX( Var % Dofs, comps )
        END IF
      END IF
    END IF

#if 0
    ! Debugging code
    IF( .NOT. EdgeBasis .AND. ASSOCIATED(Var % Perm) .AND. No > 0 .AND. ASSOCIATED( PToIndexes ) ) THEN
      BLOCK
        REAL(KIND=dp) :: vmin, vmax
        INTEGER :: kmax
        PRINT *,'VarInfo:',TRIM(Var % Name), nd, n, pToIndexes(1:nd)

        kmax = 1
        DO i=1,nd
          IF( PtoBasis(i) > pToBasis(kmax)) kmax = i
        END DO

        IF( ASSOCIATED(pToIndexes) ) THEN
          DofIndexes(1:nd) = Var % Perm(PToIndexes(1:nd))
        ELSE
          DofIndexes(1:nd) = PtoIndexes(1:nd)
        END IF
          
        vmin = MINVAL(Var % Values(DofIndexes(1:nd)),DofIndexes(1:nd)>0)
        vmax = MAXVAL(Var % Values(DofIndexes(1:nd)),DofIndexes(1:nd)>0)

        PRINT *,'PtoIndexes:',PtoIndexes(1:nd)
        PRINT *,'PtoBasis',SUM(PToBasis(1:nd)),PToBasis(1:nd)
        PRINT *,'ElemRange:',vmin< Values(1) .AND. vmax > Values(1), &
            vmin, vmax, Var % Values(Var % Perm(PToIndexes(kmax))), Values(1:No)
      END BLOCK
    END IF
#endif
      
  CONTAINS

    FUNCTION PickDgIndexes(Element) RESULT ( PToInds) 
      TYPE(Element_t), POINTER :: Element
      INTEGER, POINTER :: PtoInds(:)

      TYPE(Element_t), POINTER :: Parent
      INTEGER :: i,j,lr
            
      IF( ASSOCIATED( Element % DgIndexes ) ) THEN
        PToInds => Element % DGIndexes
      ELSE IF( ASSOCIATED( Element % BoundaryInfo ) ) THEN
        DO lr=1,2
          IF(lr==1) THEN
            Parent => Element % BoundaryInfo % Left
          ELSE
            Parent => Element % BoundaryInfo % Right
          END IF
          IF(.NOT. ASSOCIATED( Parent ) ) CYCLE
          IF(.NOT. ASSOCIATED( Parent % DGIndexes ) ) CYCLE
          IF( ALL( Var % Perm( Parent % DGIndexes ) /= 0) ) THEN                  
            DO i=1,Element % TYPE % NumberOfNodes
              DO j=1,Parent % TYPE % NumberOfNodes
                IF( Element % NodeIndexes(i) == Parent % NodeIndexes(j) ) THEN
                  DGIndexes(i) = Parent % DGIndexes(j)
                  EXIT
                END IF
              END DO
            END DO
            EXIT
          END IF
        END DO
        PtoInds => DGIndexes
      END IF
    END FUNCTION PickDgIndexes
    
  END SUBROUTINE EvaluateVariableAtGivenPoint
    
END MODULE SaveUtils
  
