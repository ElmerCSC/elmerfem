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
  USE ElementUtils, ONLY: FindParentUVW
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
    INTEGER :: ElmerCode, i,j,k,n,hits,right
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
        n = Element % TYPE % NumberOfNodes 
                
        DO right=0,1
          hits = 0
          IF(right==0) THEN 
            Parent => Element % BoundaryInfo % Left
          ELSE
            Parent => Element % BoundaryInfo % Right        
          END IF
          IF(.NOT. ASSOCIATED(Parent)) CYCLE
          
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
            DO j=1,n
              DO k=1,Parent % TYPE % NumberOfNodes
                IF(Element % NodeIndexes(j) == Parent % NodeIndexes(k)) THEN
                  BCIndexes(j) = Parent % DGIndexes(k) 
                  hits = hits + 1
                  EXIT
                END IF
              END DO
            END DO            
          END IF

          
          IF(Hits == n ) THEN
            UseIndexes => BCIndexes
            EXIT
          END IF
        END DO
          
        IF( Hits < n ) THEN
          PRINT *,'Element:',n, Element % TYPE % ElementCode, Element % NodeIndexes
          PRINT *,'Parent:',Hits,Parent % TYPE % ElementCode, Parent % NodeIndexes
          CALL Fatal('Elmer2VtkIndexes','Could not determine DG boundary indexes')
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


   !> Defines and potentially creates output directory.
   !> The output directory may given in different ways, and even be part of the
   !> filename, or be relative to home directory. We try to parse every possible
   !> scenario here that user might have in mind.
   !-----------------------------------------------------------------------------
  SUBROUTINE SolverOutputDirectory( Solver, Filename, OutputDirectory, &
      MakeDir, UseMeshDir  )

    USE ModelDescription 
    
    TYPE(Solver_t) :: Solver
    LOGICAL, OPTIONAL :: MakeDir, UseMeshDir
    CHARACTER(*) :: Filename
    CHARACTER(:), ALLOCATABLE :: OutputDirectory

    LOGICAL :: Found, AbsPathInName, DoDir, PartitioningSubDir
    INTEGER :: nd, nf, n
    CHARACTER(LEN=MAX_NAME_LEN) :: Str

    IF( PRESENT( MakeDir ) ) THEN
      DoDir = MakeDir
    ELSE
      DoDir = ( Solver % TimesVisited == 0 ) .AND. ( ParEnv % MyPe == 0 )
    END IF

    ! Output directory is obtained in order
    ! 1) solver section
    ! 2) simulation section
    ! 3) header section
    OutputDirectory = ListGetString( Solver % Values,'Output Directory',Found) 
    IF(.NOT. Found) OutputDirectory = ListGetString( CurrentModel % Simulation,&
        'Output Directory',Found) 

    IF(.NOT. Found) OutputDirectory = TRIM(OutputPath)          
    nd = LEN_TRIM(OutputDirectory)

    ! If the path is just working directory then that is not an excude
    ! to not use the mesh name, or directory that comes with the filename 
    IF(.NOT. Found .AND. nd == 1 .AND. OutputDirectory(1:1)=='.') nd = 0

    ! If requested by the optional parameter use the mesh directory when
    ! no results directory given. This is an old convection used in some solvers. 
    IF( nd == 0 .AND. PRESENT( UseMeshDir ) ) THEN
      IF( UseMeshDir ) THEN
        OutputDirectory = TRIM(CurrentModel % Mesh % Name)
        nd = LEN_TRIM(OutputDirectory)       
      END IF
    END IF

    ! Use may have given part or all of the path in the filename.
    ! This is not preferred, but we cannot trust the user.
    nf = LEN_TRIM(Filename)        
    n = INDEX(Filename(1:nf),'/')
    AbsPathInName = INDEX(FileName,':')>0 .OR. (Filename(1:1)=='/') &
        .OR. (Filename(1:1)==Backslash)

    IF( nd > 0 .AND. .NOT. AbsPathInName ) THEN
      ! Check that we have not given the path relative to home directory
      ! because the code does not understand the meaning of tilde.
      IF(nd>=2) THEN
        IF( OutputDirectory(1:2) == '~/') THEN
          CALL get_environment_variable('HOME',Str)
          OutputDirectory = TRIM(Str)//'/'//OutputDirectory(3:nd)
          nd = LEN_TRIM(OutputDirectory)
        END IF
      END IF
      ! To be on the safe side create the directory. If it already exists no harm done.
      ! Note that only one directory may be created. Hence if there is a path with many subdirectories
      ! that will be a problem. Fortran does not have a standard ENQUIRE for directories hence
      ! we just try to make it. 
      IF( DoDir ) THEN
        CALL Info('SolverOutputDirectory','Creating directory: '//TRIM(OutputDirectory(1:nd)),Level=8)
        CALL MakeDirectory( OutputDirectory(1:nd) // CHAR(0) )      
      END IF
    END IF

    ! In this case the filename includes also path and we remove it from there and
    ! add it to the directory. 
    IF( n > 2 ) THEN    
      CALL Info('SolverOutputDirectory','Parcing path from filename: '//TRIM(Filename(1:n)),Level=10)
      IF( AbsPathInName .OR. nd == 0) THEN
        ! If the path is absolute then it overruns the given path!
        OutputDirectory = Filename(1:n-1)
        nd = n-1
      ELSE
        ! If path is relative we add it to the OutputDirectory and take it away from Filename
        OutputDirectory = OutputDirectory(1:nd)//'/'//Filename(1:n-1)        
        nd = nd + n 
      END IF
      Filename = Filename(n+1:nf)      

      IF( DoDir ) THEN
        CALL Info('SolverOutputDirectory','Creating directory: '//TRIM(OutputDirectory(1:nd)),Level=8)
        CALL MakeDirectory( OutputDirectory(1:nd) // CHAR(0) )
      END IF
    END IF

    ! Finally, on request save each partitioning to different directory.
    PartitioningSubDir = ListGetLogical( Solver % Values,'Output Partitioning Directory',Found)
    IF(.NOT. Found ) THEN
      PartitioningSubDir = ListGetLogical( CurrentModel % Simulation,'Output Partitioning Directory',Found)
    END IF
    IF( PartitioningSubDir ) THEN
      OutputDirectory = TRIM(OutputDirectory)//'/np'//I2S(ParEnv % PEs)
      nd = LEN_TRIM(OutputDirectory)             
      IF( DoDir ) THEN
        CALL Info('SolverOutputDirectory','Creating directory: '//TRIM(OutputDirectory(1:nd)),Level=8)
        CALL MakeDirectory( OutputDirectory(1:nd) // CHAR(0) )
      END IF
    END IF

  END SUBROUTINE SolverOutputDirectory
  !-----------------------------------------------------------------------------


  !------------------------------------------------------------------------------
  !> Saves results in ascii format understood by the pre-/postprocessing software Gmsh.
  !------------------------------------------------------------------------------
  SUBROUTINE SaveGmshOutput( Model,Solver,dt,Transient )
    !------------------------------------------------------------------------------
    USE Types
    USE Lists
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    TYPE(Model_t) :: Model
    REAL(KIND=dp) :: dt
    LOGICAL :: Transient
    !------------------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------------------
    TYPE(Element_t),POINTER :: Element
    INTEGER, POINTER :: Perm(:)
    REAL(KIND=dp), POINTER :: Values(:),Values2(:),Values3(:)
    REAL(KIND=dp) :: Vector(3), Time
    COMPLEX(KIND=dp), POINTER :: CValues(:)
    TYPE(Variable_t), POINTER :: Solution, TimeVariable
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Mesh_t), POINTER :: Mesh

    LOGICAL :: Found, GotField, FileAppend, AlterTopology, MaskExists
    LOGICAL :: EigenAnalysis = .FALSE., EigenActive, ComponentVector, Parallel

    INTEGER :: VisitedTimes = 0, ExtCount
    INTEGER :: i,j,k,l,m,n,nsize,dim,dofs,ElmerCode, GmshCode,body_id, Vari, Rank, truedim
    INTEGER :: Tag, NumberOfAllElements, BCOffSet
    INTEGER, PARAMETER :: MaxElemCode = 827
    INTEGER :: ElmerToGmshType(MaxElemCode), GmshToElmerType(21), &
        ElmerIndexes(27), GmshIndexes(27) 
    INTEGER, POINTER :: NodeIndexes(:)

    INTEGER, ALLOCATABLE :: NodePerm(:),DgPerm(:)
    INTEGER, ALLOCATABLE, TARGET :: InvDgPerm(:), InvNodePerm(:)
    LOGICAL, ALLOCATABLE :: ActiveElem(:)
    LOGICAL :: NoPermutation, Numbering 
    INTEGER :: NumberOfGeomNodes, NumberOfDofNodes,NumberOfElements, ElemFirst, ElemLast,bc_id
    INTEGER, POINTER :: InvFieldPerm(:), DGInvPerm(:)

    INTEGER, PARAMETER :: LENGTH = 1024
    CHARACTER(LEN=LENGTH) :: Txt, FieldName, CompName
    CHARACTER(MAX_NAME_LEN) :: OutputFile
    CHARACTER(:), ALLOCATABLE :: OutputDirectory
    INTEGER :: GmshUnit
    CHARACTER(*), PARAMETER :: Caller = 'SaveGmshOutput'

    SAVE VisitedTimes

    !------------------------------------------------------------------------------

    CALL Info(Caller,'Saving results in Gmsh format')

    Mesh => Model % Mesh
    Params => Solver % Values
    Parallel = ( ParEnv % PEs > 1 )

    ExtCount = ListGetInteger( Params,'Output Count',Found)
    IF( Found ) THEN
      VisitedTimes = ExtCount
    ELSE
      VisitedTimes = VisitedTimes + 1
    END IF

    Numbering = ListGetLogical( Params,'Filename Numbering',Found ) 
    IF(.NOT. Found) Numbering = .TRUE.
    
    GmshToElmerType = (/ 202, 303, 404, 504, 808, 706, 605, 203, 306, 409, &
        510, 827, 0, 0, 101, 408, 820, 715, 613, 0, 310 /)
    ElmerToGmshType = 0

    DO i=1,SIZE(GmshToElmerType)
      j = GmshToElmerType(i)
      IF( j > 0 ) ElmerToGmshType(j) = i
    END DO

    EigenAnalysis = ListGetLogical( Params, 'Eigen Analysis', Found )
    FileAppend = ListGetLogical( Params,'File Append',Found)
    IF(.NOT. Found) FileAppend = .TRUE.
    AlterTopology = ListGetLogical( Params,'Alter Topology',Found)

    OutputFile = ListGetString( Params, 'Output File Name', Found )
    IF( Found ) THEN
      IF(INDEX(OutputFile,'.') == 0) WRITE( OutputFile,'(A,A)') TRIM(OutputFile),".msh"
    ELSE
      OutputFile = 'Output.msh'
    END IF

    CALL SolverOutputDirectory( Solver, OutputFile, OutputDirectory, UseMeshDir = .TRUE. )
    OutputFile = TRIM(OutputDirectory)// '/' //TRIM(OutputFile)

    !------------------------------------------------------------------------------
    ! Initialize stuff for masked saving
    !------------------------------------------------------------------------------
    CALL GenerateSaveMask(Mesh,Params,Parallel,0,.FALSE.,&
        NodePerm,ActiveElem,NumberOfGeomNodes,NumberOfElements,&
        ElemFirst,ElemLast)

    IF( ParEnv % PEs > 1 ) THEN
      IF( NumberOfElements == 0 ) THEN
        CALL Info(Caller,'Nothing to save in partition: '//TRIM(I2S(ParEnv % MyPe)),Level=8)
        RETURN
      ELSE
        OutputFile = TRIM(OutputFile)//'_'//I2S(ParEnv % PEs)//'np'//I2S(ParEnv % MyPe+1)
      END IF
    ELSE
      IF( NumberOfElements == 0 ) THEN
        CALL Warn(Caller,'Notging to save, this is suspicious')
        RETURN      
      END IF
    END IF

    CALL GenerateSavePermutation(Mesh,.FALSE.,.FALSE.,0,.FALSE.,&
        ActiveElem,NumberOfGeomNodes,NoPermutation,NumberOfDofNodes,&
        DgPerm,InvDgPerm,NodePerm,InvNodePerm)

    InvFieldPerm => InvNodePerm

    dim = CoordinateSystemDimension()
    IF( VisitedTimes > 1 ) THEN
      IF( AlterTopology ) THEN
        IF( Numbering ) THEN
          OutputFile = NextFreeFilename( OutputFile )
        END IF        
        CALL Info(Caller,'Writing mesh and data to a new file: '//TRIM(OutputFile))
      ELSE IF( FileAppend ) THEN      
        CALL Info(Caller,'Appending data to the same file: '//TRIM(OutputFile))
        OPEN(NEWUNIT=GmshUnit, FILE=OutputFile, POSITION='APPEND' )      
        GOTO 10
      ELSE
        IF( Numbering ) THEN
          OutputFile = NextFreeFilename( OutputFile )
        END IF        
        CALL Info(Caller,'Writing data to a new file: '//TRIM(OutputFile))
        OPEN(NEWUNIT=GmshUnit, FILE=OutputFile )
        WRITE(GmshUnit,'(A)') '$MeshFormat'
        WRITE(GmshUnit,'(A)') '2.0 0 8'
        WRITE(GmshUnit,'(A)') '$EndMeshFormat'          
        GOTO 10    
      END IF
    END IF


    ! Save the header
    !-------------------------------------------------
    CALL Info(Caller,'Saving results to file: '//TRIM(OutputFile))
    OPEN(NEWUNIT=GmshUnit, FILE=OutputFile )

    WRITE(GmshUnit,'(A)') '$MeshFormat'
    WRITE(GmshUnit,'(A)') '2.0 0 8'
    WRITE(GmshUnit,'(A)') '$EndMeshFormat'    


    ! Save the mesh nodes
    !-------------------------------------------------
    CALL Info(Caller,'Writing the mesh nodes')
    CALL WriteGmshNodes()

    ! Save the mesh elements
    !-------------------------------------------------
    CALL Info(Caller,'Writing the mesh elements')
    CALL WriteGmshElements() 

    ! With a mask the list of physical entities should be checked
    !-------------------------------------------------------------
    IF(.NOT. MaskExists ) THEN
      !    CALL WritePhysicalNames() 
    END IF

10  CONTINUE

    CALL Info(Caller,'Writing the nodal data')
    CALL WriteGmshData()

    IF(.FALSE.) THEN
      WRITE(GmshUnit,'(A)') '$ElementData'
      WRITE(GmshUnit,'(A)') '$EndElementData'
    END IF

    IF(.FALSE.) THEN
      WRITE(GmshUnit,'(A)') '$ElementNodeData'
      WRITE(GmshUnit,'(A)') '$EndElementNodeData'
    END IF

    CLOSE(GmshUnit)

    IF(ALLOCATED(DgPerm)) DEALLOCATE(DgPerm)
    IF(ALLOCATED(InvDgPerm)) DEALLOCATE(InvDgPerm)
    IF(ALLOCATED(NodePerm)) DEALLOCATE(NodePerm)
    IF(ALLOCATED(InvNodePerm)) DEALLOCATE(InvNodePerm)
    IF(ALLOCATED(ActiveElem)) DEALLOCATE(ActiveElem)


    CALL Info(Caller,'Gmsh output complete')

  CONTAINS

    SUBROUTINE WriteGmshNodes()

      nsize = NumberOfGeomNodes

      WRITE(GmshUnit,'(A)') '$Nodes'
      WRITE(GmshUnit,'(I8)') nsize
      DO i = 1, nsize
        IF( NoPermutation ) THEN
          j = i 
        ELSE
          j = InvNodePerm(i)
        END IF

        IF( dim == 3 ) THEN
          WRITE(GmshUnit,'(I8,3ES16.7E3)') i,Mesh % Nodes % x(j),Mesh % Nodes % y(j), Mesh % Nodes % z(j)
        ELSE
          WRITE(GmshUnit,'(I8,2ES16.7E3,A)') i,Mesh % Nodes % x(j),Mesh % Nodes % y(j),' 0.0' 
        END IF
      END DO
      WRITE(GmshUnit,'(A)') '$EndNodes'
    END SUBROUTINE WriteGmshNodes


    SUBROUTINE WriteGmshElements()

      nsize = NumberOfElements 

      BCOffSet = 100
      DO WHILE( BCOffset <= Model % NumberOfBodies ) 
        BCOffset = 10 * BCOffset
      END DO

      WRITE(GmshUnit,'(A)') '$Elements'
      WRITE(GmshUnit,'(I8)') nsize

      l = 0
      DO i = ElemFirst, ElemLast
        IF(.NOT. ActiveElem(i) ) CYCLE

        l = l + 1
        Element => Mesh % Elements(i)
        ElmerCode = Element % TYPE % ElementCode

        n = Element % Type % NumberOfNodes
        IF( NoPermutation ) THEN
          ElmerIndexes(1:n) = Element % NodeIndexes(1:n)
        ELSE
          ElmerIndexes(1:n) = NodePerm(Element % NodeIndexes(1:n))
        END IF

        GmshCode = ElmerToGmshType(ElmerCode)
        IF( GmshCode == 0 ) THEN
          CALL Warn(Caller,'Gmsh element index not found!')
          CYCLE
        END IF

        IF( i <= Model % NumberOfBulkElements ) THEN
          Tag = Element % BodyId
        ELSE
          DO bc_id=1,Model % NumberOfBCs
            IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc_id) % Tag ) EXIT
          END DO
          Tag = bc_id + BCOffset
        END IF

        WRITE(GmshUnit,'(I8,I3,I3,I5,I5)',ADVANCE='NO') l,GmshCode,2,Tag,Tag
        k = MOD(ElmerCode,100)

        CALL ElmerToGmshIndex(ElmerCode,ElmerIndexes,GmshIndexes)

        DO j=1,k-1
          WRITE(GmshUnit,'(I8)',ADVANCE='NO') GmshIndexes(j)
        END DO
        WRITE(GmshUnit,'(I8)') GmshIndexes(k)
      END DO
      WRITE(GmshUnit,'(A)') '$EndElements'
    END SUBROUTINE WriteGmshElements


    SUBROUTINE WritePhysicalNames()
      CALL Info(Caller,'Writing the physical entity names')
      nsize = Model % NumberOfBodies + Model % NumberOfBCs
      WRITE(GmshUnit,'(A)') '$PhysicalNames'
      WRITE(GmshUnit,'(I8)') nsize
      DO i=1,Model % NumberOfBodies 
        Txt = ListGetString( Model % Bodies(i) % Values,'Name',Found)
        IF( Found ) THEN
          WRITE(GmshUnit,'(I8,A)') i,'"'//TRIM(Txt)//'"'
        ELSE
          WRITE(GmshUnit,'(I8,A,I0,A)') i,'"Body ',i,'"'       
        END IF
      END DO
      DO i=1,Model % NumberOfBCs
        Txt = ListGetString( Model % BCs(i) % Values,'Name',Found)
        IF( Found ) THEN
          WRITE(GmshUnit,'(I8,A)') i+BCOffset,'"'//TRIM(Txt)//'"'
        ELSE
          WRITE(GmshUnit,'(I8,A,I0,A)') i+BCOffset,'"Boundary Condition ',i,'"'               
        END IF
      END DO
      WRITE(GmshUnit,'(A)') '$EndPhysicalNames'
    END SUBROUTINE WritePhysicalNames


    ! In case of DG fields we average the fields on-the-fly to nodes.
    !----------------------------------------------------------------
    SUBROUTINE CreateTemporalNodalField(Mesh,Solution,Revert)
      TYPE(Mesh_t) :: Mesh
      TYPE(Variable_t) :: Solution
      LOGICAL, OPTIONAL :: revert

      REAL(KIND=dp), POINTER :: NodalVals(:), TmpVals(:)
      INTEGER, POINTER :: NodalPerm(:), TmpPerm(:), NodalCnt(:)
      INTEGER :: i,j,k,l,n,t,dofs,ElemFam

      SAVE NodalPerm, NodalVals, NodalCnt, TmpPerm, TmpVals

      IF( PRESENT( Revert ) ) THEN
        IF( Revert ) THEN
          DEALLOCATE( NodalVals, NodalPerm, NodalCnt )
          Solution % Perm => TmpPerm
          Solution % Values => TmpVals
          RETURN
        END IF
      END IF

      dofs = Solution % dofs

      n = Mesh % NumberOfNodes
      ALLOCATE( NodalPerm(n), NodalCnt(n), NodalVals(n*dofs) ) 
      NodalPerm = 0
      NodalCnt = 0
      NodalVals = 0.0_dp

      DO t=1,Mesh % NumberOfBulkElements
        Element => Mesh % Elements(t)

        ! This is just a quick hack to not consider those element in averaging that don't
        ! even have one face on the active set of nodes. 
        IF( ALLOCATED(NodePerm) ) THEN
          ElemFam = Element % TYPE % ElementCode / 100 
          l = COUNT( NodePerm(Element % NodeIndexes ) > 0 )
          SELECT CASE(ElemFam)          
          CASE(3,4)
            IF(l<2) CYCLE
          CASE(5,6,7)
            IF(l<3) CYCLE
          CASE(8)
            IF(l<4) CYCLE          
          END SELECT
        END IF

        DO i=1,Element % TYPE % NumberOfNodes
          j = Element % DGIndexes(i)
          k = Element % NodeIndexes(i)  

          NodalCnt(k) = NodalCnt(k) + 1
          NodalPerm(k) = k

          j = Solution % Perm(j)
          IF(j==0) CYCLE

          DO l=1,dofs
            NodalVals(dofs*(k-1)+l) = NodalVals(dofs*(k-1)+l) + Solution % Values(dofs*(j-1)+l)
          END DO
        END DO
      END DO

      DO i=1,dofs
        WHERE ( NodalCnt > 0 )
          NodalVals(i::dofs) = NodalVals(i::dofs) / NodalCnt 
        END WHERE
      END DO

      TmpVals => Solution % Values
      TmpPerm => Solution % Perm

      Solution % Perm => NodalPerm
      Solution % Values => NodalVals


    END SUBROUTINE CreateTemporalNodalField




    SUBROUTINE WriteGmshData()
      INTEGER :: ii
      LOGICAL :: DgVar


      ! Time is needed
      !-------------------------------------------------
      TimeVariable => VariableGet( Model % Variables, 'Time' )        
      Time = TimeVariable % Values(1)

      ! Loop over different type of variables
      !-------------------------------------------------
      DO Rank = 0,2
        DO Vari = 1, 999
          IF(Rank==0) WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
          IF(Rank==1) WRITE(Txt,'(A,I0)') 'Vector Field ',Vari
          IF(Rank==2) WRITE(Txt,'(A,I0)') 'Tensor Field ',Vari

          FieldName = ListGetString( Params, TRIM(Txt), Found )
          IF(.NOT. Found) EXIT 
          IF( Rank == 2) THEN
            CALL Warn(Caller,'Not implemented yet for tensors!')
            CYCLE
          END IF

          ComponentVector = .FALSE.
          Solution => VariableGet( Mesh % Variables, FieldName )
          DGVar = .FALSE.

          IF(ASSOCIATED(Solution)) THEN
            DGVar = ( Solution % TYPE == Variable_on_nodes_on_elements )          
            IF(DgVar) CALL CreateTemporalNodalField(Mesh,Solution)

            Values => Solution % Values
            Perm => Solution % Perm
            dofs = Solution % DOFs
          ELSE
            IF( Rank == 1 ) THEN
              Solution => VariableGet( Mesh % Variables, FieldName//' 1' )
              IF( ASSOCIATED( Solution ) ) THEN
                ComponentVector = .TRUE.
                Values => Solution % Values
                Perm => Solution % Perm
                dofs = 1
                Solution => VariableGet( Mesh % Variables, FieldName//' 2' )
                IF( ASSOCIATED(Solution)) THEN
                  Values2 => Solution % Values
                  dofs = 2
                END IF
                Solution => VariableGet( Mesh % Variables, FieldName//' 3' )
                IF( ASSOCIATED(Solution)) THEN
                  Values3 => Solution % Values
                  dofs = 3
                END IF
              END IF
            END IF
            IF( .NOT. ASSOCIATED(Solution)) THEN
              CALL Warn(Caller,'Variable not present: '//TRIM(FieldName))
              CYCLE
            END IF
          END IF

          CALL Info(Caller,'Saving nodal variable: '//TRIM(FieldName),Level=12)
                   
          IF( ASSOCIATED(Solution % EigenVectors) ) THEN
            CALL Warn(Caller,'Eigenvectors related to field: '//TRIM(FieldName))
            CALL Warn(Caller,'Eigenvectors saving yet not supported')
          END IF

          truedim = MIN(dofs, dim)
          nsize = NumberOfGeomNodes

          WRITE(GmshUnit,'(A)') '$NodeData'
          WRITE(GmshUnit,'(A)') '1'
          WRITE(GmshUnit,'(A)') '"'//TRIM(FieldName)//'"'
          WRITE(GmshUnit,'(A)') '1'

          ! Gmsh starts steady state indexes from zero, hence deductions by one
          IF( Transient ) THEN
            WRITE(GmshUnit,'(ES16.7E3)') Time
          ELSE
            WRITE(GmshUnit,'(ES16.7E3)') Time - 1.0_dp
          END IF
          WRITE(GmshUnit,'(A)') '3'
          WRITE(GmshUnit,'(I8)') VisitedTimes-1
          IF(Rank == 0) THEN
            WRITE(GmshUnit,'(A)') '1'
          ELSE IF(Rank == 1) THEN
            WRITE(GmshUnit,'(A)') '3'
          ELSE 
            WRITE(GmshUnit,'(A)') '9'
          END IF
          WRITE(GmshUnit,'(I8)') nsize

          DO ii = 1, NumberOfGeomNodes
            IF( NoPermutation ) THEN
              i = ii 
            ELSE
              i = InvFieldPerm(ii) 
            END IF

            IF( ASSOCIATED( Perm ) ) THEN
              j = Perm(i)
            ELSE
              j = i
            END IF

            IF( Rank == 0 ) THEN
              WRITE(GmshUnit,'(I8,ES16.7E3)') ii,Values(j)
            ELSE IF(Rank == 1) THEN
              IF( j == 0 ) THEN
                WRITE(GmshUnit,'(I8,A)') ii,' 0.0 0.0 0.0'                
              ELSE IF( ComponentVector ) THEN
                IF( truedim == 2 ) THEN
                  WRITE(GmshUnit,'(I8,2ES16.7E3,A)') ii,&
                      Values(j),Values2(j),' 0.0'
                ELSE
                  WRITE(GmshUnit,'(I8,3ES16.7E3)') ii,&
                      Values(j),Values2(j),Values3(j)
                END IF
              ELSE
                IF( truedim == 2 ) THEN
                  WRITE(GmshUnit,'(I8,2ES16.7E3,A)') ii,&
                      Values(dofs*(j-1)+1),Values(dofs*(j-1)+2),' 0.0'
                ELSE
                  WRITE(GmshUnit,'(I8,3ES16.7E3)') ii,&
                      Values(dofs*(j-1)+1),Values(dofs*(j-1)+2),Values(dofs*(j-1)+3)
                END IF
              END IF
            END IF
          END DO
          WRITE(GmshUnit,'(A)') '$EndNodeData'

          IF(DgVar) CALL CreateTemporalNodalField(Mesh,Solution,Revert=.TRUE.)

        END DO
      END DO
    END SUBROUTINE WriteGmshData

!------------------------------------------------------------------------------
  END SUBROUTINE SaveGmshOutput
!------------------------------------------------------------------------------

END MODULE SaveUtils
  
