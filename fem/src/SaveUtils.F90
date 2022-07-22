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
  USE Lists
  USE Messages
  USE MeshUtils, ONLY: GetLagrangeIndexes
  
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

    !PRINT *,'Code:',ElmerCode, VtkCode   

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
          IF (ASSOCIATED(Parent % DGIndexes) ) THEN
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
        PRINT *,'Problematic BC elem:',Element % BodyId, Element % ElementIndex, Element % NodeIndexes, &
            ASSOCIATED( Element % DgIndexes ), ASSOCIATED( Element % BoundaryInfo ), DGelem, &
            Element % TYPE % ElementCode
        CALL Fatal('Elmer2VtkIndexes','Could not set indexes for boundary element!')        
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
    CHARACTER(MAX_NAME_LEN) :: Str, MaskName
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
      CALL Info(Caller,'Allocating NodePerm of size: '//TRIM(I2S(n)),Level=15)
      ALLOCATE(NodePerm(n))
    END IF
    NodePerm = 0
    
    IF(.NOT. ALLOCATED(ActiveElem) ) THEN
      n = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      CALL Info(Caller,'Allocating ActiveElem of size: '//TRIM(I2S(n)),Level=15)
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
      CALL Info(Caller,'Mask is positive for nodes: '//TRIM(I2S(NumberOfGeomNodes)),Level=15)
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
      CALL Info(Caller,'Number of active elements '//TRIM(I2S(NumberOfElements))//&
          ' out of '//TRIM(I2S(Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements)),Level=10)      
      CALL Info(Caller,'Number of geometry nodes '//TRIM(I2S(NumberOfGeomNodes))//&
          ' out of '//TRIM(I2S(Mesh % NumberOfNodes)),Level=10)
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
        CALL Info(Caller,'Maximum number of dofs in DG: '//TRIM(I2S(k)),Level=12)
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
          CALL Info(Caller,'Independent dofs in discontinuous mesh: '//TRIM(I2S(l)),Level=10)
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
      CALL Info(Caller,'Creating permutation for order '//TRIM(I2S(LagN))//' Lagrange nodes!', Level=12)

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
      
      CALL Info(Caller,'Number of dofs for higher order Lagrange elements: '//TRIM(I2S(m)),Level=12)
    ELSE
      NoPermutation = ( NumberOfGeomNodes == Mesh % NumberOfNodes )    
      IF( NoPermutation ) THEN
        DEALLOCATE( NodePerm ) 
      ELSE
        CALL Info(Caller,'Not saving all nodes, creating permutation!',Level=12)
        IF( ALLOCATED( InvNodePerm ) ) DEALLOCATE( InvNodePerm ) 
        
        ALLOCATE( InvNodePerm( NumberOfGeomNodes ) ) 
        CALL Info(Caller,'Allocating InvNodePerm of size: '//TRIM(I2S(NumberOfGeomNodes)),Level=15)
        
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


  
END MODULE SaveUtils
  
