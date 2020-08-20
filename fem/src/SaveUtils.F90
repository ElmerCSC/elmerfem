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
  IMPLICIT NONE

CONTAINS

  ! Given different criteria fos saving create a geometrical mask for elements
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
    
    IF(.NOT. ALLOCATED( NodePerm ) ) ALLOCATE(NodePerm(Mesh % NumberOfNodes))
    NodePerm = 0
    
    IF(.NOT. ALLOCATED(ActiveElem) ) &
        ALLOCATE(ActiveElem(Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements))
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
      END IF
    END IF

    SaveBoundariesOnly = ListGetLogical( Params,'Save Boundaries Only',GotIt ) 
    IF( SaveBoundariesOnly ) CALL Info(Caller,'Saving only boundary elements!',Level=8)
    
    SaveBulkOnly = ListGetLogical( Params,'Save Bulk Only',GotIt ) 
    IF( SaveBulkOnly ) CALL Info(Caller,'Saving only bulk elements!',Level=8)
    
    NumberOfGeomNodes = Mesh % NumberOfNodes
    IF( MaskExists ) THEN
      NumberOfGeomNodes = COUNT( MaskPerm(1:NumberOfGeomNodes) > 0 ) 
    END IF
    NumberOfElements = 0

    IF( NumberOfGeomNodes > 0 ) THEN
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
        IF (.NOT. IsBoundaryElement .AND. Element % BodyId < 1) CYCLE

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
          IF( i <= Mesh % NumberOfBulkElements ) THEN
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
        END IF

        IF( GotMaskCond ) THEN
          n = Element % TYPE % NumberOfNodes
          Indexes => Element % NodeIndexes

          IF( i <= Mesh % NumberOfBulkElements ) THEN
            l = Element % BodyId
            k = ListGetInteger( Model % Bodies(l) % Values,'Body Force',GotIt)
            IF( GotIt ) THEN
              MaskCond(1:n) = ListGetReal( Model % BodyForces(k) % Values, TRIM(MaskName), &
                  n, Indexes, GotIt)
            END  IF

            IF( .NOT. Hit ) THEN
              k = ListGetInteger( Model % Bodies(l) % Values,'Equation',GotIt)
              IF( GotIt ) THEN
                MaskCond(1:n) = ListGetReal( Model % Equations(k) % Values, TRIM(MaskName), &
                    n, Indexes, GotIt)
              END IF
            END IF
          ELSE
            GotIt = .FALSE.
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
          IF( .NOT. ALL(MaskCond(1:n) > 0.0_dp ) ) CYCLE
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
          NodePerm( Element % NodeIndexes ) = 1
        END IF

      END DO

      CALL Info(Caller,'Number of active elements '//TRIM(I2S(NumberOfElements))//&
          ' out of '//TRIM(I2S(Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements)),Level=7)

      NumberOfGeomNodes = COUNT( NodePerm > 0 ) 

      CALL Info(Caller,'Number of geometry nodes '//TRIM(I2S(NumberOfGeomNodes))//&
          ' out of '//TRIM(I2S(Mesh % NumberOfNodes)),Level=7)
    END IF

  END SUBROUTINE GenerateSaveMask
    

  ! Given the geometric permutation, create the dof permutation used in saving
  ! the different parts.
  !-----------------------------------------------------------------------------  
  SUBROUTINE GenerateSavePermutation(Mesh,DG,DN,SaveLinear,ActiveElem,NumberOfGeomNodes,&
      NoPermutation,NumberOfDofNodes,DgPerm,InvDgPerm,NodePerm,InvNodePerm)
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL :: DG, DN, SaveLinear
    LOGICAL, ALLOCATABLE :: ActiveElem(:)
    INTEGER :: NumberOfGeomNodes,NumberOfDofNodes
    LOGICAL :: NoPermutation
    INTEGER, ALLOCATABLE :: DgPerm(:),InvDgPerm(:),NodePerm(:),InvNodePerm(:)

    INTEGER, ALLOCATABLE :: BodyVisited(:)
    INTEGER :: i,j,k,l,m,n
    INTEGER :: Sweep
    INTEGER, POINTER :: NodeIndexes(:) 
    TYPE(Element_t), POINTER :: Element 
    TYPE(Model_t), POINTER :: Model
    CHARACTER(*), PARAMETER :: Caller = 'GenerateSavePermutation'

    
    Model => CurrentModel
        
    NumberOfDofNodes = 0
    IF( DG .OR. DN ) THEN
      NoPermutation = .FALSE.

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

    ELSE
      NoPermutation = ( NumberOfGeomNodes == Mesh % NumberOfNodes )    
      IF( NoPermutation ) THEN
        DEALLOCATE( NodePerm ) 
      ELSE
        CALL Info(Caller,'Not saving all nodes, creating permutation!',Level=12)
        IF( ALLOCATED( InvNodePerm ) ) DEALLOCATE( InvNodePerm ) 
        ALLOCATE( InvNodePerm( NumberOfGeomNodes ) ) 
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
  
