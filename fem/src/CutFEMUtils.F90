!/*****************************************************************************/
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
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 01 Oct 1996
! *
! ****************************************************************************/

!> \ingroup ElmerLib 
!> \{

!------------------------------------------------------------------------------
!>  Module containing utilities for CutFEM style of strategies.
!------------------------------------------------------------------------------

MODULE CutFemUtils
  USE Types
  USE Lists
  USE ElementUtils, ONLY : FreeMatrix
  USE Interpolation, ONLY : CopyElementNodesFromMesh
  USE ElementDescription
  USE MatrixAssembly
  USE MeshUtils, ONLY : AllocateMesh, FindMeshEdges, MeshStabParams
  USE ModelDescription, ONLY : FreeMesh
  USE SolverUtils, ONLY : GaussPointsAdapt, SolveLinearSystem, VectorValuesRange
  USE ParallelUtils
  USE MeshUtils, ONLY : PointInMesh
  
  IMPLICIT NONE

  PRIVATE
  
  LOGICAL :: CutExtend, CutExtrapolate
  LOGICAL, ALLOCATABLE :: CutDof(:)
  INTEGER, POINTER :: ExtendPerm(:) => NULL(), OrigMeshPerm(:) => NULL(), &
      CutPerm(:) => NULL(), PhiPerm(:) => NULL()
  REAL(KIND=dp), POINTER :: OrigMeshValues(:) => NULL(), CutValues(:) => NULL(), &
      ExtendValues(:) => NULL(), PhiValues(:) => NULL(), &
      OrigPrevMeshValues(:,:) => NULL(), PrevCutValues(:,:) => NULL()
  INTEGER, POINTER :: OrigActiveElements(:), AddActiveElements(:), UnsplitActiveElements(:)
  REAL(KIND=dp), ALLOCATABLE, TARGET :: CutInterp(:)
  TYPE(Matrix_t), POINTER :: NodeMatrix
  INTEGER :: CutFemBody
  CHARACTER(:), ALLOCATABLE :: CutStr
  INTEGER :: CutDofs = 0
  INTEGER :: nCase(20)

#define DEBUG_ORIENT 0
#if DEBUG_ORIENT
  REAL(KIND=dp) :: CutFEMCenter(3)
#endif

    
  PUBLIC :: CreateCutFEMMatrix, CreateCutFEMMesh, CreateCutFEMPerm, CreateCutFEMAddMesh, &
      CutFEMVariableFinalize, CutFEMSetOrigMesh, CutFEMSetAddMesh, LevelSetUpdate, &
      CutInterfaceBC, CutInterfaceBulk, CutInterfaceCheck

  PUBLIC :: CutInterp
  
  TYPE(Mesh_t), POINTER :: CutFEMOrigMesh => NULL(), CutFEMAddMesh => NULL()
  
    
CONTAINS


  ! Given a levelset function create a permutation that tells which
  ! edges and which nodes are being cut by the zero levelset.
  ! Optionally also create a permutation for the outside mesh. 
  !------------------------------------------------------------------
  SUBROUTINE CreateCutFEMPerm(Solver,UpdateCoords)
    TYPE(Solver_t) :: Solver
    LOGICAL :: UpdateCoords

    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(ValueList_t), POINTER :: Params

    INTEGER :: i,j,k,nn,ne,body_in,body_out,body_cut,InsideCnt(3),dofs
    REAL(KIND=dp) :: h1,h2,hprod,Eps,r,MaxRat
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Variable_t), POINTER :: Var, PhiVar
    TYPE(Element_t), POINTER :: Element, pElement
    CHARACTER(:), ALLOCATABLE :: str       
    LOGICAL :: Found, PassiveInside, PassiveOutside, isCut, isMore, UseAbsEps, Hit
    REAL(KIND=dp), POINTER :: xtmp(:)    
    LOGICAL :: UpdateOrigCoords
    CHARACTER(*), PARAMETER :: Caller = 'CreateCutFEMPerm'


    Params => Solver % Values    
    Mesh => Solver % Mesh

    ! Memorize original nodal variable and matrix.
    OrigMeshValues => NULL()
    OrigPrevMeshValues => NULL()
    OrigMeshPerm => NULL()
    IF(ASSOCIATED(Solver % Variable ) ) THEN
      IF(ASSOCIATED(Solver % Variable % Perm ) ) THEN
        OrigMeshValues => Solver % Variable % Values
        OrigMeshPerm => Solver % Variable % Perm
        OrigPrevMeshValues => Solver % Variable % PrevValues
      END IF
      CutDofs = Solver % Variable % dofs
      dofs = CutDofs
    END IF
    
    CutFEMOrigMesh => Solver % Mesh 
    OrigActiveElements => Solver % ActiveElements

    NodeMatrix => Solver % Matrix

    CutExtend = ListGetLogical( Params,'CutFEM extend',Found )
    CutExtrapolate = ListGetLogical( Params,'CutFEM extrapolate',Found )
    UpdateOrigCoords = ListGetLogical( Params,'CutFEM bodyfitted',Found )

    ! We always need mesh edges since the new dofs are created in intersections of levelset and edge. 
    IF(ASSOCIATED(Mesh % edges)) THEN
      CALL Info(Caller,'Mesh edges already created!',Level=12)
    ELSE
      CALL Info(Caller,'Create element edges',Level=10)
      CALL FindMeshEdges( Mesh )
      
      ! We need global numbering for the edges that we use for the unique numbering of new nodes
      IF( ParEnv % PEs > 1 ) THEN
        CALL Info(Caller,'Numbering Mesh edges in parallel')
        CALL SParEdgeNumbering(Mesh)
      END IF
    END IF

    nn = Mesh % NumberOfNodes
    ne = Mesh % NumberOfEdges

    IF( UpdateCoords ) THEN
      i = SIZE(Mesh % Nodes % x)
      IF(i < nn + ne ) THEN
        CALL Info(Caller,'Enlarging node coordinates for edge cuts from '&
            //I2S(i)//' to '//I2S(nn+ne),Level=7)
        ALLOCATE(xtmp(nn))
        xtmp = Mesh % Nodes % x(1:nn)
        DEALLOCATE(Mesh % Nodes % x)
        ALLOCATE(Mesh % Nodes % x(nn+ne))
        Mesh % Nodes % x(1:nn) = xtmp

        xtmp = Mesh % Nodes % y(1:nn)
        DEALLOCATE(Mesh % Nodes % y)
        ALLOCATE(Mesh % Nodes % y(nn+ne))
        Mesh % Nodes % y(1:nn) = xtmp

        xtmp = Mesh % Nodes % z(1:nn)
        DEALLOCATE(Mesh % Nodes % z)
        ALLOCATE(Mesh % Nodes % z(nn+ne))
        Mesh % Nodes % z(1:nn) = xtmp
        DEALLOCATE(xtmp)
      END IF
    END IF

    IF(.NOT. ALLOCATED(CutDof) ) THEN
      CALL Info(Caller,'Allocating "CutDof" field to indicate levelset cuts!',Level=20)
      ALLOCATE( CutDof(nn+ne) ) 
    END IF      
    CutDof = .FALSE.  

    ! We store the cut for future interpolation. 
    IF(.NOT. ALLOCATED(CutInterp)) THEN
      CALL Info(Caller,'Allocating "CutInterp" for edge related interpolation!',Level=20)
      ALLOCATE(CutInterp(ne))
    END IF
    CutInterp = 0.0_dp
          
    CutStr = ListGetString( Params,'Levelset Variable', Found)
    IF( .NOT. Found ) CutStr = "surface"

    PhiVar => VariableGet(Mesh % Variables, CutStr, ThisOnly=.TRUE.)
    IF(.NOT. ASSOCIATED(PhiVar) ) THEN
      CALL Fatal(Caller,'"Levelset Variable" not available: '//TRIM(CutStr))
    END IF
    PhiValues => PhiVar % Values
    PhiPerm => PhiVar % Perm

    body_in = ListGetInteger( Params,'CutFEm Inside Body',Found )
    IF(.NOT. Found) body_in = CurrentModel % NumberOfBodies
    body_out = ListGetInteger( Params,'CutFem Outside Body',Found )
    IF(.NOT. Found) body_out = body_in+1
    body_cut = MAX(body_in,body_out)

    ! This is a little dirty, we set the interface elements so we recognize them.
    IF(CutExtend) body_cut = body_cut + 1

    Eps = ListGetCReal(Params,'CutFem Epsilon',Found )
    IF(.NOT. Found) Eps = 1.0e-3
    UseAbsEps = ListGetLogical(Params,'CutFEM Epsilon Absolute',Found ) 


    ! First mark the cut nodes.
    ! These could maybe be part of the same loop as well but I separated when testing something.
    DO i=1, Mesh % NumberOfEdges
      NodeIndexes => Mesh % Edges(i) % NodeIndexes
      IF(ANY(PhiPerm(NodeIndexes) == 0)) CYCLE
      h1 = PhiValues(PhiPerm(NodeIndexes(1)))
      h2 = PhiValues(PhiPerm(NodeIndexes(2)))
      hprod = h1*h2            
      IF( hprod < 0.0_dp ) THEN
        r = ABS(h2)/(ABS(h1)+ABS(h2))        
        Hit = .FALSE.
        IF( UseAbsEps ) THEN
          IF(ABS(h2) < Eps ) THEN
            CutDof(NodeIndexes(2)) = .TRUE.
            Hit = .TRUE.
          END IF
          IF(ABS(h1) < Eps ) THEN
            CutDof(NodeIndexes(1)) = .TRUE.
            Hit = .TRUE.
          END IF
        ELSE
          IF( r <= Eps ) THEN
            CutDof(NodeIndexes(2)) = .TRUE.
            Hit = .TRUE.
          END IF
          IF((1.0-r < Eps) ) THEN
            CutDof(NodeIndexes(1)) = .TRUE.
            Hit = .TRUE.
          END IF
        END IF
      ELSE IF( ABS(hprod) < 1.0d-20 ) THEN
        IF(ABS(h1) < 1.0e-20) CutDof(NodeIndexes(1)) = .TRUE. 
        IF(ABS(h2) < 1.0e-20) CutDof(NodeIndexes(2)) = .TRUE.
      END IF
    END DO

    
    IF(ParEnv % PEs > 1 ) THEN
      BLOCK 
        INTEGER, POINTER :: Perm(:)
        INTEGER :: ni
        REAL(KIND=dp), POINTER :: CutDofR(:)

        ni = COUNT( CutDof(1:nn) .AND. Mesh % ParallelInfo % GInterface(1:nn) )
        ni = ParallelReduction( ni ) 

        IF( ni > 0 ) THEN
          ALLOCATE(CutDofR(nn),Perm(nn))
          CutDofR = 0.0_dp
          DO i=1,nn
            Perm(i) = i
          END DO

          WHERE( CutDof(1:nn) )
            CutDofR = 1.0_dp
          END WHERE
          CALL ExchangeNodalVec( Mesh % ParallelInfo, Perm, CutDofR, op = 2)
          DO i=1,nn
            IF(CutDofR(i) > 0.5_dp ) CutDof(i) = .TRUE.
          END DO
          DEALLOCATE(CutDofR, Perm )
        END IF
      END BLOCK
    END IF
    
    
    ! Then mark the edges trying to avoid nearby cuts.  
    InsideCnt = 0
    j = 0

    ! This is an add'hoc value that represents the maximum aspect ratio of elements in the mesh.
    MaxRat = 2.0
    
    DO i=1, Mesh % NumberOfEdges
      NodeIndexes => Mesh % Edges(i) % NodeIndexes
      IF(ANY(PhiPerm(NodeIndexes)==0)) CYCLE
      h1 = PhiValues(PhiPerm(NodeIndexes(1)))
      h2 = PhiValues(PhiPerm(NodeIndexes(2)))
      hprod = h1*h2            
      IF( hprod < 0.0_dp ) THEN
        r = ABS(h2)/(ABS(h1)+ABS(h2))        
        Hit = .FALSE.

        ! We may have a sloppier rule if the dof is already cut?
        ! If the rule is exactly the same then no need for separate loop.
        IF( r <= MaxRat * Eps ) THEN
          IF(CutDof(NodeIndexes(2))) CYCLE
        ELSE IF((1.0-r < MaxRat * Eps) ) THEN
          IF(CutDof(NodeIndexes(1))) CYCLE
        END IF

        j = j+1 
        CutDof(nn+i) = .TRUE.

        ! The interpolation weight should always be [0,1]
        IF(r < 0.0 .OR. r > 1.0) THEN
          PRINT *,'Invalid cutinterp:',i,j,r
        END IF    
        
        CutInterp(i) = r

        ! We update nodes so that the element on-the-fly can point to then using NodeIndexes. 
        IF( UpdateCoords ) THEN
          Mesh % Nodes % x(nn+i) = (1-r) * Mesh % Nodes % x(NodeIndexes(2)) + &
              r * Mesh % Nodes % x(NodeIndexes(1))
          Mesh % Nodes % y(nn+i) = (1-r) * Mesh % Nodes % y(NodeIndexes(2)) + &
              r * Mesh % Nodes % y(NodeIndexes(1))
          Mesh % Nodes % z(nn+i) = (1-r) * Mesh % Nodes % z(NodeIndexes(2)) + &
              r * Mesh % Nodes % z(NodeIndexes(1))
        END IF
      END IF
    END DO

    ! Should we update the original coords for nodes which closely match the levelset but not exactly.
    ! This would be the case if we want to follow the body fitted shape of a object as closely as possible. 
    ! We would not want to do in transient cases 
    IF(UpdateOrigCoords) THEN
      BLOCK 
        LOGICAL, ALLOCATABLE :: MovedNode(:)
        REAL(KIND=dp), ALLOCATABLE :: TmpCoords(:,:)
        ALLOCATE(MovedNode(nn), TmpCoords(nn,3))  
        MovedNode = .FALSE.
        TmpCoords = 0.0_dp
        DO i=1, Mesh % NumberOfEdges
          NodeIndexes => Mesh % Edges(i) % NodeIndexes
          IF(.NOT. ANY(CutDOF(NodeIndexes))) CYCLE

          h1 = PhiValues(PhiPerm(NodeIndexes(1)))
          h2 = PhiValues(PhiPerm(NodeIndexes(2)))
          hprod = h1*h2                    
          IF( hprod >= 0.0_dp ) CYCLE

          r = ABS(h2)/(ABS(h1)+ABS(h2))                        
          IF( r <= Eps ) THEN
            j = 2
          ELSE IF((1.0-r < Eps) ) THEN
            j = 1
          ELSE
            CYCLE
          END IF

          k = NodeIndexes(j)
          IF(.NOT. CutDof(k)) CYCLE
          IF(MovedNode(k)) CYCLE

          TmpCoords(k,1) = (1-r) * Mesh % Nodes % x(NodeIndexes(2)) + &
              r * Mesh % Nodes % x(NodeIndexes(1))
          TmpCoords(k,2) = (1-r) * Mesh % Nodes % y(NodeIndexes(2)) + &
              r * Mesh % Nodes % y(NodeIndexes(1))
          TmpCoords(k,3) = (1-r) * Mesh % Nodes % z(NodeIndexes(2)) + &
              r * Mesh % Nodes % z(NodeIndexes(1))
          MovedNode(k) = .TRUE.
        END DO
        k = COUNT(MovedNode)
        CALL Info(Caller,'Moved cut nodes to be exactly at zero levelset!')

        WHERE(MovedNode(1:nn))
          Mesh % Nodes % x(1:nn) = TmpCoords(1:nn,1)
          Mesh % Nodes % y(1:nn) = TmpCoords(1:nn,2)
          Mesh % Nodes % z(1:nn) = TmpCoords(1:nn,3)        
        END WHERE
        DEALLOCATE(MovedNode, TmpCoords)
      END BLOCK
    END IF

    
    IF(InfoActive(25)) THEN
      PRINT *,'CutInterp interval:',MINVAL(CutInterp),MAXVAL(CutInterp)    
      PRINT *,'Nodes split',COUNT(CutDof(1:nn))
      PRINT *,'Edges cut',COUNT(CutDof(nn+1:nn+ne))
    END IF
      
    ! Set the material for inside/outside or interface material.
    ! This way we do not need to have too complicated material sections.
    CutFEMBody = 0
    DO i=1,Mesh % NumberOfBulkElements
      Element => Mesh % Elements(i)

      NodeIndexes => Element % NodeIndexes
      IF(ANY(PhiPerm(NodeIndexes) == 0)) CYCLE

      ! So far we assume that there is only one body index used to define the CutFEM region.
      ! We are tampering with the index, so we need to store it. 
      IF(CutFEMBody == 0) THEN
        CutFEMBody = Element % BodyId
      ELSE
        IF(CutFemBody /= Element % BodyId ) THEN
          CALL Fatal(Caller,'Modify code to deal with several bodies!')
        END IF
      END IF
      
      j = -1
      IF(ANY(CutDof(nn + Element % EdgeIndexes)) ) THEN
        ! Some edge is split => interface element
        j = body_cut        
        InsideCnt(3) = InsideCnt(3)+1
      ELSE
        ! Also at interface element if we have diagonal split in a quad. 
        IF(Element % TYPE % ElementCode / 100 == 4 ) THEN
          IF(ALL(CutDof(NodeIndexes([1,3])))) THEN
            j = body_cut            
          ELSE IF(ALL(CutDof(NodeIndexes([2,4])))) THEN
            j = body_cut
          END IF
        END IF

        ! Ok, no interface. Use the min/max value to indicate whether this is inside/outside. 
        IF(j<0) THEN
          h1 = MAXVAL(PhiValues(PhiPerm(NodeIndexes)))
          h2 = MINVAL(PhiValues(PhiPerm(NodeIndexes)))
          IF(h1 > -h2) THEN
            InsideCnt(2) = InsideCnt(2)+1
            j = body_out
          ELSE
            InsideCnt(1) = InsideCnt(1)+1
            j = body_in
          END IF
        ELSE
          InsideCnt(3) = InsideCnt(3)+1
        END IF
      END IF

      Element % BodyId = j
    END DO

    IF(InfoActive(25)) THEN
      PRINT *,'Inside/Outside count:',InsideCnt
    END IF
          
    ! CutPerm is the reordered dofs for the CutFEM mesh. 
    IF(.NOT. ASSOCIATED(CutPerm)) THEN
      ALLOCATE(CutPerm(nn+ne))
      CALL info(Caller,'Allocated CutPerm of size: '//I2S(nn+ne),Level=20)
    END IF
    CutPerm = 0

    PassiveOutside = ListGetLogical( Params,'CutFEM Passive Outside',Found ) 
    IF(.NOT. Found ) PassiveOutside = (body_out == 0)
    PassiveInside = ListGetLogical( Params,'CutFEM Passive Inside',Found ) 
    IF(.NOT. Found) PassiveInside = (body_in == 0)

    ! Set all cut dofs to exist.
    WHERE(CutDof)
      CutPerm = 1
    END WHERE
    IF( PassiveOutside ) THEN
      DO i=1,nn
        j = PhiPerm(i)
        IF(j==0) CYCLE
        IF(PhiValues(j) < 0) CutPerm(i) = 1
      END DO
    ELSE IF( PassiveInside ) THEN
      DO i=1,nn
        j = PhiPerm(i)
        IF(j==0) CYCLE
        IF(PhiValues(j) > 0) CutPerm(i) = 1
      END DO
    ELSE
      ! We both inside and outside exist. 
      CutPerm(1:nn) = 1
    END IF

    j = 0
    DO i=1,nn+ne
      IF(CutPerm(i) > 0) THEN
        j = j+1
        CutPerm(i) = j
      END IF
    END DO
    k = COUNT(CutPerm(1:nn)>0) 
    CALL Info(Caller,'CutFEM number of nodes: '//I2S(j)//' (original '//I2S(k)//')',Level=7)


    ! If there is a primary variable associated to the original mesh copy that to the new mesh.
    IF(ASSOCIATED(OrigMeshValues)) THEN
      IF(ASSOCIATED(CutValues)) DEALLOCATE(CutValues)
      ALLOCATE(CutValues(dofs*j))
      CutValues = 0.0_dp

      DO i=1,dofs
        WHERE(CutPerm(1:nn) > 0 )        
          CutValues(dofs*(CutPerm-1)+i) = OrigMeshValues(dofs*(OrigMeshPerm-1)+i) 
        END WHERE
      END DO
        
      ! Point the permutation and values to the newly allocated vectors.
      ! This way 
      Solver % Variable % Perm => CutPerm
      Solver % Variable % Values => CutValues
      
      ! For transient problems do the same for PrevValues
      IF(ASSOCIATED(OrigPrevMeshValues)) THEN
        IF(ASSOCIATED(PrevCutValues)) DEALLOCATE(PrevCutValues)      
        i = SIZE(OrigPrevMeshValues,2)
        ALLOCATE(PrevCutValues(dofs*j,i))
        PrevCutValues = 0.0_dp

        ! Copy nodal values as initial guess to cut fem values. 
#if 0
        ! fix this
        DO l=1,dofs
          DO i=1,nn
            j = CutPerm(i)
            k = OrigMeshPerm(i)
            IF(j==0 .OR. k==0) CYCLE
            OrigMeshValues(dofs*(k-1)+l) = CutValues(dofs*(j-1)+l)
          END DO
        END DO
#endif

        DO i=1,SIZE(OrigPrevMeshValues,2)
          DO j=1,dofs
            WHERE(CutPerm(1:nn) > 0 )        
              PrevCutValues(dofs*(CutPerm(1:nn)-1)+j,i) = &
                  OrigPrevMeshValues(dofs*(OrigMeshPerm(1:nn)-1)+j,i) 
            END WHERE
          END DO
        END DO
        Solver % Variable % PrevValues => PrevCutValues
      END IF
    END IF
        
    ! This in an optional routine if we want to extend the field values outside
    ! active domain. The reason might be to provide better initial values for the new territory. 
    IF(CutExtend) THEN
      CALL Info(Caller,'Extending field outside the active domain!',Level=20)
      r = ListGetCReal( Params,'CutFEM extend width',Found )

      IF(.NOT. ASSOCIATED(ExtendPerm)) THEN
        ALLOCATE(ExtendPerm(nn+ne))
      END IF
      ExtendPerm = 0
      
      r = ListGetCReal( Params,'CutFEM extend width',Found )
      
      ! Set the material for inside/outside.
      DO i=1,Mesh % NumberOfBulkElements
        Element => Mesh % Elements(i)                
        IF(ANY(PhiPerm(Element % NodeIndexes) == 0)) CYCLE
        
        IF( Element % BodyId == body_cut ) THEN
          ! Mark dofs to extend on elements which lack CutFEM dofs. 
10        pElement => CutInterfaceBulk(Element,isCut,isMore)        
          IF(ANY(CutPerm(pElement % NodeIndexes) == 0) ) THEN
            ExtendPerm( pElement % NodeIndexes ) = 1          
          END IF
          IF(IsMore) GOTO 10
          ! Ok, revert the dirty flag. 
          Element % BodyId = body_cut-1
        ELSE          
          IF( ALL( CutPerm( Element % NodeIndexes ) == 0) ) THEN
            IF( Found ) THEN
              ! Mark all dofs within a defined width.
              IF(MINVAL(ABS(PhiVar % Values(PhiVar % Perm(Element % NodeIndexes )))) < r ) THEN
                ExtendPerm( Element % NodeIndexes ) = 1
              END IF
            ELSE
              ! Mark all elements in the outside region.
              ExtendPerm( Element % NodeIndexes ) = 1
            END IF
          END IF
        END IF
      END DO
              
      j = 0
      DO i=1,nn+ne
        IF(ExtendPerm(i) == 0) CYCLE
        j = j+1
        ExtendPerm(i) = j
      END DO

      k = COUNT(ExtendPerm > 0 .AND. CutPerm > 0 )
      CALL Info(Caller,'Interface dofs '//I2S(j)//' (shared '//I2S(k)//')')      
    END IF


    ! This is a dirty way to halt the progress when levelset goes beyond the planned
    ! outer boundaries.
    BLOCK
      TYPE(ValueList_t), POINTER :: BC    
      INTEGER :: bc_id
      k = Mesh % NumberOfBulkElements
      DO i=1,Mesh % NumberOfBoundaryElements
        Element => Mesh % Elements(k+i)
        NodeIndexes => Element % NodeIndexes      

        DO bc_id=1,CurrentModel % NumberOfBCs
          IF ( Element % BoundaryInfo % Constraint == CurrentModel % BCs(bc_id) % Tag ) EXIT
        END DO
        IF ( bc_id > CurrentModel % NumberOfBCs ) CYCLE     

        BC => CurrentModel % BCs(bc_id) % Values        
        IF(ListGetLogical(BC,'CutFem Forbidden Boundary',Found ) ) THEN
          IF(ANY(CutPerm(nn+Element % EdgeIndexes)>0)) THEN
            CALL Fatal(Caller,'CutFEM extends beyond forbidden boundaries!')
          END IF
        END IF
      END DO
    END BLOCK

#if DEBUG_ORIENT 
    r = HUGE(r)
    DO i=1,Mesh % NumberOfNodes
      j = PhiPerm(i)
      IF(j==0) CYCLE
      IF(PhiValues(j) < r) THEN
        r = PhiValues(j)
        CutFEMCenter(1) = Mesh % Nodes % x(i)
        CutFEMCenter(2) = Mesh % Nodes % y(i)
        CutFEMCenter(3) = Mesh % Nodes % z(i)
      END IF
    END DO
    PRINT *,'CutFEMCenter:',CutFEMCenter
#endif
    
    ! This is just counter for different split cases while developing the code. 
    nCase = 0

    Solver % CutInterp => CutInterp 
    
  END SUBROUTINE CreateCutFEMPerm


  ! Given a permutation, create a matrix. We assume simple nodal elements.
  ! Some extra dofs are created since at the interface we assume that
  ! there can be all possible connections. 
  !-----------------------------------------------------------------------
  FUNCTION CreateCutFemMatrix(Solver,Perm,MimicMat) RESULT ( A ) 
    TYPE(Solver_t) :: Solver
    INTEGER :: Perm(:)
    TYPE(Matrix_t), POINTER :: A
    TYPE(Matrix_t), POINTER, OPTIONAL :: MimicMat

    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: i,j,k,l,t,m,n,dofs,nn,active
    INTEGER, ALLOCATABLE :: BlockInds(:),DofInds(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER, SAVE :: AllocVecs(3)
    CHARACTER(*), PARAMETER :: Caller = 'CreateCutFemMatrix'

    Mesh => Solver % Mesh
    CutDofs = Solver % Variable % Dofs
    dofs = CutDofs
    
    ! Create new matrix
    A => AllocateMatrix()
    A % FORMAT = MATRIX_LIST

    ! Add extreme entry since list matrix likes to be allocated at once. 
    n = dofs * MAXVAL(Perm)
    IF(n==0) THEN
      CALL Warn(Caller,'CutFEM matrix size is zero?')
      A % NumberOfRows = 0
      RETURN
    END IF
      
    CALL Info(Caller,'Size of CutFEM matrix with '//I2S(dofs)//' dofs is: '//I2S(n),Level=10)
    

    CALL List_AddToMatrixElement(A % ListMatrix, n, n, 0.0_dp ) 

    n = 2*Mesh % MaxElementNodes
    ALLOCATE(BlockInds(n),DofInds(n*dofs))
    BlockInds = 0

    nn = Mesh % NumberOfNodes

    DO t=1,Mesh % NumberOfBulkElements
      Element => Mesh % Elements(t)
      IF(ANY(PhiPerm(Element % NodeIndexes) == 0)) CYCLE

      ! Add active node indexes.
      m = 0
      n = Element % TYPE % NumberOfNodes
      DO i=1,n
        j = Perm(Element % NodeIndexes(i))
        IF(j==0) CYCLE
        m = m+1
        BlockInds(m) = j
      END DO

      ! Add active edge indexes after node indexes.  
      n = Element % TYPE % NumberOfEdges
      DO i=1,n
        j = Perm(nn + Element % EdgeIndexes(i))
        IF(j==0) CYCLE
        m = m+1
        BlockInds(m) = j
      END DO

      ! For vector valued problems add the number of dof indeces.
      IF( dofs == 1 ) THEN
        DofInds(1:m) = BlockInds(1:m)
      ELSE
        DO i=0,dofs-1
          DofInds(m*i+1:m*(i+1)) = dofs*(BlockInds(1:m)-1)+(i+1)
        END DO
        m = m*dofs
      END IF

      ! Add locations to matrix. We add zeros since we are only creating the topology, not assembling. 
      DO i=1,m
        DO j=1,m
          CALL List_AddToMatrixElement(A % ListMatrix,DofInds(i),DofInds(j),0.0_dp)
        END DO
      END DO
    END DO

    ! Make a CRS matrix that has now a topology to account for all entries coming from cutfem. 
    CALL List_toCRSMatrix(A)
    CALL CRS_SortMatrix(A,.FALSE.)

    IF(.NOT. ASSOCIATED(A % rhs)) THEN
      ALLOCATE(A % rhs(A % NumberOfRows))
    END IF
    A % rhs = 0.0_dp


    ! MimicMat is the matrix which we should replace. So if it is transient, it will have MassValues etc. 
    IF(PRESENT(MimicMat)) THEN
      ! If the matrix does not exist, do not update.
      IF(ASSOCIATED(MimicMat)) THEN 
        AllocVecs = 0
        IF(ASSOCIATED(MimicMat % MassValues)) AllocVecs(1) = 1
        IF(ASSOCIATED(MimicMat % DampValues)) AllocVecs(2) = 1
        IF(ASSOCIATED(MimicMat % Force)) AllocVecs(3) = SIZE(MimicMat % Force,2)
      END IF
      IF(AllocVecs(1) > 0 ) THEN
        ALLOCATE(A % MassValues(SIZE(A % Values)))
        A % MassValues = 0.0_dp
      END IF
      IF(AllocVecs(2) > 0) THEN
        ALLOCATE(A % DampValues(SIZE(A % Values)))
        A % DampValues = 0.0_dp
      END IF
      n = AllocVecs(3) 
      IF(n > 0 ) THEN
        ALLOCATE(A % Force(A % NumberOfRows,n))
        A % Force = 0.0_dp
      END IF
    END IF
    
  END FUNCTION CreateCutFemMatrix


  ! This is a routine that just checks whether an element is cut.
  !----------------------------------------------------------------
  SUBROUTINE CutInterfaceCheck( Element, IsCut, IsActive, ExtPerm )
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: IsCut, IsActive
    INTEGER, POINTER, OPTIONAL :: ExtPerm(:)

    INTEGER, SAVE :: n_split, n_cut, n_act
    INTEGER :: j,j2,j3,nn,ne
    LOGICAL :: Visited = .FALSE.
    INTEGER, POINTER :: Perm(:)
    TYPE(Mesh_t), POINTER :: Mesh
    CHARACTER(*), PARAMETER :: Caller = 'CutInterfaceCheck'

    SAVE Visited, Mesh, nn

    IF(.NOT. Visited) THEN
      Mesh => CurrentModel % Solver % Mesh
      nn = Mesh % NumberOfNodes
      Visited = .TRUE.
    END IF

    IF( PRESENT(ExtPerm) ) THEN
      Perm => ExtPerm
    ELSE      
      Perm => CurrentModel % Solver % Variable % Perm     
    END IF
      
    n_split = COUNT( CutDof(nn + Element % EdgeIndexes) )
    n_cut = COUNT( CutDof(Element % NodeIndexes) )

    n_act = COUNT( Perm(Element % NodeIndexes) > 0 )

    IsCut = ( n_split > 0 .OR. n_cut > 1 )
    IsActive = (n_act == Element % TYPE % numberOfNodes ) .OR. IsCut

  END SUBROUTINE CutInterfaceCheck
      
  
  ! Given Element, levelset function and the CutDof field return information whether the element
  ! is cut and if it, should we call the routine again for the next split. 
  !----------------------------------------------------------------------------------------------
  FUNCTION CutInterfaceBulk( Element, IsCut, IsMore ) RESULT ( pElement )
    TYPE(Element_t), POINTER :: Element, pElement
    LOGICAL :: IsCut
    LOGICAL :: IsMore

    TYPE(Element_t), TARGET :: Elem303, Elem404, Elem706, Elem808
    TYPE(Element_t), POINTER :: prevElement
    INTEGER :: SgnNode, i, n, nn, ElemType, body_out, body_in, CutCnt
    LOGICAL :: Found
    REAL(KIND=dp), POINTER :: x(:), y(:), z(:)
    CHARACTER(:), ALLOCATABLE :: str       
    TYPE(Variable_t), POINTER :: PhiVar !Var
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Solver_t), POINTER :: Solver => NULL()
    CHARACTER(*), PARAMETER :: Caller = 'CutInterfaceBulk'
    TYPE(Nodes_t) :: ElemNodes
    INTEGER, ALLOCATABLE :: LocalInds(:), ElemInds(:)
    LOGICAL, ALLOCATABLE :: ElemCut(:)
    
        
    SAVE Mesh, Solver, x, y, z, Elem303, Elem404, body_in, body_out, &
        nn, CutCnt, PhiVar, ElemNodes, ElemInds, ElemCut, ElemType, LocalInds, &
        prevElement
    
    IF(.NOT. ASSOCIATED( Solver, CurrentModel % Solver ) ) THEN
      Mesh => CurrentModel % Solver % Mesh
      Solver => CurrentModel % Solver

      IF(.NOT. ASSOCIATED(ElemNodes % x)) THEN
        n = 8
        ALLOCATE( ElemNodes % x(n), ElemNodes % y(n), ElemNodes % z(n), ElemInds(n), &
            ElemCut(n), LocalInds(4))        
      END IF
              
      nn = Mesh % NumberOfNodes
      x => Mesh % Nodes % x
      y => Mesh % Nodes % y
      z => Mesh % Nodes % z

      ! Create empty element skeletons that are filled when splitting elements. 
      Elem303 % TYPE => GetElementType(303)
      ALLOCATE(Elem303 % NodeIndexes(3))      
      Elem303 % NodeIndexes = 0
      
      Elem404 % TYPE => GetElementType(404)
      ALLOCATE(Elem404 % NodeIndexes(4))      
      Elem404 % NodeIndexes = 0

      PhiVar => VariableGet(Mesh % Variables, CutStr, ThisOnly=.TRUE.)
      IF(.NOT. ASSOCIATED(PhiVar) ) THEN
        CALL Fatal(Caller,'"Levelset Variable" not available: '//TRIM(CutStr))
      END IF

      body_in = ListGetInteger( Solver % Values,'CutFEm Inside Body',Found )
      IF(.NOT. Found) body_in = CurrentModel % NumberOfBodies
      body_out = ListGetInteger( Solver % Values,'CutFem Outside Body',Found )
      IF(.NOT. Found) body_out = body_in+1
    END IF

    
    ! This is the counter for splitting.
    IF(.NOT. ASSOCIATED(prevElement,Element)) THEN
      CutCnt = 1
      prevElement => Element
      ElemType = Element % Type % ElementCode
      n = ElemType / 100

      ! For triangles & quads these are true, not for others...
      ElemInds(1:n) = Element % NodeIndexes(1:n)
      ElemInds(n+1:2*n) = nn + Element % EdgeIndexes(1:n)

      ElemCut(1:2*n) = CutDof(ElemInds(1:2*n))

      ElemNodes % x(1:2*n) = x(ElemInds(1:2*n))
      ElemNodes % y(1:2*n) = y(ElemInds(1:2*n))
      ElemNodes % z(1:2*n) = z(ElemInds(1:2*n))      
    ELSE
      CutCnt = CutCnt+1
    END IF

    pElement => Element
    CALL SplitSingleElement(Element, ElemCut, ElemNodes, CutCnt, &
        IsCut, IsMore, LocalInds, SgnNode )
    IF(.NOT. IsCut) RETURN
    
    i = COUNT(LocalInds > 0)
    SELECT CASE(i)
    CASE(3) 
      pElement => Elem303
    CASE(4)
      pElement => Elem404
    CASE DEFAULT
      CALL Fatal('CutInterfaceBulk','Impossible number of nodes!')
    END SELECT
    pElement % NodeIndexes(1:i) = ElemInds(LocalInds(1:i)) 
          
    ! This circumwents some rare case when node is cut.
    IF( body_out == 0 ) THEN
      IF( ALL( CutPerm(pElement % NodeIndexes) > 0) ) THEN
        pElement % BodyId = body_in
      ELSE
        pElement % BodyId = body_out
      END IF
    ELSE
      i = PhiVar % Perm(pElement % NodeIndexes(sgnNode))
      IF( PhiVar % Values(i) > 0.0_dp ) THEN
        pElement % BodyId = body_out      
      ELSE
        pElement % BodyId = body_in
      END IF
    END IF
      
  END FUNCTION CutInterfaceBulk
    


  ! Given Element, levelset function and the CutDof field return the elements created at the interface.
  ! In 2D and also for straight cuts in 3D we should expect to have just one element but it could
  ! be changed. This includes also wedges even if the above does not just to model how they could
  ! work in 3D. 
  !----------------------------------------------------------------------------------------------
  FUNCTION CutInterfaceBC( Element, IsCut, IsMore ) RESULT ( pElement )
    TYPE(Element_t), POINTER :: Element, pElement
    LOGICAL :: IsCut, IsMore

    TYPE(Element_t), TARGET :: Elem202, Elem303, Elem404
    TYPE(Element_t), POINTER, SAVE :: prevElement
    INTEGER :: m, n, n_split, n_cut, i, j, j2, j3, j4, nn, SplitCase
    INTEGER, POINTER :: nIndexes(:), eIndexes(:)
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL :: Visited = .FALSE., Found, VerticalCut
    REAL(KIND=dp), POINTER :: x(:), y(:), z(:)
    TYPE(Solver_t), POINTER :: Solver
    CHARACTER(*), PARAMETER :: Caller = 'CutInterfaceBC'

    SAVE Visited, Mesh, Solver, nn, x, y, z, n_split, n_cut, &
        Elem202, Elem303, Elem404, VerticalCut, m


    IF(.NOT. Visited) THEN      
      Solver => CurrentModel % Solver
      Mesh => Solver % Mesh
      nn = Mesh % NumberOfNodes
      x => Mesh % Nodes % x
      y => Mesh % Nodes % y
      z => Mesh % Nodes % z                 

      n = ListGetInteger( Solver % Values,'CutFem Interface BC',Found )

      IF( Mesh % MeshDim == 3 ) THEN
        Elem303 % TYPE => GetElementType(303)
        ALLOCATE(Elem303 % NodeIndexes(3))      
        Elem303 % NodeIndexes = 0
        ALLOCATE( Elem303 % BoundaryInfo )      
        Elem303 % BoundaryInfo % Constraint = n

        Elem404 % TYPE => GetElementType(404)
        ALLOCATE(Elem404 % NodeIndexes(4))      
        Elem404 % NodeIndexes = 0
        ALLOCATE( Elem404 % BoundaryInfo )      
        Elem404 % BoundaryInfo % Constraint = n

        VerticalCut = ListGetLogical(Solver % Values,'CutFEM vertical cut',Found ) 
      ELSE
        Elem202 % TYPE => GetElementType(202)
        ALLOCATE(Elem202 % NodeIndexes(2))      
        Elem202 % NodeIndexes = 0
        ALLOCATE( Elem202 % BoundaryInfo )      
        Elem202 % BoundaryInfo % Constraint = n

        VerticalCut = .FALSE.
      END IF
      
      Visited = .TRUE.
    END IF

    nIndexes => Element % NodeIndexes
    eIndexes => Element % EdgeIndexes
    
    
    ! This is the counter for splitting.
    IF(.NOT. ASSOCIATED(prevElement,Element)) THEN
      m = 1
      prevElement => Element
      IF( VerticalCut ) THEN
        n = SIZE(eIndexes) / 2
        n_split = COUNT( CutDof(nn + eIndexes(1:n)) )
        n = SIZE(nIndexes) / 2
        n_cut = COUNT( CutDof(nIndexes(1:n)) )
      ELSE                
        n_split = COUNT( CutDof(nn + eIndexes) )
        n_cut = COUNT( CutDof(nIndexes) )
      END IF
    ELSE
      m = m+1
    END IF

    IsMore = .FALSE.

    IF( n_split == 0 .AND. n_cut <= 1 ) THEN
      isCut = .FALSE.
      pElement => NULL()
      RETURN
    END IF

    IsCut = .TRUE.
    
    SELECT CASE( Element % TYPE % ElementCode )
    CASE( 808 ) 
      pElement => Elem303
    CASE( 706 ) 
      pElement => Elem404      
    CASE( 504 )
      pElement => Elem303
    CASE( 303, 404 ) 
      pElement => Elem202      
    CASE DEFAULT
      CALL Fatal(Caller,'Unknown element type to split: '//I2S(Element % TYPE % ElementCode)//'!')
    END SELECT   
    pElement % NodeIndexes = 0
    
    
    ! This allows use case to deal with element types, edge splits and node splits at the same time. 
    SplitCase = 100 * Element % TYPE % ElementCode + 10 * n_split + n_cut
    

    SELECT CASE( SplitCase ) 
      
    CASE( 30320 )    
      ! Find the both cut edges
      DO j=1,3
        IF( CutDof( nn + eIndexes(j) ) ) EXIT
      END DO
      DO j2=1,3
        IF(j2==j) CYCLE
        IF( CutDof( nn + eIndexes(j2) ) ) EXIT
      END DO
      pElement % NodeIndexes(1) = nn + eIndexes(j)
      pElement % NodeIndexes(2) = nn + eIndexes(j2)

    CASE( 30321 )      
      IF(m==1) THEN
        DO j=1,3
          IF( CutDof( nn + eIndexes(j) ) ) EXIT
        END DO
        DO j2=1,3
          IF(j2==j) CYCLE
          IF( CutDof( nn + eIndexes(j2) ) ) EXIT
        END DO
        pElement % NodeIndexes(1) = nn + eIndexes(j)
        pElement % NodeIndexes(2) = nn + eIndexes(j2)
        IsMore = .TRUE.
      ELSE
        DO j=1,3
          IF(CutDof(nIndexes(j))) EXIT
        END DO
        j2 = j
        IF( .NOT. CutDof(nn + eIndexes(j2) ) ) THEN
          j2 = MODULO(j-2,3)+1          
          IF(.NOT. CutDof(nn + eIndexes(j2))) THEN
            CALL Fatal('Caller','Could not imagine this 303 case!')
          END IF
        END IF
        pElement % NodeIndexes(1) = nIndexes(j)
        pElement % NodeIndexes(2) = nn + eIndexes(j2)        
        IsMore = .FALSE.        
      END IF
        
    CASE( 30311 ) 
      ! Find the edge and node that is cut
      DO j=1,3
        IF( CutDof( nn + eIndexes(j) ) ) EXIT
      END DO
      DO j2=1,3
        IF( CutDof( nIndexes(j2) ) ) EXIT
      END DO
      pElement % NodeIndexes(1) = nn + eIndexes(j)
      pElement % NodeIndexes(2) = nIndexes(j2)

    CASE( 30302 )       
      ! Find the two nodes that are cut.
      DO j=1,3
        IF( CutDof( nIndexes(j) ) ) EXIT
      END DO
      DO j2=1,3
        IF(j2==j) CYCLE
        IF( CutDof( nIndexes(j2) ) ) EXIT
      END DO
      pElement % NodeIndexes(1) = nIndexes(j)
      pElement % NodeIndexes(2) = nIndexes(j2)

    CASE( 40420 )      
      DO j=1,4
        IF(CutDof( nn + eIndexes(j))) EXIT
      END DO
      DO j2=1,4
        IF(j2==j) CYCLE
        IF(CutDof( nn + eIndexes(j2))) EXIT
      END DO
      pElement % NodeIndexes(1) = nn + eIndexes(j)
      pElement % NodeIndexes(2) = nn + eIndexes(j2)

    CASE( 40421 )      

      IF(m==1) THEN      
        DO j=1,4
          IF(CutDof( nn + eIndexes(j))) EXIT
        END DO
        DO j2=1,4
          IF(j2==j) CYCLE
          IF(CutDof( nn + eIndexes(j2))) EXIT
        END DO
        pElement % NodeIndexes(1) = nn + eIndexes(j)
        pElement % NodeIndexes(2) = nn + eIndexes(j2)
        IsMore = .TRUE.
      ELSE
        DO j=1,4
          IF(CutDof(nIndexes(j))) EXIT
        END DO
        j2 = j
        IF( .NOT. CutDof(nn + eIndexes(j2) ) ) THEN
          j2 = MODULO(j-2,4)+1          
          IF(.NOT. CutDof(nn + eIndexes(j2))) THEN
            CALL Fatal('Caller','Could not imagine this 404 case!')
          END IF
        END IF
        pElement % NodeIndexes(1) = nIndexes(j)
        pElement % NodeIndexes(2) = nn + eIndexes(j2)        
        IsMore = .FALSE.
      END IF
        
    CASE( 40411 )      
      DO j=1,4
        IF(CutDof( nn + eIndexes(j))) EXIT
      END DO
      DO j2=1,4
        IF(CutDof( nIndexes(j2))) EXIT
      END DO
      pElement % NodeIndexes(1) = nn + eIndexes(j)
      pElement % NodeIndexes(2) = nIndexes(j2)

    CASE( 40402 )      
      DO j=1,4
        IF(CutDof( nIndexes(j))) EXIT
      END DO
      DO j2=1,4
        IF(j2==j) CYCLE
        IF(CutDof( nIndexes(j2))) EXIT
      END DO
      pElement % NodeIndexes(1) = nIndexes(j)
      pElement % NodeIndexes(2) = nIndexes(j2)

    CASE( 50440 )      
      DO j=1,6
        IF(CutDof( nn + eIndexes(j))) EXIT
      END DO
      DO j2=1,6
        IF(j2==j) CYCLE
        IF(CutDof( nn + eIndexes(j2))) EXIT
      END DO
      DO j3=1,6
        IF(j3==j .OR. j3==j2) CYCLE
        IF(CutDof( nn + eIndexes(j3))) EXIT
      END DO
      DO j4=1,6
        IF(j4==j .OR. j4==j2 .OR. j4==j3) CYCLE
        IF(CutDof( nn + eIndexes(j4))) EXIT
      END DO

      IF(m==1) THEN
        pElement % NodeIndexes(1) = nn + eIndexes(j)
        pElement % NodeIndexes(2) = nn + eIndexes(j2)
        pElement % NodeIndexes(3) = nn + eIndexes(j3)
        IsMore = .TRUE.
      ELSE
        pElement % NodeIndexes(1) = nn + eIndexes(j)
        pElement % NodeIndexes(2) = nn + eIndexes(j2)
        pElement % NodeIndexes(3) = nn + eIndexes(j4)
        IsMore = .FALSE.        
      END IF

    CASE( 50430 )      
      DO j=1,6
        IF(CutDof( nn + eIndexes(j))) EXIT
      END DO
      DO j2=1,6
        IF(j2==j) CYCLE
        IF(CutDof( nn + eIndexes(j2))) EXIT
      END DO
      DO j3=1,6
        IF(j3==j .OR. j3==j2) CYCLE
        IF(CutDof( nn + eIndexes(j3))) EXIT
      END DO
      pElement % NodeIndexes(1) = nn + eIndexes(j)
      pElement % NodeIndexes(2) = nn + eIndexes(j2)
      pElement % NodeIndexes(3) = nn + eIndexes(j3)

    CASE( 50421 )      
      DO j=1,6
        IF(CutDof( nn + eIndexes(j))) EXIT
      END DO
      DO j2=1,6
        IF(j2==j) CYCLE
        IF(CutDof( nn + eIndexes(j2))) EXIT
      END DO
      DO j3=1,4
        IF(CutDof( nIndexes(j3))) EXIT
      END DO
      pElement % NodeIndexes(1) = nn + eIndexes(j)
      pElement % NodeIndexes(2) = nn + eIndexes(j2)
      pElement % NodeIndexes(3) = nIndexes(j3)

    CASE( 50412 )      
      DO j=1,6
        IF(CutDof( nn + eIndexes(j))) EXIT
      END DO
      DO j2=1,4
        IF(CutDof( nIndexes(j2))) EXIT
      END DO
      DO j3=1,4
        IF(j3==j2) CYCLE
        IF(CutDof( nIndexes(j3))) EXIT
      END DO
      pElement % NodeIndexes(1) = nn + eIndexes(j)
      pElement % NodeIndexes(2) = nIndexes(j2)
      pElement % NodeIndexes(3) = nIndexes(j3)

    CASE( 50403 )      
      i = 0
      DO j=1,4
        ! We cut all other edges expect "j"
        IF(.NOT. CutDof( nIndexes(j))) EXIT
        i = i+1
        pElement % NodeIndexes(i) = nIndexes(j)
      END DO
      
      
    CASE( 70620 )      
      ! For prisms we currently assumes that the field is cut vertically such that
      ! we can split the bottom triangle and copy the same split in 3D.
      DO j=1,3
        IF( CutDof( nn + eIndexes(j) ) ) EXIT
      END DO
      DO j2=1,3
        IF(j2==j) CYCLE
        IF( CutDof( nn + eIndexes(j2) ) ) EXIT
      END DO
      pElement % NodeIndexes(1) = nn + eIndexes(j)
      pElement % NodeIndexes(2) = nn + eIndexes(j2)
      pElement % NodeIndexes(3) = nn + eIndexes(3+j2)
      pElement % NodeIndexes(4) = nn + eIndexes(3+j)

    CASE( 70611 ) 
      ! Find the edge and node that is cut
      DO j=1,3
        IF( CutDof( nn + eIndexes(j) ) ) EXIT
      END DO
      DO j2=1,3
        IF( CutDof( nIndexes(j2) ) ) EXIT
      END DO
      pElement % NodeIndexes(1) = nn + eIndexes(j)
      pElement % NodeIndexes(2) = nIndexes(j2)
      pElement % NodeIndexes(3) = nIndexes(j2+3)
      pElement % NodeIndexes(4) = nn + eIndexes(j+3)

    CASE( 70602 )       
      ! Find the two nodes that are cut.
      DO j=1,3
        IF( CutDof( nIndexes(j) ) ) EXIT
      END DO
      DO j2=1,3
        IF(j2==j) CYCLE
        IF( CutDof( nIndexes(j2) ) ) EXIT
      END DO
      pElement % NodeIndexes(1) = nIndexes(j)
      pElement % NodeIndexes(2) = nIndexes(j2)
      pElement % NodeIndexes(2) = nIndexes(j2+3)
      pElement % NodeIndexes(1) = nIndexes(j+3)

    CASE( 80830 )      
      DO j=1,12
        IF( CutDof( nn + eIndexes(j) ) ) EXIT
      END DO
      DO j2=1,12
        IF(j2==j) CYCLE
        IF( CutDof( nn + eIndexes(j2) ) ) EXIT
      END DO
      DO j3=1,12
        IF(j3==j .OR. j3==j2) CYCLE
        IF( CutDof( nn + eIndexes(j3) ) ) EXIT
      END DO
      pElement % NodeIndexes(1) = nn + eIndexes(j)
      pElement % NodeIndexes(2) = nn + eIndexes(j2)
      pElement % NodeIndexes(3) = nn + eIndexes(j3)

    CASE DEFAULT
      PRINT *,'Unknown SplitCase:',SplitCase
      PRINT *,'EdgeCut:',CutDof(nn + Element % EdgeIndexes) 
      PRINT *,'NodeCut:',CutDof(Element % NodeIndexes)
      PRINT *,'Phi:',PhiValues(PhiPerm(Element % NodeIndexes))
      CALL Fatal(Caller,'Unknown split case in bc element divisions: '//I2S(SplitCase))
    END SELECT


    ! This is just a tentative routine where we orient the nodes of the segment such that
    ! inside/outside is always consistently on left/right. 
    IF(pElement % TYPE % ElementCode == 202 ) THEN
      BLOCK
        REAL(KIND=dp) :: pmax, p, x0, x1, xp, y0, y1, yp, dir1, dir2
        INTEGER :: i,j,imax

        ! The most trustworthy point to define the sign of the levelset is the one with extreme value. 
        pmax = 0.0_dp
        imax = 0
        DO i=1,Element % TYPE % NumberOfNodes
          j = PhiPerm(Element % NodeIndexes(i))
          IF(j==0) CYCLE
          p = PhiValues(j)
          IF(ABS(p) > ABS(pmax)) THEN
            pmax = p
            imax = Element % NodeIndexes(i)
          END IF
        END DO
          
        ! Dir is an indicator one which side of the line segment the point lies. 
        x0 = Mesh % Nodes % x(pElement % NodeIndexes(1))
        y0 = Mesh % Nodes % y(pElement % NodeIndexes(1))
        x1 = Mesh % Nodes % x(pElement % NodeIndexes(2))
        y1 = Mesh % Nodes % y(pElement % NodeIndexes(2))
        xp = Mesh % Nodes % x(imax)
        yp = Mesh % Nodes % y(imax)

        dir1 = (x1 - x0) * (yp - y0) - (y1 - y0) * (xp - x0)

        ! Switch the signs so that if the point was found from left/right side of the
        ! line segment the sign stays the same as previously.
        IF(dir1 * pmax < 0.0_dp) THEN
          j = pElement % NodeIndexes(1)
          pElement % NodeIndexes(1) = pElement % NodeIndexes(2)
          pElement % NodeIndexes(2) = j
        END IF

#if DEBUG_ORIENT
        ! Here we can check that the orientation of the edges is consistent.
        ! This check only applies to convex geometries where the centermost node
        ! can be used to check the orientation. 
        x0 = Mesh % Nodes % x(pElement % NodeIndexes(1))
        y0 = Mesh % Nodes % y(pElement % NodeIndexes(1))
        x1 = Mesh % Nodes % x(pElement % NodeIndexes(2))
        y1 = Mesh % Nodes % y(pElement % NodeIndexes(2))
        dir2 = (x1 - x0) * (CutFEMCenter(2) - y0) - (y1 - y0) * (CutFEMCenter(1) - x0)

        IF( dir2 > 0.0 ) THEN
          PRINT *,'WrongDirection:',SplitCase,m,x0,x1,y0,y1
          PRINT *,'WrongDirIndexes:',CutDof(nIndexes),'e',CutDof(nn+eIndexes)
          PRINT *,'WrongDirPhi:',PhiValues(PhiPerm(nIndexes))
          PRINT *,'WrongDirX:',Mesh % Nodes % x(nIndexes)
          PRINT *,'WrongDirY:',Mesh % Nodes % y(nIndexes)
          PRINT *,'WrongDirImax:',imax,pmax

          STOP
          j = pElement % NodeIndexes(1)
          pElement % NodeIndexes(1) = pElement % NodeIndexes(2)
          pElement % NodeIndexes(2) = j
        END IF
#endif
        
      END BLOCK
    END IF
    
    pElement % BoundaryInfo % Left => NULL() ! Element
    pElement % BodyId = 0
    
  END FUNCTION CutInterfaceBC


  ! This is currently not used. 
  !-------------------------------------------------------
  SUBROUTINE CutFEMElementCount(Solver, Perm, nBulk, nBC ) 
    TYPE(Solver_t) :: Solver
    INTEGER, POINTER :: Perm(:)
    INTEGER :: nBulk, nBC
    
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: Active, t, n, nBulk0, nBC0, t0
    TYPE(Element_t), POINTER :: Element, pElement
    LOGICAL :: isCut, isMore, isActive
    
    nBulk = 0
    nBulk0 = 0
    nBC = 0
    nBC0 = 0
    Mesh => Solver % Mesh
    
    DO t=1,Mesh % NumberOfBulkElements
      Element => Mesh % Elements(t)
      IF(ANY(PhiPerm(Element % NodeIndexes)==0)) CYCLE
      CALL CutInterfaceCheck( Element, IsCut, IsActive, Perm )
      IF(.NOT. IsActive) CYCLE      
      IF(IsCut) THEN
10      pElement => CutInterfaceBulk(Element,isCut,isMore)        
        IF(ALL(Perm(pElement % NodeIndexes) > 0) ) nBulk = nBulk + 1
        IF(IsMore) GOTO 10
      ELSE        
        nBulk0 = nBulk0 + 1
      END IF
    END DO
    

    ! Additional BC elements created on the interface. 
    DO t=1,Mesh % NumberOfBulkElements
      Element => Mesh % Elements(t)
      IF(ANY(PhiPerm(Element % NodeIndexes)==0)) CYCLE            
      CALL CutInterfaceCheck( Element, IsCut, IsActive, Perm )
      IF(.NOT. IsActive) CYCLE
20    pElement => CutInterfaceBC(Element,isCut,isMore)        
      IF(ASSOCIATED(pElement)) THEN          
        IF(ALL(Perm(pElement % NodeIndexes) > 0) ) nBC = nBC + 1
        IF(IsMore) GOTO 20
      END IF
    END DO

    ! Remaining original boundary element.
    t0 = Mesh % NumberOfBulkElements
    DO t=1,Mesh % NumberOfBoundaryElements
      Element => Mesh % Elements(t0+t)
      IF(ANY(PhiPerm(Element % NodeIndexes)==0)) CYCLE
      IF(ALL(Perm(Element % NodeIndexes) > 0) ) nBC0 = nBC0 + 1
    END DO
        
    CALL Info('CutFEMElementCount','Bulk elements remaining '//I2S(nBulk0)//' & splitted '//I2S(nBulk),Level=7)
    CALL Info('CutFEMElementCount','BC elements remaining '//I2S(nBC0)//' & splitted '//I2S(nBC),Level=7)
    
    nBC = nBC0 + nBC
    nBulk = nBulk0 + nBulk
    
  END SUBROUTINE CutFEMElementCount
    
  
  SUBROUTINE CreateCutFEMAddMesh(Solver) 
    TYPE(Solver_t) :: Solver
    
    INTEGER :: Sweep, t, n, i
    LOGICAL :: IsActive, IsCut
    TYPE(Element_t), POINTER :: Element

    CALL Info('CreateCutFEMAddMesh','Creating interface mesh from split element',Level=10)
    
    DO Sweep = 0,1 
      n = 0      
      DO t=1,CutFEMOrigMesh % NumberOfBulkElements
        Element => CutFEMOrigMesh % Elements(t)
        IF(ANY(PhiPerm(Element % NodeIndexes)==0)) CYCLE
        CALL CutInterfaceCheck( Element, IsCut, IsActive, CutPerm )
        IF(IsActive .AND. .NOT. IsCut) THEN
          n = n+1
          IF(Sweep==1) UnsplitActiveElements(n) = t
        END IF
      END DO
      IF(Sweep == 0) THEN
        ALLOCATE(UnsplitActiveElements(n))
      END IF
    END DO

    IF(ASSOCIATED(CutFEMAddMesh)) THEN
      NULLIFY(CutFEMAddMesh % Nodes % x)
      NULLIFY(CutFEMAddMesh % Nodes % y)
      NULLIFY(CutFEMAddMesh % Nodes % z)
      CALL FreeMesh(CutFEMAddMesh)
    END IF
    CutFEMAddMesh => CreateCutFEMMesh(Solver,CutFEMOrigMesh,CutPerm,&
        .TRUE.,.TRUE.,.TRUE.,Solver % Values,'dummy variable') 
    
    CALL MeshStabParams( CutFEMAddMesh )
    
    n = CutFEMAddMesh % NumberOfBulkElements
    ALLOCATE(AddActiveElements(n))
    DO i=1,n
      AddActiveElements(i) = i
    END DO

    Solver % ActiveElements => UnsplitActiveElements
    Solver % NumberOfActiveElements = SIZE(UnsplitActiveElements)

    CALL Info('CreateCutFEMAddMesh','Add mesh created with '//I2S(i)//' active elements!',Level=10)
    
  END SUBROUTINE CreateCutFEMAddMesh
    
  SUBROUTINE CutFEMSetAddMesh(Solver)
    TYPE(Solver_t) :: Solver

    Solver % Mesh => CutFEMAddMesh
    CurrentModel % Mesh => CutFEMAddMesh
    Solver % ActiveElements => AddActiveElements
    Solver % NumberOfActiveElements = SIZE(Solver % ActiveElements)
    Solver % Mesh % Edges => CutFemOrigMesh % Edges

    CALL Info('CutFEMSetAddMesh','Swapping CutFEM original mesh to interface mesh!',Level=10)
    
  END SUBROUTINE CutFEMSetAddMesh
    
  SUBROUTINE CutFEMSetOrigMesh(Solver)
    TYPE(Solver_t) :: Solver
    
    Solver % Mesh => CutFEMOrigMesh
    CurrentModel % Mesh => CutFEMOrigMesh 
    Solver % ActiveElements => UnsplitActiveElements
    Solver % NumberOfActiveElements = SIZE(Solver % ActiveElements)

    NULLIFY(CutFEMAddMesh % Edges) 

    CALL Info('CutFEMSetOrigMesh','Swapping CutFEM interface mesh to original mesh!',Level=10)
    
  END SUBROUTINE CutFEMSetOrigMesh
    


  
  ! Assembly a matrix for extrapolating values outside the active domain.
  ! Currently this is just diffusion matrix. We could perhaps use convection also. 
  !------------------------------------------------------------------------------
  SUBROUTINE LocalFitMatrix( Element, n )
    !------------------------------------------------------------------------------
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: weight, dcoeff 
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ
    REAL(KIND=dp) :: STIFF(n,n), FORCE(n), LOAD(n)
    LOGICAL :: Stat,Found,CutElem
    INTEGER :: i,j,t,p,q
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes

    SAVE Nodes
    !------------------------------------------------------------------------------

!    CALL GetElementNodes( Nodes, Element )
    CALL CopyElementNodesFromMesh( Nodes, CurrentModel % Solver % Mesh, n, Element % NodeIndexes)

    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0.0_dp

    dcoeff = 1.0_dp

    ! Numerical integration:
    !-----------------------
    IP = GaussPointsAdapt( Element )
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx )      
      IF(.NOT. Stat) CYCLE
      Weight = IP % s(t) * DetJ

      STIFF(1:n,1:n) = STIFF(1:n,1:n) + dcoeff * Weight * &
          MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )      
    END DO
    
    CutElem = .TRUE.
!    CALL DefaultUpdateEquations(STIFF,FORCE,CutElem,Element)

    !------------------------------------------------------------------------------
  END SUBROUTINE LocalFitMatrix
  !------------------------------------------------------------------------------



  ! This takes a CutFEM variable and either extrapolates it using a FE equation
  ! (for many element layers) or just extends it by extrapolating on the cut edges.  
  !---------------------------------------------------------------------------------  
  SUBROUTINE CutFEMVariableFinalize( Solver ) 
    TYPE(Solver_t) :: Solver
    
    TYPE(Matrix_t), POINTER :: B
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: pElement, Element
    INTEGER :: i,j,k,l,n,t,active,nn,ne,i1,i2,dofs
    LOGICAL :: IsCut, IsMore, Found
    REAL(KIND=dp) :: s, r, dval, norm
    REAL(KIND=dp), ALLOCATABLE :: NodeWeigth(:)
    

    Mesh => Solver % Mesh
    nn = Mesh % NumberOfNodes
    ne = Mesh % NumberOfEdges
    dofs = CutDofs 
    
    ! If we solve some other equation in between store the original norm.
    Norm = Solver % Variable % Norm

    ! Set values at shared nodes that have been computed. 
    CALL Info('CutFEMVariableFinalize','Copying values at shared nodes to the original mesh!',Level=10)
    DO l=1,dofs
      DO i=1,nn
        j = CutPerm(i)
        k = OrigMeshPerm(i)
        IF(j==0 .OR. k==0) CYCLE
        OrigMeshValues(dofs*(k-1)+l) = CutValues(dofs*(j-1)+l)
      END DO
    END DO

    
    ! We can only extrapolate using the edges that are cut since they have also one
    ! known nodal value. 
    IF(CutExtrapolate) THEN
      CALL Info('CutFEMVariableFinalize','Extrapolating values with split elements!',Level=10)
      
      ! Extrapolated nodes may have more than one hit. Hence use weigted average.
      ! This is the weight.
      ALLOCATE(NodeWeigth(nn))
      NodeWeigth = 0.0_dp

      k = 0
      DO i=1,Solver % Mesh % NumberOfEdges
        j = CutPerm(nn+i)      
        IF(j==0) CYCLE
        r = CutInterp(i)
        
        i1 = Mesh % Edges(i) % NodeIndexes(1)
        i2 = Mesh % Edges(i) % NodeIndexes(2)

        IF(CutPerm(i1) > 0 .AND. CutPerm(i2) == 0 ) THEN
          s = (1-r)
          DO k=1,dofs
            OrigMeshValues(dofs*(OrigMeshPerm(i2)-1)+k) = OrigMeshValues(dofs*(OrigMeshPerm(i2)-1)+k) + &
                s*CutValues(dofs*(CutPerm(i1)-1)+k) + (CutValues(dofs*(j-1)+k)-CutValues(dofs*(CutPerm(i1)-1)+k))
          END DO
          NodeWeigth(OrigMeshPerm(i2)) = NodeWeigth(OrigMeshPerm(i2)) + s
        ELSE IF(CutPerm(i1) == 0 .AND. CutPerm(i2) > 0) THEN
          s = r
          DO k=1,dofs
            OrigMeshValues(dofs*(OrigMeshPerm(i1)-1)+k) = OrigMeshValues(dofs*(OrigMeshPerm(i1)-1)+k) + &
                s*CutValues(dofs*(CutPerm(i2)-1)+k) + (CutValues(dofs*(j-1)+k)-CutValues(dofs*(CutPerm(i2)-1)+k))
          END DO
          NodeWeigth(OrigMeshPerm(i1)) = NodeWeigth(OrigMeshPerm(i1)) + s
        END IF
      END DO
      
      DO k=1,dofs
        WHERE( NodeWeigth(1:nn) > EPSILON(s)) 
          OrigMeshValues(k::dofs) = OrigMeshValues(k::dofs) / NodeWeigth(1:nn)
        END WHERE
      END DO
    END IF

    
    ! Extend values using FEM strategies beyond value set above. 
    ! We can extrapolate much but the extrapolation method is nonhysical.
    IF( CutExtend ) THEN
      CALL Info('CutFEMVariableFinalize','Extending values from inside to outside using FEM!',Level=10)
      IF(dofs > 1) THEN
        CALL Fatal('CutFEMVariableFinalize','Extending values only coded for one dofs!')        
      END IF
      
      B => CreateCutFEMMatrix(Solver,ExtendPerm)      
      ALLOCATE(ExtendValues(B % NumberOfRows))
      ExtendValues = 0.0_dp
      Solver % Matrix => B
      Solver % Variable % Values => ExtendValues
      Solver % Variable % Perm => ExtendPerm

      DO t=1,Mesh % NumberOfBulkElements
        Element => Mesh % Elements(t)
        IF(ANY(PhiPerm(Element % NodeIndexes)==0)) CYCLE
        
        n  = Element % Type % NumberOfNodes

30      pElement => CutInterfaceBulk(Element,isCut,isMore)        
        IF(isCut) THEN          
          n  = pElement % Type % NumberOfNodes
          IF(ALL(ExtendPerm(pElement % NodeIndexes) > 0) ) THEN
            CALL LocalFitMatrix( pElement, n )
          END IF
          IF(IsMore) GOTO 30
        ELSE
          IF(ALL(ExtendPerm(Element % NodeIndexes) > 0) ) THEN
            CALL LocalFitMatrix( Element, n )
          END IF
        END IF
      END DO


      ! On the shared nodes of the "inside" and "outside" regions. 
      ! Set Dirichlet conditions in somewhat dirty way for now. 
      DO i=1,Mesh % NumberOfNodes + Mesh % NumberOfEdges
        j = CutPerm(i)
        k = ExtendPerm(i)
        IF(j==0 .OR. k==0) CYCLE
        
        dval = CutValues(j)
        s = B % Values(B % diag(k))
        CALL ZeroRow(B, k)
        CALL SetMatrixElement(B,k,k,s)
        B % rhs(k) = s * dval
      END DO
      
      CALL SolveLinearSystem( B, B % rhs, ExtendValues, norm, Solver % Variable % dofs, Solver ) 
      
      CALL FreeMatrix(B)
      Solver % Matrix => NULL()
      
      DO i=1,nn
        j = ExtendPerm(i)
        k = OrigMeshPerm(i)
        IF(j==0 .OR. k==0) CYCLE
        OrigMeshValues(k) = ExtendValues(j)
      END DO
    END IF
    
    ! Revert to the original field that is present everywhere.           
    Solver % Variable % Values => OrigMeshValues
    Solver % Variable % Perm => OrigMeshPerm
    Solver % Variable % PrevValues => OrigPrevMeshValues
    Solver % Variable % Norm = Norm
    
    ! Revert to original body id's.
    ! If we don't do this then ActiveElements is spoiled. 
    DO t=1,Mesh % NumberOfBulkElements        
      Element => Mesh % Elements(t)
      IF(ALL(PhiPerm(Element % NodeIndexes)>0)) THEN
        Element % BodyId = CutFemBody
      END IF
    END DO
    
  END SUBROUTINE CutFEMVariableFinalize


!------------------------------------------------------------------------------
!> Split a mesh at zero levelset by adding new nodes at the interface.
!> The idea is to be able to better represent shapes that are not initially
!> presented by body fitted finite element mesh. This is a modifieid version
!> of similar routine in MeshUtils that utilizes the CutInterface* routines.
!------------------------------------------------------------------------------
  FUNCTION CreateCutFEMMesh(Solver,Mesh,Perm,CreateBC,CreateBulk,&
      AddMeshMode, Vlist,ProjectPrefix) RESULT( NewMesh )
!------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    TYPE(Mesh_t) :: Mesh
    INTEGER, POINTER :: Perm(:)
    LOGICAL :: CreateBC, CreateBulk, AddMeshMode
    TYPE(ValueList_t), POINTER :: Vlist
    CHARACTER(*) :: ProjectPrefix
    TYPE(Mesh_t), POINTER :: NewMesh
!------------------------------------------------------------------------------
    INTEGER :: i, j, k, n
    INTEGER :: NodeCnt
    INTEGER :: nn, ne, nBC, nBulk, t, ntot, Sweep, InterfaceBC
    LOGICAL :: Found, isActive, isMore, isCut
    TYPE(Element_t), POINTER :: pElement,Element
    REAL(KIND=dp) :: r
    INTEGER, POINTER :: MeshPerm(:) => NULL()
    REAL(KIND=dp), POINTER :: Values(:)
    CHARACTER(:), ALLOCATABLE :: VarName       
    CHARACTER(*), PARAMETER :: Caller = 'CreateCutFEMMesh'

    SAVE MeshPerm
    
!------------------------------------------------------------------------------
    IF(.NOT. (CreateBC .OR. CreateBulk)) THEN
      CALL Info(Caller,'Nothing to do!?')
      RETURN
    END IF

    IF( AddMeshMode ) THEN
      CALL Info( Caller, 'Creating mesh including splitted elements only!')
    ELSE IF(.NOT. CreateBulk ) THEN
      CALL Info( Caller, 'Creating mesh including isoline boundary elements only!')
    ELSE 
      CALL Info( Caller, 'Creating actual mesh splitted by zero levelset!')
    END IF
      

    CALL ResetTimer(Caller)

    ! Define the nodes to be included in the new mesh.
    !----------------------------------
    nn = Mesh % NumberOfNodes
    ne = Mesh % NumberOfEdges

    IF(.NOT. ( CreateBulk .OR. AddMeshMode ) ) THEN
      ALLOCATE(MeshPerm(nn+ne))
      MeshPerm = 0    
      j = 0
      DO i=1,nn+ne
        IF(CutDof(i)) THEN
          IF(AddMeshMode) THEN
            MeshPerm(i) = i
          ELSE
            j=j+1
            MeshPerm(i) = j
          END IF
        END IF
      END DO
    ELSE
      MeshPerm => Perm      
    END IF

    NewMesh => AllocateMesh()    
    NewMesh % SingleMesh = Mesh % SingleMesh
    NewMesh % MaxNDofs = Mesh % MaxNDofs
    NewMesh % MeshDim = Mesh % MeshDim
    NewMesh % MaxElementNodes = Mesh % MaxElementNodes
    
    IF( AddMeshMode ) THEN
      ! In add mesh mode we retain the nodes and coordinates of the original mesh
      ! and just create the elements and their topologies.
      ! This saves work and memory. 
      NewMesh % Name = TRIM(Mesh % Name)//'-addmesh'
      NodeCnt = Mesh % NumberOfNodes
      NewMesh % Nodes % x => Mesh % Nodes % x
      NewMesh % Nodes % y => Mesh % Nodes % y
      NewMesh % Nodes % z => Mesh % Nodes % z
    ELSE
      ! This mode is intended for the isoline that is created at the levelset.
      NewMesh % Name = TRIM(Mesh % Name)//'-cutfem'
      NodeCnt = MAXVAL(MeshPerm)
      NewMesh % OutputActive = .TRUE.

      CALL AllocateVector( NewMesh % Nodes % x, NodeCnt ) 
      CALL AllocateVector( NewMesh % Nodes % y, NodeCnt ) 
      CALL AllocateVector( NewMesh % Nodes % z, NodeCnt ) 

      ! This includes already the nodes created at the intersections. 
      DO i=1,nn+ne
        j = MeshPerm(i)
        IF(j==0) CYCLE      
        NewMesh % Nodes % x(j) = Mesh % Nodes % x(i)
        NewMesh % Nodes % y(j) = Mesh % Nodes % y(i)
        NewMesh % Nodes % z(j) = Mesh % Nodes % z(i)
      END DO
    END IF
      
    CALL Info(Caller,'Number of nodes in CutFEM mesh: '//I2S(NodeCnt),Level=6)

    NewMesh % NumberOfNodes = NodeCnt
    NewMesh % Nodes % NumberOfNodes = NodeCnt

    InterfaceBC = ListGetInteger( Solver % Values,'CutFEM Interface BC',Found )    
        
    ! The 1st cycle just compute the number of elements.
    ! In between allocate the mesh elements.
    ! The 2nd cycle add the detected elements to the list.
    !------------------------------------------------------
    DO Sweep=0,1
      nBulk = 0
      nBC = 0
      
      IF(CreateBulk) THEN
        DO t=1, Mesh % NumberOfBulkElements       
          Element => Mesh % Elements(t)      
          CALL CutInterfaceCheck( Element, IsCut, IsActive, Perm )          
          !IF(.NOT. IsActive) CYCLE      
          IF(IsCut) THEN
10          pElement => CutInterfaceBulk(Element,isCut,isMore)        
            IF(ALL(Perm(pElement % NodeIndexes) > 0) ) THEN
              nBulk = nBulk + 1
              IF(Sweep==1) CALL AddElementData(pElement,nBulk)
            END IF
            IF(IsMore) GOTO 10
          ELSE IF(.NOT. AddMeshMode ) THEN       
            ! We we create only interface then the standard bulk elements are not included!
            IF(ANY(Perm(Element % NodeIndexes) == 0) ) CYCLE
            nBulk = nBulk + 1
            IF(Sweep==1) CALL AddElementData(Element,nBulk)
          END IF
        END DO
      END IF

      IF(CreateBC) THEN
        DO t=1,Mesh % NumberOfBulkElements
          Element => Mesh % Elements(t)
          !CALL CutInterfaceCheck( Element, IsCut, IsActive, Perm )
          !IF(.NOT. IsActive) CYCLE
20        pElement => CutInterfaceBC(Element,isCut,isMore)        
          IF(isCut) THEN
            IF(ASSOCIATED(pElement)) THEN          
              IF(ASSOCIATED(Perm)) THEN
                IF(ALL(Perm(pElement % NodeIndexes) > 0) ) THEN
                  nBC = nBC + 1
                  IF(Sweep==1) CALL AddElementData(pElement,nBulk+nBC,InterfaceBC)
                END IF
              END IF
              IF(IsMore) GOTO 20
            END IF
          END IF
        END DO
      END IF

      
      IF( Sweep == 0 ) THEN
        IF( CreateBulk ) THEN
          NewMesh % NumberOfBulkElements = nBulk
          NewMesh % NumberOfBoundaryElements = nBC
        ELSE
          NewMesh % NumberOfBulkElements = nBC
          NewMesh % NumberOfBoundaryElements = 0
        END IF

        IF(InfoActive(25)) THEN
          PRINT *,'Old Element Counts:',Mesh % NumberOfBulkElements, Mesh % NumberOfBoundaryElements
          PRINT *,'New Element Counts:',nBulk, nBC
        END IF
          
        CALL AllocateVector( NewMesh % Elements, nBulk+nBC )
        CALL Info(Caller,'New mesh allocated for '//I2S(nBulk+nBc)//' elements', Level=10 )
      END IF

    END DO

#if 1
    ! if add mesh mode we can just use oldparallel structures
    IF( ParEnv % PEs > 1 .AND. .NOT. AddMeshMode ) CALL CutFEMParallelMesh()
#endif
      
    
    ! If we create interface only then we have original numbering and may use 
    IF(.NOT. AddMeshMode ) THEN 
      CALL InterpolateLevelsetVariables()
    END IF
      
    IF(.NOT. ( CreateBulk .OR. AddMeshMode) ) DEALLOCATE(MeshPerm)

    CALL CheckTimer(Caller,Delete=.TRUE.)     
    CALL Info(Caller,'Zero levelset mesh was created',Level=8)

  CONTAINS

#if 1
    ! We do not need to update the Mesh % ParallelInfo, only Matrix % ParallelInfo!
    
    SUBROUTINE CutFEMParallelMesh()

      INTEGER :: istat,n0,n
          
      CALL Info(Caller,'Creating ParallelInfo for CutFEM mesh structures!',Level=10)      
      IF(.NOT. ASSOCIATED(Mesh % ParallelInfo % GlobalDOFS) ) THEN
        CALL Fatal(Caller,'Original mesh has no GlobalDOFs numbering!')
      END IF
      IF(.NOT. ASSOCIATED(Mesh % Edges) ) THEN
        CALL Fatal(Caller,'Original mesh requires edges!')
      END IF
      
      ! Use maximum nodal index as the offset for nodes defined on cut edges.
      n0 = MAXVAL( Mesh % ParallelInfo % GlobalDOFs )
      n0 = ParallelReduction(n0,2)

      n = NewMesh % NumberOfNodes
      CALL Info(Caller,'Allocating parallel structures for '//I2S(n)//' nodes',Level=10)

      ALLOCATE(NewMesh % ParallelInfo % GlobalDOFs(n), STAT=istat )
      IF ( istat /= 0 ) &
          CALL Fatal( Caller, 'Unable to allocate NewMesh % ParallelInfo % NeighbourList' )
      NewMesh % ParallelInfo % GlobalDOFs = 0
      ALLOCATE(NewMesh % ParallelInfo % GInterface(n), STAT=istat )
      IF ( istat /= 0 ) &
          CALL Fatal( Caller, 'Unable to allocate NewMesh % ParallelInfo % NeighbourList' )
      NewMesh % ParallelInfo % GInterface = .FALSE.

      ALLOCATE(NewMesh % ParallelInfo % NeighbourList(n), STAT=istat )
      IF ( istat /= 0 ) &
          CALL Fatal( Caller, 'Unable to allocate NewMesh % ParallelInfo % NeighbourList' )
      DO i=1,n
        NULLIFY(NewMesh % ParallelInfo % NeighbourList(i) % Neighbours)
      END DO      
      
      DO i=1,nn+ne
        j = MeshPerm(i)
        IF(j<=0) CYCLE
        
        IF(i<=nn) THEN
          NewMesh % ParallelInfo % GInterface(j) = Mesh % ParallelInfo % GInterface(i)
          NewMesh % ParallelInfo % GlobalDOFs(j) = Mesh % ParallelInfo % GlobalDOFs(i)
          k = SIZE(Mesh % ParallelInfo % NeighbourList(i) % Neighbours)
          ALLOCATE(NewMesh % ParallelInfo % NeighbourList(j) % Neighbours(k))
          NewMesh % ParallelInfo % NeighbourList(j) % Neighbours = &
              Mesh % ParallelInfo % NeighbourList(i) % Neighbours            
        ELSE
          NewMesh % ParallelInfo % GInterface(j) = Mesh % ParallelInfo % EdgeInterface(i-nn)
          NewMesh % ParallelInfo % GlobalDOFs(j) = n0 + Mesh % Edges(i-nn) % GElementIndex         

          k = SIZE(Mesh % ParallelInfo % EdgeNeighbourList(i-nn) % Neighbours)
          !PRINT *,'ass1 vals:',ParEnv % MyPe, k,j,Mesh % ParallelInfo % EdgeNeighbourList(i-nn) % Neighbours
          ALLOCATE(NewMesh % ParallelInfo % NeighbourList(j) % Neighbours(k))
          NewMesh % ParallelInfo % NeighbourList(j) % Neighbours = &
              Mesh % ParallelInfo % EdgeNeighbourList(i-nn) % Neighbours                                
        END IF
      END DO


      DO i = 1, NewMesh % NumberOfNodes
        IF(.NOT. ASSOCIATED(NewMesh % ParallelInfo % NeighbourList(i) % Neighbours)) THEN
          !PRINT *,ParEnv % MYPE, 'nn:',nn,ne, MAXVAL(MeshPerm), NewMesh % NumberOfNodes, &
          !    SIZE(MeshPerm)
          CALL Fatal('CutFEMParallelMesh','Neighbours not associated: '//I2S(i))
        END IF
      END DO
      
    END SUBROUTINE CutFEMParallelMesh

#endif        

    
    ! We can easily interpolate any variable on the new nodes created on the edge. 
    ! if we know values on both nodes. 
    !-----------------------------------------------------------------------------
    SUBROUTINE InterpolateLevelsetVariables()

      INTEGER :: iVar, dofs, l
      TYPE(Variable_t), POINTER :: Var
      REAL(KIND=dp), POINTER :: Values(:)
      INTEGER, POINTER :: Perm(:)
      LOGICAL :: IsCutVar
      
      DO iVar = -1,100    
        
        IF(iVar == -1) THEN
          ! We want to always interpolate the primary variable!
          ! This is the only variable living in the "CutFEM" universe.
          Var => Solver % Variable
          VarName = Solver % Variable % name
          IsCutVar = .TRUE.
        ELSE
          IF(iVar == 0) THEN
            ! We also want to interpolate the levelset variable. 
            VarName = CutStr
          ELSE                       
            VarName = ListGetString( Vlist,TRIM(ProjectPrefix)//' '//I2S(iVar), Found )
            IF(.NOT. Found ) EXIT    

            ! These are cases "-1" and "-0" that are always done!
            IF(VarName == Solver % Variable % Name ) CYCLE
            IF(VarName == CutStr ) CYCLE
          END IF
            
          Var => VariableGet( Mesh % Variables, VarName, ThisOnly = .TRUE. )
          IF(.NOT. ASSOCIATED(Var)) THEN
            CALL Fatal('InterpolateLevelsetVariable','Could not find variable in 2D mesh: '//TRIM(VarName))
          END IF
          IsCutVar = .FALSE.
        END IF
          
        CALL Info('InterpolateLevelsetVariable','Doing field: '//TRIM(Var % Name))

        dofs = Var % dofs
        NULLIFY(Values)        
        ALLOCATE(Values(dofs*NodeCnt))
        Values = 0.0_dp

        ! Create unity permutation vector. 
        NULLIFY(Perm)
        ALLOCATE(Perm(NodeCnt))
        DO i=1,NodeCnt
          Perm(i) = i
        END DO
          
        ! If the size of permutation is nn+ne it is a sign that the field is associated
        ! to the CutFEM field.
        IF(IsCutVar) THEN
          ntot = nn + ne
        ELSE
          ntot = nn 
        END IF
        
        DO i=1,ntot
          j = Var % Perm(i)
          k = MeshPerm(i)
          IF(j==0 .OR. k==0) CYCLE
          DO l=1,dofs
            Values(dofs*(k-1)+l) = Var % Values(dofs*(j-1)+l)
          END DO
        END DO

        ! We do not want to interpolate the cut var that has been computed exactly to the new virtual nodes.
        ! The interpolated values would be worse than the computed ones. 
        IF(.NOT. IsCutVar ) THEN
          CALL Info(Caller,'Interpolating values: '//TRIM(VarName),Level=10)
          DO i=1,ne
            k = MeshPerm(nn+i)
            IF(k==0) CYCLE
            r = CutInterp(i)

            Values(k) = 0.0_dp

            j = Var % Perm(Mesh % Edges(i) % NodeIndexes(1))
            IF(j>0) THEN
              DO l=1,dofs
                Values(dofs*(k-1)+l) = r*Var % Values(dofs*(j-1)+l)
              END DO
            END IF
              
            j = Var % Perm(Mesh % Edges(i) % NodeIndexes(2))
            IF(j>0) THEN
              DO l=1,dofs
                Values(dofs*(k-1)+l) = Values(dofs*(k-1)+l) + (1-r)*Var % Values(dofs*(j-1)+l)
              END DO
            END IF
          END DO
        END IF

        CALL Info(Caller,'Projected variable: '//TRIM(VarName),Level=10)
        CALL VariableAddVector( NewMesh % Variables, NewMesh, Solver, VarName, Var % Dofs, Values, Perm )
        
        IF(InfoActive(25)) THEN
          PRINT *,'Range:',MINVAL(Values),MAXVAL(Values),SIZE(Values), Var % Dofs, SIZE(Values)
        END IF
      END DO
      
    END SUBROUTINE InterpolateLevelsetVariables


    ! Here we just add some data to the new mesh.
    !--------------------------------------------
    SUBROUTINE AddelementData(pElement,ElemInd,BCtag)
      TYPE(Element_t), POINTER :: pElement
      INTEGER :: ElemInd
      INTEGER, OPTIONAL :: BCTag
      
      TYPE(Element_t), POINTER :: Enew
      INTEGER :: n

      Enew => NewMesh % Elements(ElemInd)        
      Enew % PartIndex = pElement % PartIndex
      Enew % BodyId = pElement % BodyId
      Enew % ElementIndex = ElemInd

      n = pElement % TYPE % NumberOfNodes
      Enew % TYPE => GetElementType(pElement % TYPE % ElementCode)
      
      CALL AllocateVector( ENew % NodeIndexes, n)
      IF( AddMeshMode ) THEN
        Enew % NodeIndexes(1:n) = pElement % NodeIndexes(1:n)
      ELSE
        Enew % NodeIndexes(1:n) = MeshPerm(pElement % NodeIndexes(1:n))
      END IF
      Enew % NDOFs = n
      Enew % EdgeIndexes => NULL()
      Enew % FaceIndexes => NULL()

      IF(PRESENT(BCTag)) THEN
        ALLOCATE(Enew % BoundaryInfo)
        Enew % BoundaryInfo % Constraint = BCTag
      END IF

      ! This effects only parallel runs, but testing parallel costs more...      
      Enew % PartIndex = ParEnv % MyPe
      
    END SUBROUTINE AddelementData

  END FUNCTION CreateCutFEMMesh
!------------------------------------------------------------------------------



  ! Here we create a 1D mesh on the zero level-set and convent it to new location,
  ! and compute new signed distance. 
  !--------------------------------------------  
  SUBROUTINE LevelSetUpdate(Solver,Mesh)

    TYPE(Solver_t) :: Solver
    TYPE(Mesh_t) :: Mesh

    TYPE(Variable_t), POINTER :: PhiVar1D, PhiVar2D, pVar
    TYPE(Mesh_t), POINTER :: IsoMesh => NULL()
    REAL(KIND=dp), POINTER :: x(:), y(:)
    REAL(KIND=dp) :: val, Vx, Vy, dt, VPhi, PhiMax, BW, cosphi0
    REAL(KIND=dp), ALLOCATABLE :: NodeDiff(:),ValStore(:)
    CHARACTER(:), ALLOCATABLE :: str       
    LOGICAL :: Found, Nonzero, MovingLevelset, NormalMove, NodeHistory, Positive, CheckLS
    LOGICAL, ALLOCATABLE :: Trust(:), DoesElemIntersect(:)
    INTEGER :: nVar,i,j,iAvoid,iSolver,k,l,counter,NNeighbours,MyPE
    INTEGER, ALLOCATABLE :: LocalPerm(:),Neighbours(:),nSend(:),nRecv(:)
    INTEGER, POINTER :: NodeIndexes(:)
    
    TYPE PolylineData_t
      INTEGER :: nLines = 0, nNodes = 0
      REAL(KIND=dp), ALLOCATABLE :: Vals(:,:)
      INTEGER, ALLOCATABLE :: Prev(:,:), Next(:,:)
      LOGICAL, ALLOCATABLE :: Intersect(:)
      REAL(KIND=dp) :: IsoLineBB(4), MeshBB(4)
    END TYPE PolylineData_t
    TYPE(PolylineData_t),  ALLOCATABLE, TARGET, SAVE :: PolylineData(:)

    TYPE Intersects_t
      INTEGER :: nIntersects = 0, Size = 0
      REAL(KIND=dp), ALLOCATABLE :: Intersection(:,:), LS(:,:)
      INTEGER, ALLOCATABLE :: Elements(:,:)
      LOGICAL, ALLOCATABLE :: Found(:,:)
    END TYPE Intersects_t
    TYPE(Intersects_t),  ALLOCATABLE, TARGET, SAVE :: Intersects(:)
    
    
    SAVE IsoMesh

    IsoMesh => CreateCutFEMMesh(Solver,Mesh,Solver % Variable % Perm,&
        .TRUE.,.FALSE.,.FALSE.,Solver % Values,'isoline variable')     
    IsoMesh % Name = TRIM(Mesh % Name)//'-isomesh'
    
    pVar => VariableGet( Mesh % Variables,'timestep size' )
    dt = pVar % Values(1)

    phiVar1D => VariableGet( IsoMesh % Variables, CutStr, ThisOnly = .TRUE.)
    IF(.NOT. ASSOCIATED(PhiVar1D)) THEN
      CALL Fatal('LevelSetUpdate','Levelset function ("'//TRIM(CutStr)//'") needed in 1D mesh!')
    END IF

    ! This should be ok by construction but some testing does not hurt...
    DO i=1,SIZE(phiVar1D % Perm)
      IF(i /= PhiVar1D % Perm(i) ) THEN
        CALL Fatal('LevelSetUpdate','PhiVar1D permutation not unity map')
      END IF
    END DO
    
    phiVar2D => VariableGet( Mesh % Variables, CutStr, ThisOnly = .TRUE.)
    IF(.NOT. ASSOCIATED(PhiVar2D)) THEN
      CALL Fatal('LevelSetUpdate','Levelset function needed in 2D mesh!')
    END IF

    ! We can move the levelset either by moving its coordinates, or
    ! by adding values to the levelset function. The first should be used if we
    ! know the direction of the movement and the latter if we know the normal movement.
    x => Isomesh % Nodes % x    
    y => Isomesh % Nodes % y
    
    MovingLevelset = .FALSE.

    ! It turns out that if the polyline is not a zero levelset but something else
    ! we need to try to ensure that the direction is almost normal to the element segment.
    VPhi = ListGetCReal( Solver % Values,'CutFEM critical angle',Found )
    IF(.NOT. Found) Vphi = 60.0_dp
    cosphi0 = COS(Vphi*PI/180.0_dp)

    Nonzero = ListGetLogical( Solver % Values,'CutFEM signed distance nonzero',Found )

    NormalMove = ListGetLogical( Solver % Values,'CutFEM normal move',Found )

    NodeHistory = ListGetLogical( Solver % Values,'CutFEM node history',Found )

    CheckLS = ListGetLogical( Solver % Values,'CutFEM check ls field',Found )

    ALLOCATE(NodeDiff(Mesh % NumberOfNodes), Trust(Mesh % NumberOfNodes))
    NodeDiff = 0.0_dp; Trust = .FALSE.
    IF(NodeHistory) THEN
      IF(.NOT. ALLOCATED(PolylineData)) THEN
        ALLOCATE(PolylineData(ParEnv % PEs))
      END IF
      CALL PopulatePolyline()

      DO i=1, Mesh % NumberOfNodes
        j = PhiVar2D % Perm(i)
        IF(j==0 .AND. .NOT. NonZero) CYCLE
        IF(j==0) STOP

        val = SignedDistance(i, Trust(i))

        IF(Trust(i)) THEN
          NodeDiff(i) = PhiVar2D % Values(j) - val
        ELSE
          NodeDiff(i) = 0.0_dp
        END IF
      END DO

      ! Deallocate data, next time this will be different.
      DO i=1,ParEnv % PEs
        IF(PolylineData(i) % nLines > 0) THEN
          DEALLOCATE(PolylineData(i) % Vals)
          DEALLOCATE(PolylineData(i) % Prev, &
              PolylineData(i) % Next)
        END IF
      END DO
    END IF


    ! This assumes constant levelset convection. Mainly for testing.
    Vx = ListGetCReal( Solver % Values,'Levelset Velocity 1',Found )
    IF(Found) THEN
      IF(ABS(Vx) > EPSILON(Vx)) THEN
        MovingLevelset = .TRUE.
        x = x + Vx * dt
      END IF
    END IF
      
    Vy = ListGetCReal( Solver % Values,'Levelset Velocity 2',Found )
    IF(Found) THEN
      IF(ABS(Vy) > EPSILON(Vy)) THEN
        MovingLevelset = .TRUE.
        y = y + Vy * dt
      END IF
    END IF

    
    ! This assumes constant calving speed. Mainly for testing.
    VPhi = ListGetCReal( Solver % Values,'Levelset Speed',Found )
    IF(Found) THEN
      IF(ABS(VPhi) > EPSILON(VPhi)) THEN
        PhiVar1D % Values = PhiVar1D % Values + VPhi * dt 
        MovingLevelset = .TRUE.
        NonZero = .TRUE.
      END IF
    END IF
    
    ! Position dependent levelset velocity & calving speed.
    str = ListGetString( Solver % Values,'Levelset Velocity Variable',Found )
    IF(Found) THEN
      pVar => VariableGet( IsoMesh % Variables,TRIM(str)//' 1',UnfoundFatal=.TRUE.) 
      x = x + pVar % Values * dt
      pVar => VariableGet( IsoMesh % Variables,TRIM(str)//' 2',UnfoundFatal=.TRUE.) 
      y = y + pVar % Values * dt
      MovingLevelset = .TRUE.      
    END IF
      
    str = ListGetString( Solver % Values,'Levelset Speed Variable',Found )
    IF(Found) THEN
      VPhi = ListGetCReal( Solver % Values,'Levelset Speed Multiplier',Found )
      IF(.NOT. Found) VPhi = 1.0_dp
      pVar => VariableGet( IsoMesh % Variables,TRIM(str),UnfoundFatal=.TRUE.) 
      !PRINT *,'Levelset Speed range:',MINVAL(pVar % Values), MAXVAL(pVar % Values)
      PhiVar1D % Values = PhiVar1D % Values - pVar % Values * Vphi * dt
      !PRINT *,'Levelset range:',MINVAL(PhiVar1D % Values), MAXVAL(PhiVar1D % Values)
      MovingLevelset = .TRUE.
      Nonzero = .TRUE.
    END IF
              
    PhiMax = MAXVAL(ABS(PhiVar1D % Values)) 
    PhiMax = 1.01 * ( PhiMax + SQRT(Vx**2+Vy**2)*dt )

    IF(NormalMove .OR. NonZero) THEN
      CALL LevelsetNormalMove()
      NonZero = .FALSE.
    END IF

    IF(.NOT. ALLOCATED(PolylineData)) THEN
      ALLOCATE(PolylineData(ParEnv % PEs))
    END IF
    CALL PopulatePolyline()

    Trust = .FALSE.
    DO i=1, Mesh % NumberOfNodes
      j = PhiVar2D % Perm(i)
      IF(j==0 .AND. .NOT. NonZero) CYCLE
      IF(j==0) STOP
#if 0      
      val = PhiVar2D % Values(j)
      IF(val > BW ) THEN
        val = val - BW
      ELSE IF(val < -BW ) THEN
        val = val + BW
      ELSE        
        val = SignedDistance(i)
      END IF
#else
      val = SignedDistance(i, Trust(i))
#endif

      IF(MovingLevelset) THEN
        PhiVar2D % Values(j) = val
      END IF
    END DO

    IF(CheckLS) CALL CheckLSField()

    IF(NodeHistory) THEN
      PhiVar2D % Values(PhiVar2D % Perm) = PhiVar2D % Values(PhiVar2D % Perm) + NodeDiff
    END IF

    ! Deallocate data, next time this will be different.
    DO i=1,ParEnv % PEs 
      IF(PolylineData(i) % nLines > 0) THEN
        DEALLOCATE(PolylineData(i) % Vals)
        DEALLOCATE(PolylineData(i) % Prev, &
              PolylineData(i) % Next)
      END IF
    END DO
    
    Solver % Mesh % Next => IsoMesh
    
  CONTAINS

    !------------------------------------------------------------------------------
    !> Moves levelset coordinates in normal direction.
    !> This eliminates the need to compute distance relative to non-zero levelset.
    !------------------------------------------------------------------------------
    SUBROUTINE LevelsetNormalMove()
      !------------------------------------------------------------------------------      
      REAL(KIND=dp) :: x0,y0,x1,y1,ss,phi0,phi1,Normal(2),wei
      INTEGER :: i,j,k,n,i0,i1
      
      REAL(KIND=dp), ALLOCATABLE :: dx(:), dy(:), NodeWeight(:)
      
      n = SIZE(PhiVar1D % Values)
      ALLOCATE(dx(n),dy(n),NodeWeight(n))

      ! Allocate temporal space since we cannot make the computation in-place since the
      ! coordinates and levelset will be needed more than one times. 
      dx = 0.0_dp
      dy = 0.0_dp
      NodeWeight = 0.0_dp
            
      DO i=1,IsoMesh % NumberOfBulkElements        
        i0 = IsoMesh % Elements(i) % NodeIndexes(1)
        i1 = IsoMesh % Elements(i) % NodeIndexes(2)

        x0 = x(i0); y0 = y(i0)
        x1 = x(i1); y1 = y(i1)

        ss = (x0-x1)**2 + (y0-y1)**2
        IF(ss < EPSILON(ss) ) CYCLE        

        ss = SQRT(ss)
        Normal(1) = (y1-y0)
        Normal(2) = -(x1-x0)
        Normal = Normal / ss

        wei = 1.0_dp
        ! wei = ss

        phi0 = PhiVar1D % Values(i0)
        phi1 = PhiVar1D % Values(i1)

        dx(i0) = dx(i0) + Normal(1) * phi0          
        dy(i0) = dy(i0) + Normal(2) * phi0          
        NodeWeight(i0) = NodeWeight(i0) + wei

        dx(i1) = dx(i1) + Normal(1) * phi1          
        dy(i1) = dy(i1) + Normal(2) * phi1          
        NodeWeight(i1) = NodeWeight(i1) + wei
      END DO

      ! Finally move the coordinates to the direction of the averaged normal.
      WHERE(NodeWeight > EPSILON(wei))      
        x = x + dx / NodeWeight
        y = y + dy / NodeWeight 
      END WHERE
      PhiVar1D % Values = 0.0_dp
      
      DEALLOCATE(dx,dy,NodeWeight)

    END SUBROUTINE LevelsetNormalMove


! post processing for errors when line segs intersect either at nodes or arbitrary
! here we assume if previous elem was all +ve or all -ve then no intersection
! therefore any node that has a different sign to remaining nodes this a dir error
! modify to use bounding box for triangle-polyline intersects

    SUBROUTINE CheckLSField()

      ALLOCATE(DoesElemIntersect(Mesh % NumberOfBulkElements))
      CALL GetElemIntersectLogical(DoesElemIntersect)

      DO i=1, Mesh % NumberOfNodes
        IF(Trust(i)) CYCLE

        FOund = .FALSE.
        DO j=1, Mesh % NumberOfBulkElements
          NodeIndexes => Mesh % Elements(j) % NodeIndexes
          IF(.NOT. ANY(NodeIndexes == i)) CYCLE
          IF(DoesElemIntersect(j)) CYCLE
          Found=.TRUE.
          EXIT
        END DO

        IF(.NOT. Found) Trust(i) = .TRUE. ! we have no choice but to trust
      END DO

      counter = 0

      ! allocate boundary nodes

      IF(ParEnv % PEs > 1) THEN
        NNeighbours = COUNT(ParEnv % IsNeighbour)
        Neighbours = PACK( (/ (i,i=1,ParEnv % PEs) /), ParEnv % IsNeighbour)

        ALLOCATE(nSend(NNeighbours), nRecv(NNeighbours))
        MyPE = ParEnv % MyPE + 1
        ALLOCATE(LocalPerm( MAXVAL(Mesh % ParallelInfo % GlobalDOFs)))

        !local perm gdof>local
        DO i=1, Mesh % NumberOfNodes
          LocalPerm(Mesh % ParallelInfo % GlobalDOFs(i)) = i
        END DO
      END IF

      ! all reduce trust here

      DO WHILE(.TRUE.)

        IF(ALL(Trust) .AND. ParEnv % PEs == 1) EXIT

        DO i=1, Mesh % NumberOfBulkElements
          IF(DoesElemIntersect(i)) CYCLE

          NodeIndexes => Mesh % Elements(i) % NodeIndexes
          !loop elems cycle if intersect

          IF(COUNT(.NOT. Trust(NodeIndexes)) /= 1) CYCLE ! can only correct one at a time

          ! all same signs so we must trust all nodes
          ! only elems with some trusted nodes make it this for so ok
          IF(ALL(PhiVar2D % Values(PhiVar2D % Perm(NodeIndexes)) > 0.0_dp) .OR. &
            ALL(PhiVar2D % Values(PhiVar2D % Perm(NodeIndexes)) < 0.0_dp)) THEN
              Trust(NodeIndexes) = .TRUE.
          END IF

          Positive = .FALSE.
          IF(COUNT(PhiVar2D % Values(PhiVar2D % Perm(NodeIndexes)) > 0.0_dp) > 1) Positive =.TRUE.
          DO l=1, Mesh % Elements(i) % TYPE % NumberOfNodes
            IF(Positive .AND. PhiVar2D % Values(PhiVar2D % Perm(NodeIndexes(l))) < 0 .AND. &
                    .NOT. Trust(NodeIndexes(l))) THEN
              PhiVar2D % Values(PhiVar2D % Perm(NodeIndexes(l))) = &
                    -1*PhiVar2D % Values(PhiVar2D % Perm(NodeIndexes(l)))
              Trust(NodeIndexes(l)) = .TRUE.
            END IF
            IF(.NOT. Positive .AND. PhiVar2D % Values(PhiVar2D % Perm(NodeIndexes(l))) > 0 .AND. &
                    .NOT. Trust(NodeIndexes(l))) THEN
              PhiVar2D % Values(PhiVar2D % Perm(NodeIndexes(l))) = &
                    -1*PhiVar2D % Values(PhiVar2D % Perm(NodeIndexes(l)))
              Trust(NodeIndexes(l)) = .TRUE.
            END IF
          END DO
        END DO


        IF(ParEnv % PEs > 1) THEN
          !parcomm
          BLOCK
            INTEGER, ALLOCATABLE :: SendIndices(:,:),PTrustGDOFs(:)
            REAL(KIND=dp), ALLOCATABLE :: Pvals(:)
            INTEGER :: comm, ierr, status(MPI_STATUS_SIZE),Phase,cnt
            LOGICAL :: Finish

            comm = Solver % Matrix % Comm

            ! all trust reduce here
            Finish = ALL(Trust)
            CALL MPI_ALLREDUCE(MPI_IN_PLACE, Finish, 1, MPI_LOGICAL, MPI_LAND, comm, ierr)

            IF(Finish) EXIT

            CALL MPI_BARRIER(comm, ierr)

            DO Phase=0,1
              nSend = 0
              DO i=1, Mesh % NumberOfNodes
                IF(.NOT. Trust(i)) CYCLE
                DO j=1, NNeighbours
                  IF(ANY(Mesh % ParallelInfo % NeighbourList(i) % Neighbours == Neighbours(j)-1)) THEN
                    nSend(j) = nSend(j) + 1
                    IF(Phase == 1) SendIndices(j, nSend(j)) = i
                  END IF
                END DO
              END DO

              IF(Phase == 0) ALLOCATE(SendIndices(NNeighbours, MAXVAL(nSend)))
            END DO

            DO i=1,NNeighbours
              j = Neighbours(i)

              CALL MPI_BSEND( nSend(i), 1, MPI_INTEGER,j-1, &
                1001, comm, ierr )

              IF(nSend(i) == 0) CYCLE

              CALL MPI_BSEND( Mesh % ParallelInfo % GlobalDOFS(SendIndices(i, : nSend(i))), nSend(i), MPI_INTEGER,j-1, &
                1002, comm, ierr )
              CALL MPI_BSEND( PhiVar2D % Values(PhiVar2D % Perm(SendIndices(i, : nSend(i)))), nSend(i), MPI_DOUBLE_PRECISION,j-1, &
                1003, comm, ierr )

            END DO

            DO i=1,NNeighbours
              j = Neighbours(i)

              CALL MPI_RECV( nRecv(i), 1, MPI_INTEGER,j-1, &
                1001, comm, status, ierr )
            END DO

            ALLOCATE(PTrustGDOFs(SUM(nRecv)), Pvals(SUM(nRecv)))

            cnt = 0
            DO i=1,NNeighbours
              j = Neighbours(i)
              IF(nRecv(i) == 0) CYCLE

              CALL MPI_RECV( PTrustGDOFs(cnt+1:cnt+nRecv(i)), nRecv(i), MPI_INTEGER,j-1, &
                1002, comm, status, ierr )
              CALL MPI_RECV( Pvals(cnt+1:cnt+nRecv(i)), nRecv(i), MPI_DOUBLE_PRECISION,j-1, &
                1003, comm, status, ierr )
              cnt = cnt + nRecv(i)
            END DO


            CALL MPI_BARRIER( comm, ierr )

            Trust(LocalPerm(PTrustGDOFs)) = .TRUE.
            PhiVar2D % Values(PhiVar2D % Perm(LocalPerm(PTrustGDOFs))) = Pvals

            DEALLOCATE(PTrustGDOFs,Pvals,SendIndices)

          END BLOCK
        END IF

        counter = counter + 1
        IF(counter > 10 ) CALL FATAL('CheckLSField','Stuck in loop')
      END DO

    END SUBROUTINE CheckLSField

    !strategy loop through all and mark elem/polyline intersections
    ! use dir to get whether inside
    ! remove all until next intersection
    SUBROUTINE GetElemIntersectLogical(DoesElemIntersect)
      !------------------------------------------------------------------------------
      !TYPE(Mesh_t), POINTER :: Mesh
      !TYPE(Variable_t), POINTER :: PhiVar2D
      !TYPE(PolylineData_t), POINTER :: PolylineData
      LOGICAL :: DoesElemIntersect(:)
      !------------------------------------------------------------------------------
      INTEGER :: i,j,k,nLines
      INTEGER, POINTER :: NodeIndexes(:)
      REAL(KIND=dp) :: xmin_tri,xmax_tri,ymin_tri,ymax_tri,xmin_seg,xmax_seg,ymin_seg,ymax_seg,&
        a1(2),a2(2),intersect_p(2),x0,x1,y0,y1
      LOGICAL :: intersect

      DoesElemIntersect = .FALSE.
      DO i=1,Mesh % NumberOfBulkElements
        NodeIndexes => Mesh % Elements(i) % NodeIndexes

        xmin_tri = MINVAL(Mesh % Nodes % x(NodeIndexes))
        xmax_tri = MAXVAL(Mesh % Nodes % x(NodeIndexes))
        ymin_tri = MINVAL(Mesh % Nodes % y(NodeIndexes))
        ymax_tri = MAXVAL(Mesh % Nodes % y(NodeIndexes))

        IF(ALL(PhiVar2D % Values(PhiVar2D % Perm(NodeIndexes)) > 0.0_dp)) CYCLE
        IF(ALL(PhiVar2D % Values(PhiVar2D % Perm(NodeIndexes)) < 0.0_dp)) CYCLE

        ! can speed this up using local polylines rather than global
        Intersect = .FALSE.
        DO k=1,ParEnv % PEs
          DO j=1, PolylineData(k) % nLines

            x0 = PolylineData(k) % Vals(j,1)
            x1 = PolylineData(k) % Vals(j,2)
            y0 = PolylineData(k) % Vals(j,3)
            y1 = PolylineData(k) % Vals(j,4)

            xmin_seg = MIN(x0, x1)
            xmax_seg = MAX(x0, x1)
            ymin_seg = MIN(y0, y1)
            ymax_seg = MAX(y0, y1)

            ! Quick rejection
            IF(xmax_seg < xmin_tri .OR. xmax_tri < xmin_seg) CYCLE
            IF(ymax_seg < ymin_tri .OR. ymax_tri < ymin_seg) CYCLE

            ! do we actually intersect?
            a1(1) = Mesh % Nodes % x(NodeIndexes(1))
            a1(2) = Mesh % Nodes % y(NodeIndexes(1))
            a2(1) = Mesh % Nodes % x(NodeIndexes(2))
            a2(2) = Mesh % Nodes % y(NodeIndexes(2))
            CALL LineSegmentsIntersect(a1,a2,(/x0,y0/),(/x1,y1/),intersect_p,intersect)

            IF(.NOT. intersect) THEN !1,3
              a2(1) = Mesh % Nodes % x(NodeIndexes(3))
              a2(2) = Mesh % Nodes % y(NodeIndexes(3))
              CALL LineSegmentsIntersect(a1,a2,(/x0,y0/),(/x1,y1/),intersect_p,intersect)
            END IF

            IF(.NOT. intersect) THEN !2,3
              a1(1) = Mesh % Nodes % x(NodeIndexes(2))
              a1(2) = Mesh % Nodes % y(NodeIndexes(2))
              CALL LineSegmentsIntersect(a1,a2,(/x0,y0/),(/x1,y1/),intersect_p,intersect)
            END IF

            IF(Intersect) THEN
              DoesElemIntersect(i) = .TRUE.
              EXIT
            END IF

          END DO
          IF(Intersect) EXIT
        END DO
      END DO


    END SUBROUTINE GetElemIntersectLogical

    !------------------------------------------------------------------------------
    !> Computes the signed distance to zero levelset. 
    !------------------------------------------------------------------------------
    SUBROUTINE PopulatePolyline()
      !------------------------------------------------------------------------------      
      REAL(KIND=dp) :: x0,y0,x1,y1,ss,TotLineLen,phi0,phi1,&
        x2,x3,y2,y3,xmin0,xmax0,ymin0,ymax0,xmin1,xmax1,ymin1,ymax1,intersect_p(2)
      INTEGER :: i,j,k,n,m,i0,i1,nCol,dofs,k2,m2,i2,i3,NNeighbours,neighbour,p0,p1
      INTEGER, ALLOCATABLE :: PL_local(:,:),NodeToElem(:,:),NodeElemCount(:),RecvElems(:,:),&
        rdisps(:),RecvFrom(:),Neighbours(:)
      TYPE(Variable_t), POINTER :: Var1D
      INTEGER :: iVar,MyPe,PEs,Phase,nLines,Next,counter,NInter,SharedCount,NextP,nMax,NPInter,Nexts(2)
      LOGICAL, ALLOCATABLE :: Skip(:),RemoveLines(:,:)
      LOGICAL :: intersect
      !------------------------------------------------------------------------------

      nCol = 6

      nVar = 0
      iAvoid = 0
      iSolver = 0
      TotLineLen = 0.0_dp
      
      DO k = 1,100    
        str = ListGetString( Solver % Values,'isoline variable '//I2S(k), Found )
        IF(.NOT. Found ) EXIT            

        ! The levelset is really computed, do not interpolate it. 
        IF(str == CutStr) iAvoid = k 
        IF(str == Solver % Variable % Name) iSolver = k

        Var1D => VariableGet( IsoMesh % Variables, str, ThisOnly = .TRUE. )
        IF(.NOT. ASSOCIATED(Var1D)) EXIT
        nVar = k
        nCol = nCol + 2 * Var1D % Dofs
        
        IF(InfoActive(25)) THEN
          PRINT *,'Mesh1D Range for '//TRIM(Var1D % Name)//': ',&
              MINVAL(Var1D % Values), MAXVAL(Var1D % Values) 
        END IF
      END DO

      m = Isomesh % NumberOfBulkElements
      MyPe = ParEnv % MyPe + 1
      PEs = ParEnv % PEs

      IF(PEs > 1) THEN
        BLOCK
          INTEGER, ALLOCATABLE :: SendCount(:),SendElems(:,:)
          INTEGER :: comm, ierr, status(MPI_STATUS_SIZE)

          comm = Solver % Matrix % Comm

          ALLOCATE(SendCount(PEs))

          ! if a node is at an intersect can be shared with 3+ parts
          SendCount = 0
          DO k=1, PEs
            IF(k==MyPE) CYCLE
            DO i=1, IsoMesh % NumberOfBulkElements
              i0 = IsoMesh % Elements(i) % NodeIndexes(1)
              i1 = IsoMesh % Elements(i) % NodeIndexes(2)
              IF(.NOT. (IsoMesh % ParallelInfo % GInterface(i0) .OR. &
                IsoMesh % ParallelInfo % GInterface(i1)) ) CYCLE
              IF(IsoMesh % ParallelInfo % GInterface(i0)) THEN
                IF(ANY(IsoMesh % ParallelInfo % NeighbourList(i0) % Neighbours == k-1)) THEN
                  SendCount(k) = SendCount(k) + 1
                END IF
              END IF
              IF(IsoMesh % ParallelInfo % GInterface(i1)) THEN
                IF(ANY(IsoMesh % ParallelInfo % NeighbourList(i1) % Neighbours == k-1)) THEN
                  SendCount(k) = SendCount(k) + 1
                END IF
              END IF
            END DO
          END DO

          ALLOCATE(SendElems(SUM(SendCount),2))

          counter = 0
          DO k=1, PEs
            IF(k==MyPE) CYCLE
            DO i=1, IsoMesh % NumberOfBulkElements
              i0 = IsoMesh % Elements(i) % NodeIndexes(1)
              i1 = IsoMesh % Elements(i) % NodeIndexes(2)
              IF(.NOT. (IsoMesh % ParallelInfo % GInterface(i0) .OR. &
                IsoMesh % ParallelInfo % GInterface(i1)) ) CYCLE
              IF(IsoMesh % ParallelInfo % GInterface(i0)) THEN
                IF(ANY(IsoMesh % ParallelInfo % NeighbourList(i0) % Neighbours == k-1)) THEN
                  counter = counter + 1
                  SendElems(counter,1) = i
                  SendElems(counter,2) = IsoMesh % ParallelInfo % GlobalDOFS(i0)
                END IF
              END IF
              IF(IsoMesh % ParallelInfo % GInterface(i1)) THEN
                IF(ANY(IsoMesh % ParallelInfo % NeighbourList(i1) % Neighbours == k-1)) THEN
                  counter = counter + 1
                  SendElems(counter,1) = i
                  SendElems(counter,2) = IsoMesh % ParallelInfo % GlobalDOFS(i1)
                END IF
              END IF
            END DO
          END DO

          ALLOCATE(rdisps(PEs))

          ! send count should equal recvcount
          !CALL MPI_ALLTOALL(SendCount, 1, MPI_INTEGER, RecvCount, 1, MPI_INTEGER, comm, ierr)

          rdisps(1) = 0
          DO i=2, PEs
            rdisps(i) = rdisps(i-1) + SendCount(i-1)
          END DO

          SharedCount = SUM(SendCount)
          ALLOCATE(RecvElems(SharedCount,2))

          ! these can be all gather?!
          CALL MPI_Alltoallv(SendElems(:,1), SendCount, rdisps, MPI_INTEGER, &
                            RecvElems(:,1), SendCount, rdisps, MPI_INTEGER, &
                            comm, ierr)
          CALL MPI_Alltoallv(SendElems(:,2), SendCount, rdisps, MPI_INTEGER, &
                            RecvElems(:,2), SendCount, rdisps, MPI_INTEGER, &
                            comm, ierr)

          ALLOCATE(RecvFrom(SharedCount))
          counter = 0
          DO i=1,PEs
            RecvFrom(counter+1:counter+SendCount(i)) = i
            counter = counter + SendCount(i)
          END DO

        END BLOCK
      END IF

      ! We may find use for bounding boxes later on. 

      PolylineData(MyPe) % IsoLineBB(1) = MINVAL(IsoMesh % Nodes % x)
      PolylineData(MyPe) % IsoLineBB(2) = MAXVAL(IsoMesh % Nodes % x)
      PolylineData(MyPe) % IsoLineBB(3) = MINVAL(IsoMesh % Nodes % y)
      PolylineData(MyPe) % IsoLineBB(4) = MAXVAL(IsoMesh % Nodes % y)
#if 0
      PolylineData(MyPe) % MeshBB(1) = MINVAL(Mesh % Nodes % x)
      PolylineData(MyPe) % MeshBB(2) = MAXVAL(Mesh % Nodes % x)
      PolylineData(MyPe) % MeshBB(3) = MINVAL(Mesh % Nodes % y)
      PolylineData(MyPe) % MeshBB(4) = MAXVAL(Mesh % Nodes % y)
#endif

      ALLOCATE(Skip(IsoMesh % NumberOfBulkElements), &
        NodeToElem(IsoMesh % NumberOfNodes,2), NodeElemCount(IsoMesh % NumberOfNodes))

      IF(.NOT. ALLOCATED(Intersects)) &
        ALLOCATE(Intersects(PEs))
      CALL InitialiseIntersects(10)

      Skip = .FALSE.
      NodeToElem = 0; NodeElemCount = 0
      NInter = 0
      DO Phase=0,1
        m = 0
        DO i=1,IsoMesh % NumberOfBulkElements        
          i0 = IsoMesh % Elements(i) % NodeIndexes(1)
          i1 = IsoMesh % Elements(i) % NodeIndexes(2)

          x0 = x(i0); y0 = y(i0)
          x1 = x(i1); y1 = y(i1)
          
          ss = (x0-x1)**2 + (y0-y1)**2
          
          ! This is too short for anything useful...
          ! Particularly difficult it is to decide on left/right if the segment is a stub.
          IF(ss < EPSILON(ss) ) THEN
            Skip(i) = .TRUE.
            CYCLE
          END IF
                    
          m = m+1

          ! assumption here is there are no divides or intersections at nodes
          ! should be fine sine this is only used for arbitray intersections
          ! shared node intersections dealt with later
          IF(Phase == 0) THEN
            NodeElemCount(i0) = NodeElemCount(i0) + 1
            NodeElemCount(i1) = NodeElemCount(i1) + 1

            IF(PEs > 1) THEN
            ! can't use SIZE(Mesh % ParallelInfo % NeighbourList(i0) % Neighbours)-1
            ! as this includes all neighbours to the bulk element from Mesh not the ones to CutFEM
            ! so we have to check against shared globaldofs
              IF(IsoMesh % ParallelInfo % GInterface(i0) .OR. &
                  IsoMesh % ParallelInfo % GInterface(i1)) THEN
                DO j=1, SharedCount
                  IF(IsoMesh % ParallelInfo % GlobalDOFs(i0) == RecvElems(j,2)) THEN
                    NodeElemCount(i0) = NodeElemCount(i0) + 1
                  END IF
                  IF(IsoMesh % ParallelInfo % GlobalDOFs(i1) == RecvElems(j,2)) THEN
                    NodeElemCount(i1) = NodeElemCount(i1) + 1
                  END IF
                END DO
              END IF
            END IF

            IF(NodeToElem(i0,1) == 0) THEN
              NodeToElem(i0,1) = m
            ELSE
              NodeToElem(i0,2) = m
            END IF

            IF(NodeToElem(i1,1) == 0) THEN
              NodeToElem(i1,1) = m
            ELSE
              NodeToElem(i1,2) = m
            END IF
          END IF


          IF(Phase==0) CYCLE

          !assign neigbours
          IF(NodeToElem(i0,1) == i) THEN
            PolylineData(MyPe) % Prev(m,1) = NodeToElem(i0,2)
          ELSE
            PolylineData(MyPe) % Prev(m,1) = NodeToElem(i0,1)
          END IF

          IF(NodeToElem(i1,1) == i) THEN
            PolylineData(MyPe) % Next(m,1) = NodeToElem(i1,2)
          ELSE
            PolylineData(MyPe) % Next(m,1) = NodeToElem(i1,1)
          END IF

          IF(PEs > 1) THEN
            IF(IsoMesh % ParallelInfo % GInterface(i0) .OR. &
                  IsoMesh % ParallelInfo % GInterface(i1)) THEN
              DO j=1, SharedCount
                IF(IsoMesh % ParallelInfo % GlobalDOFs(i0) == RecvElems(j,2)) THEN
                  PolylineData(MyPe) % Prev(m,1) = RecvElems(j,1)
                  PolylineData(MyPe) % Prev(m,2) = RecvFrom(j)
                END IF
                IF(IsoMesh % ParallelInfo % GlobalDOFs(i1) == RecvElems(j,2)) THEN
                  PolylineData(MyPe) % Next(m,1) = RecvElems(j,1)
                  PolylineData(MyPe) % Next(m,2) = RecvFrom(j)
                END IF
              END DO
            END IF
          END IF

          ! node on intersection so mark with negative
          IF(NodeElemCount(i0) > 2) PolylineData(MyPe) % Prev(m,1) = -1
          IF(NodeElemCount(i1) > 2) PolylineData(MyPe) % Next(m,1) = -1
          !IF(NodeElemCount(i0) > 2 .OR. NodeElemCount(i1) > 2) &
          !  PRINT*, ParEnv % MyPE, 'nodeelemcount', NodeElemCount(i0), NodeElemCount(i1), &
          !    'indices', i0, i1

          PL_local(m,1) = i0
          PL_local(m,2) = i1

          TotLineLen = TotLineLen + SQRT(ss)

          phi0 = PhiVar1D % Values(i0)
          phi1 = PhiVar1D % Values(i1)
          
          ! Coordinates for the polyline. 
          PolylineData(MyPe) % Vals(m,1) = x0
          PolylineData(MyPe) % Vals(m,2) = x1
          PolylineData(MyPe) % Vals(m,3) = y0
          PolylineData(MyPe) % Vals(m,4) = y1

          ! Levelset values for the polyline. 
          PolylineData(MyPe) % Vals(m,5) = phi0
          PolylineData(MyPe) % Vals(m,6) = phi1


          ! min/max for bounds check
          xmax0 = MAX(x0,x1); xmin0 = MIN(x0,x1)
          ymax0 = MAX(y0,y1); ymin0 = MIN(y0,y1)

          ! check and mark intersections
          m2 = m
          DO j=i+1, IsoMesh % NumberOfBulkElements
            IF(Skip(j)) CYCLE
            m2 = m2 + 1

            i2 = IsoMesh % Elements(j) % NodeIndexes(1)
            i3 = IsoMesh % Elements(j) % NodeIndexes(2)

            IF(i2==i0 .OR. i2==i1 .OR. i3==i0 .OR. i3==i1) CYCLE ! share node

            x2 = x(i2); y2 = y(i2)
            x3 = x(i3); y3 = y(i3)

            xmax1 = MAX(x2,x3); xmin1 = MIN(x2,x3)
            ymax1 = MAX(y2,y3); ymin1 = MIN(y2,y3)

            ! bounds check can we intersect
            IF(xmax1 < xmin0 .OR. xmax0 < xmin1) CYCLE
            IF(ymax1 < ymin0 .OR. ymax0 < ymin1) CYCLE

            ! do we actually intersect?
            CALL LineSegmentsIntersect((/x2,y2/),(/x3,y3/),(/x0,y0/),(/x1,y1/),intersect_p,intersect)

            IF(Intersect) THEN

              Intersects(MyPE) % nIntersects = Intersects(MyPE) % nIntersects + 1
              NInter = Intersects(MyPE) % nIntersects

              IF(NInter > Intersects(MyPE) % Size) &
                  CALL DoubleIntersectsAllocation()

              Intersects(MyPE) % Elements(NInter,1) = m
              Intersects(MyPE) % Elements(NInter,2) = m2
              Intersects(MyPE) % Elements(NInter,3) = MyPE
              Intersects(MyPE) % Elements(NInter,4) = MyPE

              Intersects(MyPE) % Intersection(NInter,1) = intersect_p(1)
              Intersects(MyPE) % Intersection(NInter,2) = intersect_p(2)

              Found = GetPointLevelset((/x0,y0,0.0_dp/), val)
              Intersects(MyPE) % Found(NInter,1) = Found
              IF(Found) Intersects(MyPE) % LS(NInter,1) = val

              Found = GetPointLevelset((/x1,y1,0.0_dp/), val)
              Intersects(MyPE) % Found(NInter,2) = Found
              IF(Found) Intersects(MyPE) % LS(NInter,2) = val

              Found = GetPointLevelset((/x2,y2,0.0_dp/), val)
              Intersects(MyPE) % Found(NInter,3) = Found
              IF(Found) Intersects(MyPE) % LS(NInter,3) = val

              Found = GetPointLevelset((/x3,y3,0.0_dp/), val)
              Intersects(MyPE) % Found(NInter,4) = Found
              IF(Found) Intersects(MyPE) % LS(NInter,4) = val

              EXIT ! assumption is that each line can only have one intersection
            END IF
          END DO

          j = 7

          DO k = 1,nVar
            IF(k==iAvoid) CYCLE

            str = ListGetString( Solver % Values,'isoline variable '//I2S(k), Found )
            Var1D => VariableGet( IsoMesh % Variables, str, ThisOnly = .TRUE. )

            dofs = Var1D % Dofs
            DO k2=1,dofs              
              PolylineData(MyPe) % Vals(m,j)   = Var1D % Values(dofs*(Var1D % Perm(i0)-1)+k2)
              PolylineData(MyPe) % Vals(m,j+1) = Var1D % Values(dofs*(Var1D % Perm(i1)-1)+k2)
              j = j+2
            END DO
          END DO
        END DO

        
        IF(Phase==0) THEN
          CALL Info('LevelsetUpdate','Allocating PolylineData of size '//I2S(m)//' x '//I2S(nCol),Level=8)
          PolylineData(MyPe) % nLines = m
          PolylineData(MyPe) % nNodes = Mesh % NumberOfNodes
          ALLOCATE(PolylineData(MyPe) % Vals(m,nCol), PolylineData(MyPe) % Next(m,2),&
              PolylineData(MyPE) % Prev(m,2), PL_local(m,2))
          PolylineData(MyPe) % Vals = 0.0_dp
          PolylineData(MyPe) % Next(:,1) = 0
          PolylineData(MyPe) % Prev(:,1) = 0
          PolylineData(MyPe) % Prev(:,2) = MyPE
          PolylineData(MyPe) % Next(:,2) = MyPE
        END IF
        
      END DO

      !need to get share local intersects
      IF(PEs > 1 ) THEN        
        BLOCK
          INTEGER, ALLOCATABLE :: nPar(:)
          INTEGER :: comm, ierr, status(MPI_STATUS_SIZE)
          
          ALLOCATE(nPar(PEs))
          comm = Solver % Matrix % Comm

          nPar = 0
          nPar(MyPe) = PolylineData(MyPe) % nLines
          CALL MPI_ALLREDUCE(MPI_IN_PLACE, nPar, PEs, MPI_INTEGER, MPI_MAX, comm, ierr)
          DO i=1,PEs
            PolylineData(i) % nLines = nPar(i)
          END DO
          
          nPar = 0
          nPar(MyPe) = Mesh % NumberOfNodes
          CALL MPI_ALLREDUCE(MPI_IN_PLACE, nPar, PEs, MPI_INTEGER, MPI_MAX, comm, ierr)
          DO i=1,PEs
            PolylineData(i) % nNodes = nPar(i)
          END DO
          CALL MPI_BARRIER( comm, ierr )
          
          IF( PolylineData(MyPe) % nNodes > 1) THEN
            DO i=1,PEs            
              IF(i==MyPe) CYCLE
              m = PolylineData(i) % nLines
              IF(m>0) ALLOCATE(PolylineData(i) % Vals(m,nCol), &
                  PolylineData(i) % Next(m,2), PolylineData(i) % Prev(m,2))
            END DO
          END IF

          DO i=1,PEs
            IF(i==MyPe) CYCLE              
            IF(PolylineData(MyPe) % nLines == 0 .OR. PolylineData(i) % nNodes == 0 ) CYCLE

            ! Sent data from partition MyPe to i           
            k = PolylineData(MyPe) % nLines * nCol
            CALL MPI_BSEND( PolylineData(MyPe) % Vals, k, MPI_DOUBLE_PRECISION,i-1, &
                1001, comm, ierr )
            k = PolylineData(MyPe) % nLines * 2
            CALL MPI_BSEND( PolylineData(MyPe) % Next, k, MPI_INTEGER,i-1, &
                1002, comm, ierr )
            CALL MPI_BSEND( PolylineData(MyPe) % Prev, k, MPI_INTEGER,i-1, &
                1003, comm, ierr )
          END DO
            
          DO i=1,PEs
            IF(i==MyPe) CYCLE              
            IF(PolylineData(i) % nLines == 0 .OR. PolylineData(MyPe) % nNodes == 0 ) CYCLE
            
            ! Receive data from partition i to MyPe
            k = PolylineData(i) % nLines * nCol
            CALL MPI_RECV( PolylineData(i) % Vals, k, MPI_DOUBLE_PRECISION,i-1, &
                1001, comm, status, ierr )
            k = PolylineData(i) % nLines * 2
            CALL MPI_RECV( PolylineData(i) % Next, k, MPI_INTEGER,i-1, &
                1002, comm, status, ierr )
            CALL MPI_RECV( PolylineData(i) % Prev, k, MPI_INTEGER,i-1, &
                1003, comm, status, ierr )
          END DO

          CALL MPI_BARRIER( comm, ierr )

          !comm local intersects
          nPar = 0
          nPar(MyPE) = NInter
          CALL MPI_ALLREDUCE(MPI_IN_PLACE, nPar, PEs, MPI_INTEGER, MPI_MAX, comm, ierr)

          k = SUM( PolylineData(1:PEs) % nLines ) 
          CALL Info('LevelSetUpdate','Number of line segments in parallel system: '//I2S(k),Level=7)
          
          CALL MPI_ALLREDUCE(MPI_IN_PLACE, TotLineLen, PEs, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)                
        END BLOCK
      END IF

      ! find local-parallel intersects
      IF(PEs > 1) THEN
        ! could instead search through isomesh and get neighbours?
        NNeighbours = COUNT(ParEnv % IsNeighbour)
        Neighbours = PACK( (/ (i,i=1,PEs) /), ParEnv % IsNeighbour)

        ! only local to parallel intersections
        ! then store if smaller part store intersect info
        ! store intersect direction of local line
        nLines = PolylineData(MyPE) % nLines

        NPInter = 0

        nLines = PolylineData(MyPE) % nLines

        DO i=1,nLines
          x0 = PolylineData(MyPE) % Vals(i,1)
          x1 = PolylineData(MyPE) % Vals(i,2)
          y0 = PolylineData(MyPE) % Vals(i,3)
          y1 = PolylineData(MyPE) % Vals(i,4)

          xmax0 = MAX(x0,x1); xmin0 = MIN(x0,x1)
          ymax0 = MAX(y0,y1); ymin0 = MIN(y0,y1)

          Intersect = .FALSE.
          DO neighbour = 1, NNeighbours

            IF(MyPE > Neighbours(neighbour)) CYCLE ! only confirm on lowest part

            k = neighbours(neighbour)

            ! part not in BB
            IF(PolylineData(MyPE) % IsoLineBB(MyPE) < xmin0 .OR. &
                xmax0 < PolylineData(k) % IsoLineBB(1)) CYCLE
            IF(PolylineData(MyPe) % IsoLineBB(MyPE) < ymin0 .OR. &
                ymax0 < PolylineData(k) % IsoLineBB(3)) CYCLE

            DO j=1, PolylineData(k) % NLines

              x2 = PolylineData(k) % Vals(j,1)
              x3 = PolylineData(k) % Vals(j,2)
              y2 = PolylineData(k) % Vals(j,3)
              y3 = PolylineData(k) % Vals(j,4)

              !need to discount if either node is the same
              IF(ABS(x0-x2) > AEPS .OR. ABS(x0-x3) <  AEPS) CYCLE
              IF(ABS(x1-x2) > AEPS .OR. ABS(x1-x3) <  AEPS) CYCLE
              IF(ABS(y0-y2) > AEPS .OR. ABS(y0-y3) <  AEPS) CYCLE
              IF(ABS(y1-y2) > AEPS .OR. ABS(y1-y3) <  AEPS) CYCLE

              xmax1 = MAX(x2,x3); xmin1 = MIN(x2,x3)
              ymax1 = MAX(y2,y3); ymin1 = MIN(y2,y3)

              ! bounds check can we intersect
              IF(xmax1 < xmin0 .OR. xmax0 < xmin1) CYCLE
              IF(ymax1 < ymin0 .OR. ymax0 < ymin1) CYCLE

              ! do we actually intersect?
              CALL LineSegmentsIntersect((/x2,y2/),(/x3,y3/),(/x0,y0/),(/x1,y1/),intersect_p,intersect)

              IF(Intersect) THEN

                Intersects(MyPE) % nIntersects = Intersects(MyPE) % nIntersects + 1
                NInter = Intersects(MyPE) % nIntersects

                IF(NInter > Intersects(MyPE) % Size) &
                  CALL DoubleIntersectsAllocation()

                Intersects(MyPE) % Elements(NInter,1) = i
                Intersects(MyPE) % Elements(NInter,2) = j
                Intersects(MyPE) % Elements(NInter,3) = MyPE
                Intersects(MyPE) % Elements(NInter,4) = k

                Intersects(MyPE) % Intersection(NInter,1) = intersect_p(1)
                Intersects(MyPE) % Intersection(NInter,2) = intersect_p(2)

                Found = GetPointLevelset((/x0,y0,0.0_dp/), val)
                Intersects(MyPE) % Found(NInter,1) = Found
                IF(Found) Intersects(MyPE) % LS(NInter,1) = val

                Found = GetPointLevelset((/x1,y1,0.0_dp/), val)
                Intersects(MyPE) % Found(NInter,2) = Found
                IF(Found) Intersects(MyPE) % LS(NInter,2) = val

                Found = GetPointLevelset((/x2,y2,0.0_dp/), val)
                Intersects(MyPE) % Found(NInter,3) = Found
                IF(Found) Intersects(MyPE) % LS(NInter,3) = val

                Found = GetPointLevelset((/x3,y3,0.0_dp/), val)
                Intersects(MyPE) % Found(NInter,4) = Found
                IF(Found) Intersects(MyPE) % LS(NInter,4) = val

                EXIT
              END IF
            END DO

            IF(Intersect) EXIT
          END DO
        END DO
      END IF ! PES > 1

      CALL ReduceLocalIntersectMemory()

      ! this block breaks complier
      IF(PEs > 1) THEN
        !parallel comm
        ! share intersects_t
        BLOCK
          INTEGER, ALLOCATABLE :: nPar(:),UpdateIndices(:,:),AllUpdateIndices(:)
          REAL(KIND=dp), ALLOCATABLE :: UpdateVals(:),AllUpdateVals(:)
          INTEGER :: comm, ierr, status(MPI_STATUS_SIZE), Updates

          ALLOCATE(nPar(PEs))
          comm = Solver % Matrix % Comm

          nPar = 0
          nPar(MyPe) = Intersects(MyPe) % nIntersects
          CALL MPI_ALLREDUCE(MPI_IN_PLACE, nPar, PEs, MPI_INTEGER, MPI_MAX, comm, ierr)
          DO i=1,PEs
            Intersects(i) % nIntersects = nPar(i)
          END DO

          IF( PolylineData(MyPe) % nNodes > 1) THEN
            DO i=1,PEs
              IF(i==MyPe) CYCLE
              m = Intersects(i) % nIntersects
              IF(m>0) ALLOCATE(Intersects(i) % Intersection(m,2), &
                  Intersects(i) % Elements(m,4), Intersects(i) % Found(m,4), &
                  Intersects(i) % LS(m,4))
            END DO
          END IF

          DO i=1,PEs
            IF(i==MyPe) CYCLE
            IF(Intersects(MyPe) % nIntersects == 0 .OR. PolylineData(i) % nNodes == 0 ) CYCLE

            ! Sent data from partition MyPe to i
            k = Intersects(MyPe) % nIntersects * 2
            CALL MPI_BSEND( Intersects(MyPe) % Intersection, k, MPI_DOUBLE_PRECISION,i-1, &
                1004, comm, ierr )

            k = Intersects(MyPe) % nIntersects * 4
            CALL MPI_BSEND( Intersects(MyPe) % Elements, k, MPI_INTEGER,i-1, &
                1005, comm, ierr )
            CALL MPI_BSEND( Intersects(MyPe) % Found, k, MPI_LOGICAL,i-1, &
                1006, comm, ierr )
            CALL MPI_BSEND( Intersects(MyPe) % LS, k, MPI_DOUBLE_PRECISION,i-1, &
                1007, comm, ierr )
          END DO

          DO i=1,PEs
            IF(i==MyPe) CYCLE
            IF(Intersects(i) % nIntersects == 0 .OR. PolylineData(MyPe) % nNodes == 0 ) CYCLE

            ! Receive data from partition i to MyPe
            k = Intersects(i) % nIntersects * 2
            CALL MPI_RECV( Intersects(i) % Intersection, k, MPI_DOUBLE_PRECISION,i-1, &
                1004, comm, status, ierr )

            k = Intersects(i) % nIntersects * 4
            CALL MPI_RECV( Intersects(i) % Elements, k, MPI_INTEGER,i-1, &
                1005, comm, status, ierr )
            CALL MPI_RECV( Intersects(i) % Found, k, MPI_LOGICAL,i-1, &
                1006, comm, status, ierr )
            CALL MPI_RECV( Intersects(i) % LS, k, MPI_DOUBLE_PRECISION,i-1, &
                1007, comm, status, ierr )
          END DO

          CALL MPI_BARRIER( comm, ierr )


          !once parcomm over loop through unfound on neighbours
          nMax = SUM(Intersects(Neighbours) % nIntersects) ! should be low
          ALLOCATE(UpdateIndices(nMax,3),UpdateVals(nMax))
          Updates = 0
          DO neighbour = 1, NNeighbours
            k = neighbours(neighbour)

            nInter = Intersects(k) % nIntersects

            IF(nInter == 0) CYCLE

            DO i=1, nInter

              IF(ALL(Intersects(k) % Found(i,:))) CYCLE

              i0 = Intersects(k) % Elements(i,1)
              p0 = Intersects(k) % Elements(i,3)
              i1 = Intersects(k) % Elements(i,2)
              p1 = Intersects(k) % Elements(i,4)

              DO j=1,4
                IF(Intersects(k) % Found(i,j)) CYCLE

                IF(j > 2) THEN
                  i0=i1; p0=p1
                END IF

                IF(j==1 .OR. j==3) THEN
                  x0 = PolylineData(p0) % vals(i0,1)
                  y0 = PolylineData(p0) % vals(i0,3)
                ELSE
                  x0 = PolylineData(p0) % vals(i0,2)
                  y0 = PolylineData(p0) % vals(i0,4)
                END IF

                Found = GetPointLevelset((/x0,y0,0.0_dp/), val)

                IF(Found) THEN
                  Intersects(k) % Found(i,j) = Found
                  Intersects(k) % LS(i,j) = val

                  Updates = Updates + 1
                  UpdateIndices(Updates,1) = k
                  UpdateIndices(Updates,2) = i
                  UpdateIndices(Updates,3) = j
                  UpdateVals(Updates) = val

                  ! need to comm this?
                END IF
              END DO
            END DO
          END DO

          nPar = 0
          nPar(MyPE) = Updates
          CALL MPI_ALLREDUCE(MPI_IN_PLACE, nPar, PEs, MPI_INTEGER, MPI_MAX, comm, ierr)

          rdisps(1) = 0
          DO i=2,PEs
            rdisps(i) = rdisps(i-1) + nPar(i-1)
          END DO

          Updates = SUM(nPar)
          ALLOCATE(AllUpdateVals(Updates),AllUpdateIndices(Updates*3))
          CALL MPI_ALLGATHERV(UpdateVals, nPar(MyPE), MPI_DOUBLE_PRECISION, &
                AllUpdateVals, nPar, rdisps, MPI_DOUBLE_PRECISION, comm, ierr)

          CALL MPI_ALLGATHERV(UpdateIndices, nPar(MyPE)*3, MPI_INTEGER, &
                AllUpdateIndices, nPar*3, rdisps*3, MPI_INTEGER, comm, ierr)

          CALL MPI_BARRIER(comm, ierr)

          ! update values
          DO n=1,Updates
            k = AllUpdateIndices((n-1)*3 + 1)
            i = AllUpdateIndices((n-1)*3 + 2)
            j = AllUpdateIndices((n-1)*3 + 3)
            Intersects(k) % LS(i,j) = AllUpdateVals(n)
          END DO

        END BLOCK
      END IF ! parallel

      nMax = 0
      DO i=1,PEs
        nLines = PolylineData(i) % nLines
        ALLOCATE(PolylineData(i) % Intersect(nLines))
        PolylineData(i) % Intersect = .FALSE.
        IF(nMax < PolylineData(i) % nLines) nMax = PolylineData(i) % nLines
      END DO

      DO i=1,PEs
        nInter = Intersects(i) % nIntersects
        IF(nInter == 0) CYCLE
        DO j=1,Ninter
          i0 = Intersects(i) % Elements(j,1)
          p0 = Intersects(i) % Elements(j,3)
          i1 = Intersects(i) % Elements(j,2)
          p1 = Intersects(i) % Elements(j,4)

          PolylineData(p0) % Intersect(i0) = .TRUE.
          PolylineData(p1) % Intersect(i1) = .TRUE.
        END DO
      END DO

      ALLOCATE(RemoveLines(PEs, nMax))
      RemoveLines = .FALSE.
      ! find intersections
      DO k = 1, ParEnv % PEs
        nInter = Intersects(k) % nIntersects
        IF(nInter == 0) CYCLE

        DO i = 1, nInter

          DO phase=1,2

            IF(Phase == 1) THEN
              i0 = Intersects(k) % Elements(i,1)
              p0 = Intersects(k) % Elements(i,3)
              Positive = Intersects(k) % LS(i,1) > Intersects(k) % LS(i,2)
            ELSE
              i0 = Intersects(k) % Elements(i,2)
              p0 = Intersects(k) % Elements(i,4)
              Positive = Intersects(k) % LS(i,3) > Intersects(k) % LS(i,4)
            END IF

            IF(Positive) THEN
              next = PolylineData(p0) % Next(i0,1)
              nextp = PolylineData(p0) % Next(i0,2)
              counter = 0
              DO WHILE(.TRUE.)
                ! if part of an intersect? hashmap?
                IF(PolylineData(NextP) % Next(Next,1) > 0 .AND. .NOT. PolylineData(NextP) % Intersect(Next)) THEN
                  RemoveLines(NextP, Next) = .TRUE.
                  Nexts = PolylineData(NextP) % Next(Next,:)
                  !NextP = PolylineData(NextP) % Next(Next,2)
                  Next = Nexts(1)
                  NextP = Nexts(2)
                ELSE
                  EXIT
                END IF
                counter =counter+1
                IF(counter>10) CALL FATAL('CutFEM', 'stuck in loop')
              END DO

              ! adjust node
              PolylineData(p0) % Vals(i0,2) = Intersects(k) % Intersection(i,1)
              PolylineData(p0) % Vals(i0,4) = Intersects(k) % Intersection(i,2)

              ! adjust mesh for if we want to save
              ! this shared node could present on multiple parts?
              IF(p0==MyPe) THEN
                IsoMesh % Nodes % x(PL_local(i0,2)) = Intersects(k) % Intersection(i,1)
                IsoMesh % Nodes % y(PL_local(i0,2)) = Intersects(k) % Intersection(i,2)
              END IF

            ELSE ! negative

              next = PolylineData(p0) % Prev(i0,1)
              nextp = PolylineData(p0) % Prev(i0,2)
              !RemoveLines(NextP,Next) = .TRUE.
              DO WHILE(.TRUE.)
                IF(PolylineData(NextP) % Prev(Next,1) > 0 .AND. .NOT. PolylineData(NextP) % Intersect(Next)) THEN
                  RemoveLines(NextP,Next) = .TRUE.
                  Nexts = PolylineData(NextP) % Prev(Next,:)
                  ! this is changing as Next changes
                  !NextP = PolylineData(NextP) % Prev(Next,2)
                  Next = Nexts(1)
                  NextP = Nexts(2)
                ELSE
                  EXIT
                END IF
              END DO

              ! adjust node
              PolylineData(p0) % Vals(i0,1) = Intersects(k) % Intersection(i,1)
              PolylineData(p0) % Vals(i0,3) = Intersects(k) % Intersection(i,2)

              ! adjust mesh for if we want to save
              ! what about shared nodes?
              IF(p0==MyPe) THEN
                IsoMesh % Nodes % x(PL_local(i0,1)) = Intersects(k) % Intersection(i,1)
                IsoMesh % Nodes % y(PL_local(i0,1)) = Intersects(k) % Intersection(i,2)
              END IF

            END IF
          END DO
        END DO
      END DO

      CALL DeallocateMeshLines(IsoMesh, RemoveLines)
      CALL DeallocatePolyLines(RemoveLines, nCol)

      WRITE(Message,'(A,ES12.3)') 'Cutfem isoline length:',TotLineLen
      CALL Info('LevelSetUpdate',Message,Level=6)

      CALL ListAddConstReal(CurrentModel % Simulation,'res: cutfem isoline length',TotLineLen )
      
      IF(ALLOCATED(Intersects)) THEN
        DO i=1,PEs
          IF(ALLOCATED(Intersects(i) % Intersection)) &
            DEALLOCATE(Intersects(i) % Intersection)
          IF(ALLOCATED(Intersects(i) % Elements)) &
            DEALLOCATE(Intersects(i) % Elements)
          IF(ALLOCATED(Intersects(i) % Found)) &
            DEALLOCATE(Intersects(i) % Found)
          IF(ALLOCATED(Intersects(i) % LS)) &
            DEALLOCATE(Intersects(i) % LS)
        END DO
      END IF

      DO i=1, PEs
        IF(ALLOCATED(PolylineData(i) % Intersect)) &
          DEALLOCATE(PolylineData(i) % Intersect)
      END DO

      IF(InfoActive(25)) THEN
        CALL Info('LevelSetUpdate','Polyline interval for Isoline variables')

        j = SIZE(PolylineData(MyPe) % Vals(:,1))
        CALL VectorValuesRange(PolylineData(MyPe) % Vals(:,1),j,'x')
        CALL VectorValuesRange(PolylineData(MyPe) % Vals(:,3),j,'y')
        CALL VectorValuesRange(PolylineData(MyPe) % Vals(:,5),j,'phi')
        i = 7
        DO k = 1,nVar
          str = ListGetString( Solver % Values,'isoline variable '//I2S(k), Found )          
          CALL VectorValuesRange(PolylineData(MyPe) % Vals(:,i),j,TRIM(str))
          i = i+2
        END DO
      END IF
      
    END SUBROUTINE PopulatePolyline
    !------------------------------------------------------------------------------

    !initialise local intersects
    SUBROUTINE InitialiseIntersects(Size)
      INTEGER :: Size
      !-----------------------------
      INTEGER :: MyPE

      MyPE = ParEnv % MyPE + 1

      ALLOCATE(Intersects(MyPE) % Intersection(Size,2), &
              Intersects(MyPE) % Elements(Size,4), &
              Intersects(MyPE) % Found(Size,4), &
              Intersects(MyPE) % LS(Size,4))

      Intersects(MyPE) % nIntersects = 0
      Intersects(MyPE) % Size = Size
      Intersects(MyPE) % Found = .FALSE.
    END SUBROUTINE InitialiseIntersects

    ! double local intersect memmory on the fly
    SUBROUTINE DoubleIntersectsAllocation()
      !-----------------------------
      REAL(KIND=dp), ALLOCATABLE :: WorkReal(:,:)
      LOGICAL, ALLOCATABLE :: WorkLog(:,:)
      INTEGER, ALLOCATABLE :: WorkInt(:,:)
      INTEGER :: Size,Size2,MyPE

      MyPE = ParEnv % MyPE + 1
      Size = Intersects(MyPE) % Size
      Size2 = Size*2

      ALLOCATE(WorkReal(Size,4))
      WorkReal = Intersects(MyPE) % LS
      DEALLOCATE(Intersects(MyPE) % LS)
      ALLOCATE(Intersects(myPE) % LS(Size2,4))
      Intersects(MyPE) % LS(:Size,:) = WorkReal

      WorkReal(:,1:2) = Intersects(MyPE) % Intersection
      DEALLOCATE(Intersects(MyPE) % Intersection)
      ALLOCATE(Intersects(MyPE) % Intersection(Size2,2))
      Intersects(MyPE) % Intersection(:Size,:) = WorkReal(:,1:2)
      DEALLOCATE(WorkReal)

      ALLOCATE(WorkLog(Size,4))
      WorkLog = Intersects(MyPE) % Found
      DEALLOCATE(Intersects(MyPE) % Found)
      ALLOCATE(Intersects(myPE) % Found(Size2,4))
      Intersects(MyPE) % Found(1:Size,:) = WorkLog
      DEALLOCATE(WorkLog)

      ALLOCATE(WorkInt(Size,4))
      WorkInt = Intersects(MyPE) % Elements
      DEALLOCATE(Intersects(MyPE) % Elements)
      ALLOCATE(Intersects(MyPE) % Elements(Size2,4))
      Intersects(MyPE) % Elements(1:Size,:) = WorkInt
      DEALLOCATE(WorkInt)

      Intersects(MyPE) % Size = Size2

    END SUBROUTINE DoubleIntersectsAllocation

    ! double local intersect memmory on the fly
    SUBROUTINE ReduceLocalIntersectMemory()
      !-----------------------------
      REAL(KIND=dp), ALLOCATABLE :: WorkReal(:,:)
      LOGICAL, ALLOCATABLE :: WorkLog(:,:)
      INTEGER, ALLOCATABLE :: WorkInt(:,:)
      INTEGER :: Size,Size2,MyPE

      MyPE = ParEnv % MyPE + 1
      Size = Intersects(MyPE) % Size
      Size2 = Intersects(MyPE) % nIntersects

      ALLOCATE(WorkReal(Size2,4))
      WorkReal = Intersects(MyPE) % LS(1:Size2,:)
      DEALLOCATE(Intersects(MyPE) % LS)
      ALLOCATE(Intersects(myPE) % LS(Size2,4))
      Intersects(MyPE) % LS = WorkReal

      WorkReal(:,1:2) = Intersects(MyPE) % Intersection(1:Size2,:)
      DEALLOCATE(Intersects(MyPE) % Intersection)
      ALLOCATE(Intersects(MyPE) % Intersection(Size2,2))
      Intersects(MyPE) % Intersection = WorkReal(:,1:2)
      DEALLOCATE(WorkReal)

      ALLOCATE(WorkLog(Size2,4))
      WorkLog = Intersects(MyPE) % Found(1:Size2,:)
      DEALLOCATE(Intersects(MyPE) % Found)
      ALLOCATE(Intersects(myPE) % Found(Size2,4))
      Intersects(MyPE) % Found = WorkLog
      DEALLOCATE(WorkLog)

      ALLOCATE(WorkInt(Size,4))
      WorkInt = Intersects(MyPE) % Elements(1:Size2,:)
      DEALLOCATE(Intersects(MyPE) % Elements)
      ALLOCATE(Intersects(MyPE) % Elements(Size2,4))
      Intersects(MyPE) % Elements = WorkInt
      DEALLOCATE(WorkInt)

      Intersects(MyPE) % Size = Size2

    END SUBROUTINE ReduceLocalIntersectMemory

    SUBROUTINE DeallocateMeshLines(Mesh, RemoveLines)
      TYPE(Mesh_t), POINTER :: Mesh
      LOGICAL :: RemoveLines(:,:)
      !------------------------------
      TYPE(Element_t), POINTER :: Element,WorkElements(:)
      INTEGER :: i,MyPE,n

      MyPE = ParEnv % MyPE + 1
      n = Mesh % NumberOfBulkElements

      IF( ALL(.NOT. RemoveLines(MyPe, :n)) ) RETURN ! no change

      ALLOCATE(WorkElements(COUNT(.NOT. RemoveLines(MyPe, :n))))
      WorkElements = PACK(Mesh % Elements, (.NOT. RemoveLines(MyPe,:n)))

      DO i=1, Mesh % NumberOfBulkElements
        IF(RemoveLines(MyPE, i)) THEN
          Element => Mesh % Elements(i)
          IF ( ASSOCIATED( Element % NodeIndexes ) ) &
            DEALLOCATE( Element % NodeIndexes )
          Element % NodeIndexes => NULL()
        END IF
      END DO

      IF(ASSOCIATED(Mesh % Elements)) DEALLOCATE(Mesh % Elements)
      Mesh % Elements => WorkElements
      Mesh % NumberOfBulkElements = SIZE(WorkElements)
      NULLIFY(WorkElements)

    END SUBROUTINE DeallocateMeshLines

    SUBROUTINE DeallocatePolyLines(RemoveLines, nCol)
      LOGICAL :: RemoveLines(:,:)
      INTEGER :: nCol
      !------------------------------
      INTEGER :: i,j,k,p,n,PEs,nMax,counter,nLines,NNeighbours,MyPE
      REAL(KIND=dp), ALLOCATABLE :: Vals(:,:)
      INTEGER, ALLOCATABLE :: WorkInt(:,:), WorkInt2(:,:),PToN(:),Neighbours(:),ns(:)
      LOGICAL, ALLOCATABLE :: WorkLogical(:)

      PEs = ParEnv % PEs
      MyPE = ParEnv % MyPE + 1

      IF(PEs > 1) THEN
        ParEnv % IsNeighbour(MyPE) = .TRUE.
        Neighbours = PACK( (/ (i,i=1,PEs) /), ParEnv % IsNeighbour)
        ParEnv % IsNeighbour(MyPE) = .FALSE.
      ELSE
        ALLOCATE(Neighbours(1))
        Neighbours = MyPE
      END IF

      NNeighbours = SIZE(Neighbours)


      ALLOCATE(PToN(PEs), Ns(PEs))
      PToN = 0

      nMax = 0
      DO i=1, PEs
        IF(nMax < PolylineData(i) % nLines) nMax = PolylineData(i) % nLines
        ns(i) = COUNT(.NOT. RemoveLines(i, :PolylineData(i) % nLines))
        DO j=1,NNeighbours
          IF(i == Neighbours(j)) THEN
            PToN(i) = j
            EXIT
          END IF
        END DO
      END DO

      ALLOCATE(Vals(nMax,nCol), WorkInt(nMax,2), WorkLogical(nMax), &
        WorkInt2(NNeighbours,nMax))

      WorkInt2 = 0
      DO i=1, NNeighbours
        k = Neighbours(i)
        counter = 0
        DO j=1, PolylineData(k) % nLines
          IF(RemoveLines(k,j)) CYCLE
          counter = counter + 1
          WorkInt2(k,j) = counter
        END DO
      END DO

      DO i=1, PEs
        nLines = PolylineData(i) % nLines

        n = Ns(i)

        IF(ALL(.NOT. RemoveLines(i,:nLines))) THEN
          DO j=1, n
            IF (PolylineData(i) % Prev(j,1) > ns(PolylineData(i) % Prev(j,2)) ) &
              PolylineData(i) % Prev(j,1) = 0
            IF (PolylineData(i) % Next(j,1) > ns(PolylineData(i) % Next(j,2)) ) &
              PolylineData(i) % Next(j,1) = 0
          END DO
          CYCLE ! no change
        END IF

        counter = 0
        IF(n > 0) THEN
          IF(PolylineData(i) % nLines > 0) THEN

            DO j=1, nCol
              Vals(1:n,j) = PACK(PolylineData(i) % Vals(:,j), .NOT. RemoveLines(i,:nLines))
            END DO

            DEALLOCATE(PolylineData(i) % Vals)
            ALLOCATE(PolylineData(i) % Vals(n,nCol))
            PolylineData(i) % Vals = Vals(1:n,:)

            WorkInt(1:n,1) = PACK(PolylineData(i) % Prev(:,1), .NOT. RemoveLines(i,:nLines))
            WorkInt(1:n,2) = PACK(PolylineData(i) % Prev(:,2), .NOT. RemoveLines(i,:nLines))
            DEALLOCATE(PolylineData(i) % Prev)
            ALLOCATE(PolylineData(i) % Prev(n,2))
            PolylineData(i) % Prev(:,1) = WorkInt(1:n,1)
            PolylineData(i) % Prev(:,2) = WorkInt(1:n,2)

            WorkInt(1:n,1) = PACK(PolylineData(i) % Next(:,1), .NOT. RemoveLines(i,:nLines))
            WorkInt(1:n,2) = PACK(PolylineData(i) % Next(:,2), .NOT. RemoveLines(i,:nLines))
            DEALLOCATE(PolylineData(i) % Next)
            ALLOCATE(PolylineData(i) % Next(n,2))
            PolylineData(i) % Next(:,1) = WorkInt(1:n,1)
            PolylineData(i) % Next(:,2) = WorkInt(1:n,2)

            ! update next and prev
            DO j=1, n
              IF(PolylineData(i) % Prev(j,1) == -1) CYCLE ! shared node

              IF(PolylineData(i) % Prev(j,1) /= 0) THEN
                k = PToN(PolylineData(i) % Prev(j,2))
                PolylineData(i) % Prev(j,1) = WorkInt2(k, PolylineData(i) % Prev(j,1))
                IF (PolylineData(i) % Prev(j,1) > ns(PolylineData(i) % Prev(j,2)) ) &
                  PolylineData(i) % Prev(j,1) = 0
                !PolylineData(i) % Prev(j,2) = WorkInt2(PolylineData(i) % Prev(j,2))
              END IF
              IF(PolylineData(i) % Next(j,1) /= 0) THEN
                k = PToN(PolylineData(i) % Next(j,2))
                PolylineData(i) % Next(j,1) = WorkInt2(k, PolylineData(i) % Next(j,1))
                IF (PolylineData(i) % Next(j,1) > ns(PolylineData(i) % Next(j,2)) ) &
                  PolylineData(i) % Next(j,1) = 0
                !PolylineData(i) % Next(j,2) = WorkInt2(PolylineData(i) % Next(j,2))
              END IF
            END DO

            PolylineData(i) % nLines = n
          END IF
        END IF
      END DO

    END SUBROUTINE DeallocatePolyLines



    ! use old ls field to determine which direction to remove
    ! positive a1 is +ve compared to a2 on old ls field
    FUNCTION GetPointLevelset(Point, val) RESULT(Found)
      !----------------------------------------
      REAL(KIND=dp) :: Point(3),val
      LOGICAL :: Found
      !----------------------------------------
      REAL(KIND=dp) :: LocalCoords(3),ElementValues(3)
      TYPE(Element_t), POINTER :: HitElement
      LOGICAL :: FirstTime=.TRUE.

      SAVE :: FirstTime

      Found = PointInMesh(Solver, Point, LocalCoords, HitElement, ExtInitialize=FirstTime)!, &
          !CandElement, ExtInitialize )

      IF(.NOT. Found) RETURN

      ElementValues = PhiVar2D % Values(PhiVar2D % Perm(HitElement % NodeIndexes))

      val = InterpolateInElement( HitElement, ElementValues, &
                          LocalCoords(1), LocalCoords(2), LocalCoords(3) )

      FirstTime = .FALSE.

    END FUNCTION GetPointLevelset

    !------------------------------------------------------------------------------
    !> Computes the signed distance to zero levelset. 
    !------------------------------------------------------------------------------
    FUNCTION SignedDistance(node, trusted) RESULT(phip)
      !------------------------------------------------------------------------------
      INTEGER :: node
      REAL(KIND=dp) :: phip
      LOGICAL, INTENT(out) :: trusted
      !-----------------------------------------------------------------------------
      REAL(KIND=dp) :: xp,yp,x0,y0,x1,y1,xm,ym,a,b,c,d,s,dir1,&
          dist2,mindist2,dist,mindist,smin,ss,phim,cosphi,&
          x2,x3,y2,y3,denom,ua,ub,xi,yi,dir_int,sgn_sum,seg_vx,seg_vy,&
          dx,dy,dir_final,dir0
      INTEGER :: i,i0,i1,j,k,n,sgn,m,imin,kmin,dofs,k2,count,next,nextp
      TYPE(Variable_t), POINTER :: Var1D, Var2D
      INTEGER :: nCol, nLines
      REAL(KIND=dp), POINTER :: pValues(:)
      !------------------------------------------------------------------------------
      mindist2 = HUGE(mindist2)
      mindist = HUGE(mindist)
      sgn = 1
      
      xp = Mesh % Nodes % x(node)
      yp = Mesh % Nodes % y(node)
      
      m = 0
      nCol = 7


      DO k = 1, ParEnv % PEs
        nLines = PolylineData(k) % nLines
        IF(nLines == 0) CYCLE

        DO i=1,nLines
          x0 = PolylineData(k) % Vals(i,1)
          x1 = PolylineData(k) % Vals(i,2)
          y0 = PolylineData(k) % Vals(i,3)
          y1 = PolylineData(k) % Vals(i,4)

          a = xp - x0
          b = x0 - x1
          d = y0 - y1
          c = yp - y0

          ! Find the closest distance with the line segment.
          s = -(a*b + c*d) / (b**2 + d**2)
          ! Intersection can not be beyond the element segment.
          s = MIN( MAX( s, 0.0d0), 1.0d0 )
          xm = (1-s) * x0 + s * x1
          ym = (1-s) * y0 + s * y1
          dist2 = (xp - xm)**2 + (yp - ym)**2

          IF(nonzero) THEN
            a = xp - xm
            c = yp - ym            
            ! The signed distance should in a nonzero levelset be computed rather
            ! perpendicular from the line segment.
            cosphi = ABS(a*b + c*d)/SQRT((a**2+c**2)*(b**2+d**2))
            IF(cosphi > cosphi0 ) CYCLE

            ! We need true distances since the offset cannot be added otherwise.
            dist2 = SQRT(dist2)

            ! The line segment including the zero levelset might not be exactly zero...
            ! By definition we don't have permutation here!
            phim = (1-s) * PolylineData(k) % Vals(i,5) + s * PolylineData(k) % Vals(i,6)

            ! In order to test when need to be close enough.
            IF(dist2 > ABS(mindist) + ABS(phim) ) CYCLE

            ! Dir is an indicator one which side of the line segment the point lies. 
            ! We have ordered the edges so that "dir1" should be consistent.
            dir1 = (x1 - x0) * (yp - y0) - (y1 - y0) * (xp - x0)

            ! If the control point and found point lie on the same side they are inside. 
            IF(dir1 < 0.0_dp ) THEN
              sgn = -1
            ELSE
              sgn = 1
            END IF

            dist = sgn * dist2 + phim
            ! Ok, close but no honey. 
            IF( ABS(dist) > ABS(mindist) ) CYCLE

            mindist = dist
          ELSE
            ! Here we can compare the squares saving one expensive operation.
            IF(dist2 > mindist2 ) CYCLE

            ! Dir is an indicator one which side of the line segment the point lies. 
            ! We have ordered the edges soe that "dir1" should be consistent.
            dir1 = (x1 - x0) * (yp - y0) - (y1 - y0) * (xp - x0)

            ! If the control point and found point lie on the same side they are inside. 
            IF(dir1 < 0.0_dp ) THEN
              sgn = -1
            ELSE
              sgn = 1
            END IF
          END IF
          
          ! Save these values for interpolation.
          m = m+1
          mindist2 = dist2          
          smin = s
          imin = i
          kmin = k
          dir_final = dir1
        END DO
      END DO
        
      IF(nonzero) THEN
        phip = mindist        
      ELSE
        phip = sgn * SQRT(mindist2)
      END IF

      trusted = .TRUE.
#if 1
      ! need to have think about how this works for parallel
      IF(PolylineData(kmin) % Prev(imin,1) > 0) THEN
        next = PolylineData(kmin) % Prev(imin,1)
        nextp = PolylineData(kmin) % Prev(imin,2)

        x0 = PolylineData(nextp) % Vals(next,1)
        x1 = PolylineData(nextp) % Vals(next,2)
        y0 = PolylineData(nextp) % Vals(next,3)
        y1 = PolylineData(nextp) % Vals(next,4)
        dir0 = (x1 - x0) * (yp - y0) - (y1 - y0) * (xp - x0)

        ! are the signs different?
        IF (SIGN(1.0_dp, dir0) /= SIGN(1.0_dp, dir_final)) trusted = .FALSE.
      ELSE
        trusted = .FALSE.
      END IF

      IF(trusted) THEN
        IF(PolylineData(kmin) % Next(imin,1) > 0) THEN
          next = PolylineData(kmin) % Next(imin,1)
          nextp = PolylineData(kmin) % Next(imin,2)

          x0 = PolylineData(nextp) % Vals(next,1)
          x1 = PolylineData(nextp) % Vals(next,2)
          y0 = PolylineData(nextp) % Vals(next,3)
          y1 = PolylineData(nextp) % Vals(next,4)
          dir0 = (x1 - x0) * (yp - y0) - (y1 - y0) * (xp - x0)

          ! are the signs different?
          IF (SIGN(1.0_dp, dir0) /= SIGN(1.0_dp, dir_final)) trusted = .FALSE.
        ELSE
          trusted = .FALSE.
        END IF
      END IF
#endif

      ! We can carry the fields with the zero levelset. This is like pure advection.
      ! We should make this less laborious my fetching the pointers first...
      ! Also, do not reinterpolate the nodes that are already ok!
      IF( nVar > 0 .AND. CutPerm(node) == 0 .AND. kmin == Parenv % MyPE+1) THEN !.AND. PhiVar2D % Perm(node) == 0 ) THEN

        i0 = IsoMesh % Elements(imin) % NodeIndexes(1)
        i1 = IsoMesh % Elements(imin) % NodeIndexes(2)

        DO i = 1,nVar
          IF(i==iAvoid) CYCLE
          str = ListGetString( Solver % Values,'isoline variable '//I2S(i), Found )
          
          IF(i==iSolver) THEN
            j = OrigMeshPerm(node)
            dofs = Solver % Variable % Dofs
            pValues => OrigMeshValues
          ELSE
            Var2D => VariableGet( Mesh % Variables, str, ThisOnly = .TRUE. )
            j = Var2D % Perm(node)
            dofs = Var2D % dofs
            pValues => Var2D % Values
          END IF
          
          IF(j==0) THEN
            nCol = nCol+2*dofs
            PRINT *,'We should maybe not be here?:',TRIM(str)
            STOP
          END IF
            
          ! Interpolate from the closest distance.
          ! This is done similarly as the interpolation of coordinates. 
          DO k2 = 1,dofs          
            pValues(dofs*(j-1)+k2) = &
                (1-smin) * PolylineData(kmin) % Vals(imin,nCol) + &
                smin * PolylineData(kmin) % Vals(imin,nCol+1)             
            nCol = nCol+2
          END DO
        END DO
      END IF

      !PRINT *,'phip:',phip, m      

    END FUNCTION SignedDistance
    !------------------------------------------------------------------------------

    SUBROUTINE LineSegmentsIntersect ( a1, a2, b1, b2, intersect_point, does_intersect, buffer)
      ! Find if two 2D line segments intersect
      ! Line segment 'a' runs from point a1 => a2, same for b

      IMPLICIT NONE

      REAL(KIND=dp) :: a1(2), a2(2), b1(2), b2(2), intersect_point(2)
      LOGICAL :: does_intersect
      REAL(KIND=dp), OPTIONAL :: buffer
      !-----------------------
      REAL(KIND=dp) :: r(2), s(2), rxs, bma(2), t, u, err_buffer

      does_intersect = .FALSE.
      intersect_point = 0.0_dp

      r = a2 - a1
      s = b2 - b1

      rxs = VecCross2D(r,s)

      IF(rxs == 0.0_dp) RETURN

      bma = b1 - a1

      t = VecCross2D(bma,s) / rxs
      u = VecCross2D(bma,r) / rxs

      IF(PRESENT(Buffer)) THEN
        err_buffer = Buffer/rxs
      ELSE
        err_buffer = AEPS
      END IF
      IF(t < 0.0_dp-err_buffer .OR. t > 1.0_dp+err_buffer .OR. &
          u < 0.0_dp-err_buffer .OR. u > 1.0_dp+err_buffer) RETURN

      intersect_point = a1 + (t * r)
      does_intersect = .TRUE.

    END SUBROUTINE LineSegmentsIntersect

    FUNCTION VecCross2D(a, b) RESULT (c)
      REAL(KIND=dp) :: a(2), b(2), c

      c = a(1)*b(2) - a(2)*b(1)

    END FUNCTION VecCross2D

  END SUBROUTINE LevelSetUpdate
  
END MODULE CutFemUtils

!> \} ElmerLib

