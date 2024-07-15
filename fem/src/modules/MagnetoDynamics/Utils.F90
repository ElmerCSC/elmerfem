!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
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
! *  Original Date: 08 Jun 1997
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!>  Utilities for the A-V solvers of electromagnetism
!> \ingroup Solvers
!-------------------------------------------------------------------------------
MODULE MagnetoDynamicsUtils

   USE MGDynMaterialUtils
   IMPLICIT NONE

   INTEGER :: JfixPhase
   REAL(KIND=dp), POINTER :: Jfixrhs(:)
   COMPLEX(KIND=dp), POINTER :: JfixRhsC(:)
   INTEGER, POINTER :: JfixSurfacePerm(:) 
   REAL(KIND=dp), ALLOCATABLE, TARGET :: JfixSurfaceVec(:)
   COMPLEX(KIND=dp), ALLOCATABLE :: JfixSurfaceVecC(:)
   
   COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)

   INTERFACE SetDOFtoValue
     MODULE PROCEDURE SetDOFtoValueR, SetDOFtoValueC
   END INTERFACE

   INTERFACE GetReluctivity
     MODULE PROCEDURE GetReluctivityR, GetReluctivityC, &
         GetReluctivityTensorR, GetReluctivityTensorC
   END INTERFACE

CONTAINS

  RECURSIVE FUNCTION AddConstraintFromBulk(A, M0) RESULT (M)
    USE DefUtils
    TYPE(Matrix_t), INTENT(IN) :: A
    TYPE(Matrix_t), POINTER, INTENT(INOUT) :: M0
    TYPE(Matrix_t), POINTER :: M

    integer :: n, n1

    IF(.NOT. A % FORMAT == MATRIX_CRS) CALL Fatal("AddConstraintFromBulk", "Matrix A is not a CRS matrix")

    IF (.NOT. ASSOCIATED(M0)) THEN
      CALL info("AddConstraintFromBulk", "M0 not associated, creating now.", level=10)
      M => AllocateMatrix()
      ! Rowdata copy
      n = A % NumberOfRows
      ALLOCATE(M % RHS(size(A % RHS)))
      ALLOCATE(M % Rows(size(A%rows)))
      ALLOCATE(M % ConstrainedDof(size(M % Rows))) 
      ALLOCATE(M % Diag(size(A % Diag)))

      M % RHS(:) = 0.0_dp
      M % ConstrainedDOF(:) = .FALSE.

      M % ListMatrix => NULL()
      M % NumberOfRows = A % NumberOfRows

      DO n1 = 1, SIZE(M % Rows)
        M % Rows(n1) = A % Rows(n1)
      END DO

      DO n1 = 1,SIZE(M % diag)
        M % Diag(n1) = A % Diag(n1) 
        IF (M % diag(n1) == 0) THEN
          WRITE (*,*) 'diag', n1, 'is zero'
        end if
      END DO

      ! Columndata copy
      n = size(A % Values)
      ALLOCATE(M % Values(n), M % DValues(n), M % Cols(n))
      M % Cols = A % Cols
      M % DValues = 0.0_dp
      M % Values = 0.0_dp

      !! Permutation copy
      IF(ASSOCIATED(A % Perm)) THEN
        CALL info("AddConstraintFromBulk", "Copying perm.", level=3)
        ALLOCATE(M % perm(size(A%perm)))
        M % perm = A % perm
     END IF

      !! Inverse Permutation copy
      IF(ASSOCIATED(A % InvPerm)) THEN
        CALL info("AddConstraintFromBulk", "Copying inverse perm.", level=3)
        ALLOCATE(M % InvPerm(size(A%Invperm)))
        M % InvPerm = A % InvPerm 
      END IF

      M % FORMAT = A % FORMAT
      M0 => M
    ELSE
      CALL info("AddConstraintFromBulk", "M0 is associated, recursing...", level=10)
      M => AddConstraintFromBulk(A, M0 % ConstraintMatrix)
    END IF
  END FUNCTION AddConstraintFromBulk


!------------------------------------------------------------------------------
  SUBROUTINE SetDOFToValueR(Solver,k,VALUE)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp) :: VALUE,v
    TYPE(Solver_t) :: Solver
    INTEGER :: n,k
    TYPE(Matrix_t), POINTER :: A
    TYPE(Mesh_t), POINTER :: Mesh

    Mesh => GetMesh(Solver)
    n = Solver % Variable % Perm(k+Mesh % NumberOfNodes)
    A => GetMatrix()
    CALL CRS_SetSymmDirichlet(A,A % RHS,n,VALUE)
!------------------------------------------------------------------------------
  END SUBROUTINE SetDOFToValueR
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE SetDOFToValueC(Solver,k,VALUE)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    COMPLEX(KIND=dp) :: VALUE
    TYPE(Solver_t) :: Solver
    INTEGER :: n,k
    TYPE(Matrix_t), POINTER :: A
    TYPE(Mesh_t), POINTER :: Mesh

    Mesh => GetMesh(Solver)
    n = Solver % Variable % Perm(k+Mesh % NumberOfNodes)
    A => GetMatrix()

    CALL CRS_SetSymmDirichlet(A,A % RHS,2*(n-1)+1,REAL(VALUE))
    CALL CRS_SetSymmDirichlet(A,A % RHS,2*(n-1)+2,AIMAG(VALUE))
!------------------------------------------------------------------------------
  END SUBROUTINE SetDOFToValueC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE GetReluctivityR(Material,Acoef,n)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: Material
    REAL(KIND=dp) :: Acoef(:)
    INTEGER :: n
!------------------------------------------------------------------------------
    LOGICAL :: Found, FirstTime = .TRUE.
    REAL(KIND=dp) :: Avacuum

    SAVE FirstTime, Avacuum 
!------------------------------------------------------------------------------

    IF ( FirstTime ) THEN
      Avacuum = GetConstReal( CurrentModel % Constants, &
              'Permeability of Vacuum', Found )
      IF (.NOT. Found ) Avacuum = PI * 4.0d-7
      FirstTime = .FALSE.
    END IF
  
    Acoef(1:n) = GetReal( Material, 'Relative Permeability', Found )
    IF ( Found ) THEN
      Acoef(1:n) = Avacuum * Acoef(1:n)
    ELSE
      Acoef(1:n) = GetReal( Material, 'Permeability', Found )
    END IF
    IF ( Found ) THEN
      Acoef(1:n) = 1._dp / Acoef(1:n)
    ELSE
      Acoef(1:n) = GetReal( Material, 'Reluctivity', Found )
    END IF
    IF( .NOT. Found .AND. &
        .NOT. ListCheckPresent(Material, 'H-B Curve') ) THEN
      CALL Fatal('GetReluctivityR','Give > Relative Permeability < or > Reluctivity <  for material!')
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE GetReluctivityR
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
 SUBROUTINE GetReluctivityC(Material,Acoef,n)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: Material
    COMPLEX(KIND=dp) :: Acoef(:)
    INTEGER :: n
!------------------------------------------------------------------------------
    LOGICAL :: L, Found, FirstTime = .TRUE.
    REAL(KIND=dp) :: Avacuum

    SAVE Avacuum, FirstTime
!------------------------------------------------------------------------------

    IF ( FirstTime ) THEN
      Avacuum = GetConstReal( CurrentModel % Constants, &
              'Permeability of Vacuum', Found )
      IF(.NOT. Found ) Avacuum = PI * 4.0d-7
      FirstTime = .FALSE.
    END IF

    Acoef(1:n) = GetReal( Material, 'Relative Permeability', Found )
    IF ( Found ) THEN
      Acoef(1:n) = Avacuum * Acoef(1:n)
    ELSE
      Acoef(1:n) = GetReal( Material, 'Permeability', Found )
    END IF
    IF ( Found ) THEN
      Acoef(1:n) = 1._dp / Acoef(1:n)
    ELSE
      Acoef(1:n) = GetReal( Material, 'Reluctivity', Found )
      Acoef(1:n) = CMPLX( REAL(Acoef(1:n)), &
         GetReal( Material, 'Reluctivity im', L ), KIND=dp )
      Found = Found .OR. L
    END IF
    IF( .NOT. Found .AND. &
        .NOT. ListCheckPresent(Material, 'H-B Curve') ) THEN
      CALL Fatal('GetReluctivityC','Give > Relative Permeability < or > Reluctivity <  for material!')
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE GetReluctivityC
!------------------------------------------------------------------------------

!> Get a real-valued reluctivity tensor. This subroutine seeks values which
!> are strictly given as reluctivity (giving the permeability is not an option
!> here).
!------------------------------------------------------------------------------
  SUBROUTINE GetReluctivityTensorR(Material, Acoef, n, Found)
!-------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER, INTENT(IN) :: Material
    REAL(KIND=dp), POINTER, INTENT(OUT) :: Acoef(:,:,:)
    INTEGER, INTENT(IN) :: n
    LOGICAL, INTENT(OUT) :: Found
!------------------------------------------------------------------------------
    REAL(KIND=dp), SAVE :: nu_vacuum
    LOGICAL, SAVE :: FirstTime = .TRUE.
!------------------------------------------------------------------------------

    IF ( FirstTime ) THEN
      nu_vacuum = GetConstReal( CurrentModel % Constants, &
              'Permeability of Vacuum', Found )
      IF (.NOT. Found ) THEN
        nu_vacuum = 1.0d0/(PI * 4.0d-7)
      ELSE
        nu_vacuum =  1.0d0/nu_vacuum
      END IF
      FirstTime = .FALSE.
    END IF

    CALL GetRealArray( Material, Acoef, 'Reluctivity', Found )

    IF (.NOT. Found) THEN
      CALL GetRealArray( Material, Acoef, 'Relative Reluctivity', Found )
      IF (Found) Acoef = nu_vacuum * Acoef
    END IF
!-------------------------------------------------------------------------------
  END SUBROUTINE GetReluctivityTensorR
!-------------------------------------------------------------------------------

!> Get a complex-valued reluctivity tensor. This subroutine seeks values which
!> are strictly given as reluctivity (giving the permeability is not an option
!> here).
!------------------------------------------------------------------------------
  SUBROUTINE GetReluctivityTensorC(Material, Acoef, n, Found)
!-------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER, INTENT(IN) :: Material
    COMPLEX(KIND=dp), POINTER, INTENT(OUT) :: Acoef(:,:,:)
    INTEGER, INTENT(IN) :: n                                      ! An inactive variable
    LOGICAL, INTENT(OUT) :: Found                                 
!-------------------------------------------------------------------------------
    LOGICAL :: Found_im
    REAL(KIND=dp), POINTER :: work(:,:,:) => NULL()
    INTEGER :: n1, n2, n3
    REAL(KIND=dp), SAVE :: nu_vacuum
    LOGICAL, SAVE :: FirstTime = .TRUE.
!------------------------------------------------------------------------------

    IF ( FirstTime ) THEN
      nu_vacuum = GetConstReal( CurrentModel % Constants, &
              'Permeability of Vacuum', Found )
      IF (.NOT. Found ) THEN
        nu_vacuum = 1.0d0/(PI * 4.0d-7)
      ELSE
        nu_vacuum =  1.0d0/nu_vacuum
      END IF
      FirstTime = .FALSE.
    END IF
    
    IF (ASSOCIATED(Acoef)) DEALLOCATE(Acoef)

    CALL GetRealArray( Material, work, 'Reluctivity', Found )
    IF (.NOT. Found) THEN
      CALL GetRealArray( Material, work, 'Relative Reluctivity', Found )
      IF (Found) work = nu_vacuum * work
    END IF

    IF (Found) THEN
      n1 = SIZE(work,1)
      n2 = SIZE(work,2)
      n3 = SIZE(work,3)
      ALLOCATE(Acoef(n1, n2, n3))
      Acoef(:,:,:) = CMPLX(work(:,:,:), 0.0d0, kind=dp)
    END IF

    CALL GetRealArray( Material, work, 'Reluctivity im', Found_im )
    IF (.NOT. Found_im) THEN
      CALL GetRealArray( Material, work, 'Relative Reluctivity im', Found_im )
      IF (Found_im) work = nu_vacuum * work
    END IF
    
    IF (Found_im) THEN
      n1 = SIZE(work,1)
      n2 = SIZE(work,2)
      n3 = SIZE(work,3)
      IF (.NOT. ASSOCIATED(Acoef)) THEN
        ALLOCATE(Acoef(n1, n2, n3))
        Acoef(:,:,:) = CMPLX(0.0d0, work(:,:,:), kind=dp)
      ELSE
        IF (SIZE(Acoef,1) /= n1 .OR. SIZE(Acoef,2) /= n2 .OR.  SIZE(Acoef,3) /= n3) &
            CALL Fatal('GetReluctivityTensorC', 'Reluctivity and Reluctivity im of different size')
        Acoef(1:n1,1:n2,1:n3) = CMPLX(REAL(Acoef(1:n1,1:n2,1:n3)), work(1:n1,1:n2,1:n3), kind=dp)
      END IF
    END IF
    Found = Found .OR. Found_im

    IF (ASSOCIATED(work)) DEALLOCATE(work)
!-------------------------------------------------------------------------------
  END SUBROUTINE GetReluctivityTensorC
!-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  !> Packs rows associated with edge dofs from constraint matrix and
  !> adds ColOffset to column indices of (1,1) block.
  !> If CM % ConstrainedDOF(i) is true then the ith row is treted empty.
  !! @param ColOffset add to column index of that correspond to nodal dofs
  !!                   by this amount 
  !-------------------------------------------------------------------------------
  SUBROUTINE PackEdgeRows(CM, Model, ColOffset)
    !-------------------------------------------------------------------------------
    TYPE(Matrix_t), INTENT(INOUT), POINTER :: CM 
    TYPE(Model_t), INTENT(IN) :: Model
    INTEGER, OPTIONAL :: ColOffset
    !-------------------------------------------------------------------------------
    LOGICAL, ALLOCATABLE :: Emptyrow(:)
    INTEGER :: ColOffsetF
    REAL(KIND=dp), POINTER CONTIG :: V(:), RHS(:)
    INTEGER, POINTER CONTIG :: R(:), C(:), InvPerm(:), Perm(:)
    INTEGER :: i, j, numempty, k, fileind

    ALLOCATE(EmptyRow(CM % NumberOfRows))
    EmptyRow = .TRUE.

    IF (.NOT. PRESENT(ColOffset)) THEN
      ColOffsetF = CM % NumberOfRows
    ELSE
      ColOffsetF = ColOffset
    END IF

    DO j = 1, CM % NumberOfRows
      DO i = CM % Rows(j), CM % Rows(j+1)-1
        IF (CM % Cols(i) <= Model % NumberOfNodes) &
             CM % Cols(i) = CM % Cols(i) + ColOffsetF
      END DO
    END DO

    numempty = 0
    ROW_LOOP: DO i = 1, CM % NumberOfRows
      ! If CM % ConstrainedDOF(i) is true, then this must be a row corresponding to
      ! edge dof so it must be zero. 
      IF ( CM % ConstrainedDOF(i) ) THEN
        Emptyrow(i) = .TRUE.
      ELSE ! Otherwise the row might correspond with dirichlet scalar dof
        ! Such row must not be packed but diagonal entry of the final matrix
        ! must be set to 1. (Done later?)
        COL_LOOP : DO j = CM % Rows(i), CM % Rows(i+1)-1
          IF (abs(CM % values(j)) >= 1e-12_dp) THEN
            EmptyRow(i) = .FALSE.
            EXIT COL_LOOP
          END IF
        END DO COL_LOOP
      END IF

      IF (EmptyRow(i)) THEN
        CM % RHS(i) = 0.0_dp
        IF (CM % Constraineddof(i)) numempty = numempty + 1
      END IF
    END DO ROW_LOOP


    ! Pack empty rows
    ALLOCATE(R(CM % NumberOfRows-numempty+1))
    ALLOCATE(InvPerm(CM % NumberOfRows-numempty))
    ALLOCATE(Perm(CM % NumberOfRows-numempty))
    ALLOCATE(RHS(CM % NumberOfRows-numempty))
    RHS(:) = 0.0_dp
    k = 2
    R(1) = 1
    DO i = 2, CM % NumberOfRows+1
      IF(.not. emptyrow(i-1)) THEN
        IF(ASSOCIATED(CM%InvPerm)) InvPerm(k-1) = CM % InvPerm(i-1)  !< TODO: This might be nonsense
        IF(ASSOCIATED(CM%Perm)) Perm(k-1) = CM % Perm(i-1)
        R(k) = R(k-1) + CM % Rows(i) - CM % rows(i-1)

        k = k + 1
      END IF
    END DO
    k = size(R)
    ALLOCATE(C(R(k)-1), V(R(k)-1))

    k = 1
    DO i = 2, CM % NumberOfRows+1
      IF (.not. EmptyRow(i-1)) THEN
        DO j = 0, R(k+1)-R(k)-1
          C(R(k)+j) = CM % cols(CM % rows(i-1)+j)
          V(R(k)+j) = CM % Values(CM % rows(i-1)+j)
        END DO
        k = k+1
      END IF
    END DO

    deallocate(CM % Rows)
    deallocate(CM % Cols)
    deallocate(CM % Values)
    deallocate(cm % rhs)
    if(ASSOCIATED(CM % InvPerm)) THEN
      DEALLOCATE(CM % InvPerm)
      CM % InvPerm => InvPerm
    END IF
    IF(ASSOCIATED(CM % Perm)) THEN 
      DEALLOCATE(CM % Perm)
      CM % Perm => Perm
    END IF

    CM % rows => R
    CM % cols => C
    CM % values => V
    CM % RHS => RHS
    CM % NumberOfRows = size(R,1)-1

  !-------------------------------------------------------------------------------
  END SUBROUTINE PackEdgeRows
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  SUBROUTINE CommunicateCondReg(Solver, Mesh, CondReg)
  !-------------------------------------------------------------------------------
     TYPE(Mesh_t) :: Mesh
     LOGICAL :: CondReg(:)
     TYPE(Solver_t) :: Solver

     INTEGER, ALLOCATABLE :: s_e(:,:), r_e(:)
     INTEGER :: i,j,k,l,n,ii(ParEnv % PEs), ierr, status(MPI_STATUS_SIZE)

     IF( ParEnv % PEs<=1 ) RETURN

     ALLOCATE(s_e(Mesh % NumberOfNodes, ParEnv % PEs), &
               r_e(Mesh % NUmberOfNodes) )
     ii = 0
     DO i=1,Mesh % NumberOfNodes
       IF(.NOT.CondReg(i) .AND. Mesh % ParallelInfo % GInterface(i) ) THEN
          DO j=1,SIZE(Mesh % ParallelInfo % Neighbourlist(i) % Neighbours)
            k = Mesh % ParallelInfo % Neighbourlist(i) % Neighbours(j)
            IF ( k== ParEnv % MyPE ) CYCLE
            k = k + 1
            ii(k) = ii(k) + 1
            s_e(ii(k),k) = Mesh % ParallelInfo % GlobalDOFs(i)
          END DO
       END IF
     END DO

     DO i=0, ParEnv % PEs-1
       IF(i==ParEnv % myPE) CYCLE
       k = i+1
       CALL MPI_BSEND( ii(k),1,MPI_INTEGER,i,110,Solver % matrix % comm,ierr )
       IF( ii(k) > 0 ) THEN
         CALL MPI_BSEND( s_e(1:ii(k),k),ii(k),MPI_INTEGER,i,111,solver % matrix % comm,ierr )
       END IF
     END DO

     DO i=0, ParENv % PEs-1
       IF(i==ParEnv % myPE) CYCLE
       CALL MPI_RECV( n,1,MPI_INTEGER,i,110,Solver % Matrix % Comm, status,ierr )
       IF ( n>0 ) THEN
         CALL MPI_RECV( r_e,n,MPI_INTEGER,i,111,Solver % Matrix % Comm, status,ierr )
         DO j=1,n
           k = SearchNode( Mesh % ParallelInfo, r_e(j) )
           IF ( k>0 ) CondReg(k) = .FALSE.
         END DO
       END IF
     END DO
     DEALLOCATE(s_e, r_e)
  !-------------------------------------------------------------------------------
  END SUBROUTINE CommunicateCondReg
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE RecvDoneNodesAndEdges(Solver,Mesh,Done,TreeEdges)
  !-------------------------------------------------------------------------------
     TYPE(Solver_t):: Solver
     TYPE(Mesh_t) :: Mesh
     LOGICAL :: Done(:), TreeEdges(:)
  
     INTEGER, ALLOCATABLE :: s_e(:,:), r_e(:), iperm(:)
     INTEGER :: i,j,k,l,n,ii(ParEnv % PEs), ierr, status(MPI_STATUS_SIZE)

     IF(ParEnv % myPE <= 0) RETURN

     ALLOCATE(r_e(Mesh % NumberOfEdges), iperm(SIZE(Solver % Variable % Perm)))
     iperm = 0
     DO i=1,SIZE(Solver % Variable % Perm)
       IF ( Solver % Variable % Perm(i) > 0 ) &
         iperm(Solver % Variable % Perm(i)) = i
     END DO

     DO i=0, ParEnv % MyPE-1
       CALL MPI_RECV( n,1,MPI_INTEGER,i,112,Solver % Matrix % Comm, status,ierr )
       IF ( n>0 ) THEN
         CALL MPI_RECV( r_e,n,MPI_INTEGER,i,113,Solver % Matrix % Comm, status,ierr )
         DO j=1,n
           k = SearchNode( Solver % Matrix % ParallelInfo, r_e(j), Order=Solver % Variable % Perm )
           k = iperm(k) - Mesh % NumberOfNodes
           IF ( k>0 .AND. k<=SIZE(TreeEdges) ) THEN
             TreeEdges(k) = .TRUE.
           END IF
         END DO
       END IF
     END DO

     DO i=0, ParENv % myPE-1
       CALL MPI_RECV( n,1,MPI_INTEGER,i,114,Solver % Matrix % Comm, status,ierr )
       IF ( n>0 ) THEN
         CALL MPI_RECV( r_e,n,MPI_INTEGER,i,115,Solver % Matrix % Comm, status,ierr )
         DO j=1,n
           k = SearchNode( Mesh % ParallelInfo, r_e(j) )
           IF ( k>0 ) Done(k) = .TRUE.
         END DO
       END IF
     END DO
  !-------------------------------------------------------------------------------
  END SUBROUTINE RecvDoneNodesAndEdges
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  SUBROUTINE SendDoneNodesAndEdges(Solver,Mesh,Done,TreeEdges)
  !-------------------------------------------------------------------------------
     TYPE(Solver_t):: Solver
     TYPE(Mesh_t) :: Mesh
     LOGICAL :: Done(:), TreeEdges(:)
  
     INTEGER, ALLOCATABLE :: s_e(:,:), r_e(:), iperm(:)
     INTEGER :: i,j,k,l,n,ii(ParEnv % PEs), ierr, status(MPI_STATUS_SIZE)

     IF(Parenv % myPE < ParEnv % PEs-1) THEN

      ALLOCATE(s_e(Mesh % NumberOfEdges, ParEnv % PEs) )

      ii = 0
      DO i=1,Mesh % NumberOfedges
        IF ( TreeEdges(i) .AND. Mesh % ParallelInfo % EdgeInterface(i) ) THEN
          DO j=1,SIZE(Mesh % ParallelInfo % EdgeNeighbourlist(i) % Neighbours)
            k = Mesh % ParallelInfo % EdgeNeighbourlist(i) % Neighbours(j)
            IF ( k>ParEnv % myPE ) THEN
               k = k + 1
               ii(k) = ii(k) +1
               s_e(ii(k),k) = Solver % Matrix % ParallelInfo % GlobalDOFs( & 
                    Solver % Variable % Perm(Mesh % NumberOfNodes+i) )
            END IF
          END DO
        END IF
      END DO
 
      DO i=Parenv % mype+1,Parenv % PEs-1
        k = i+1
        CALL MPI_BSEND( ii(k),1,MPI_INTEGER,i,112,Solver % matrix % comm,ierr )
        IF( ii(k) > 0 ) THEN
          CALL MPI_BSEND( s_e(1:ii(k),k),ii(k),MPI_INTEGER,i,113,solver % matrix % comm,ierr )
        END IF
      END DO

      ii = 0
      DO i=1,Mesh % NumberOfNodes
        IF ( Done(i) .AND. Mesh % ParallelInfo % GInterface(i) ) THEN
          DO j=1,SIZE(Mesh % ParallelInfo % Neighbourlist(i) % Neighbours)
            k = Mesh % ParallelInfo % Neighbourlist(i) % Neighbours(j)
            IF ( k>ParEnv % myPE ) THEN
               k = k + 1
               ii(k) = ii(k) +1
               s_e(ii(k),k) = Mesh % ParallelInfo % GlobalDOFs(i)
            END IF
          END DO
        END IF
      END DO

      DO i=Parenv % mype+1,Parenv % PEs-1
        k = i+1
        CALL MPI_BSEND( ii(k),1,MPI_INTEGER,i,114,Solver % matrix % comm,ierr )
        IF( ii(k) > 0 ) THEN
          CALL MPI_BSEND( s_e(1:ii(k),k),ii(k),MPI_INTEGER,i,115,solver % matrix % comm,ierr )
        END IF
      END DO
    END IF

    CALL SparIterBarrier()
  !-------------------------------------------------------------------------------
  END SUBROUTINE SendDoneNodesAndEdges
  !-------------------------------------------------------------------------------

  
  !-------------------------------------------------------------------------------
  ! Mark nodes that are on outer boundary using face elements and node parmutation.
  !-------------------------------------------------------------------------------
  SUBROUTINE MarkOuterNodes(Mesh,Perm,SurfaceNodes,SurfacePerm,EnsureBC) 

    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: Perm(:),SurfacePerm(:)
    INTEGER :: SurfaceNodes
    LOGICAL :: EnsureBC
    
    INTEGER :: snodes0, snodes, i,t,n,ActParents,ParParents
    TYPE(Element_t), POINTER :: Element, P1, P2, P
    LOGICAL, ALLOCATABLE :: BcNode(:)
    
    CALL Info('MarkOuterNodes','Marking outer nodes on outer boundary',Level=8)
    
    SurfaceNodes = 0

    IF( Mesh % NumberOfFaces == 0 ) THEN
      CALL Fatal('MarkOuterNodes','The faces are not created!')
    END IF
    
    n = Mesh % NumberOfNodes
    IF(.NOT. ASSOCIATED( SurfacePerm ) ) THEN
      ALLOCATE( SurfacePerm( n ) )
    END IF
    SurfacePerm = 0
    
       
    DO t=1, Mesh % NumberOfFaces 
      
      Element => Mesh % Faces(t)         
      
      IF( ParEnv % PEs > 1 ) THEN
        ! Don't set BCs on partition interfaces
        IF( Mesh % ParallelInfo % FaceInterface(t) ) CYCLE
      END IF
      
      P1 => Element % BoundaryInfo % Left
      P2 => Element % BoundaryInfo % Right
      
      ActParents = 0
      ParParents = 0

      IF( ASSOCIATED( P1 ) ) THEN
        IF (ALL(Perm(P1 % NodeIndexes)>0)) THEN
          ActParents = ActParents + 1
          IF( P1 % PartIndex == ParEnv % MyPe ) ParParents = ParParents + 1
        ELSE
          NULLIFY( P1 )
        END IF
      END IF
      IF( ASSOCIATED( P2 ) ) THEN
        IF (ALL(Perm(P2 % NodeIndexes)>0)) THEN
          ActParents = ActParents + 1
          IF( P2 % PartIndex == ParEnv % MyPe ) ParParents = ParParents + 1         
        ELSE
          NULLIFY( P2 )
        END IF
      END IF

      ! We have either none or both parents as active.
      ! The BCs will be set only to outer boundaries of the domain. 
      IF( ActParents /= 1 ) CYCLE
      
      ! The one parent is not a true one!
      ! This can happen when we have halo elements. 
      IF( ParEnv % PEs > 0 .AND. ParParents == 0 ) CYCLE
      
      SurfacePerm(Element % NodeIndexes) = 1
    END DO

    IF( EnsureBC ) THEN
      ALLOCATE( BcNode(n) )
      BcNode = .FALSE.

      DO t=Mesh % NumberOfBulkElements +1, &
          Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements      
        Element => Mesh % Elements(t)         
        BcNode(Element % NodeIndexes) = .TRUE.
      END DO

      snodes0 = COUNT( SurfacePerm > 0 )
      snodes0 = ParallelReduction(snodes0) 

      !DO i=1,n
      !  IF( SurfacePerm(i) > 0 .AND. .NOT. BcNode(i) ) THEN
      !    PRINT *,'node:',ParEnv % MyPe, i, Mesh % Nodes % x(i), Mesh % Nodes % y(i), Mesh % Nodes % z(i)
      !  END IF
      !END DO

      WHERE( .NOT. BcNode ) SurfacePerm = 0
      DEALLOCATE( BcNode ) 
    END IF
      
    ! Create numbering for the surface nodes
    snodes = 0
    DO i=1,n
      IF( SurfacePerm(i) > 0 ) THEN
        snodes = snodes + 1
        SurfacePerm(i) = snodes
      END IF
    END DO     
    
    snodes = ParallelReduction(snodes) 
    CALL Info('MarkOuterNodes','Total number of surface nodes: '//I2S(snodes),Level=6)

    IF( EnsureBC ) THEN
      IF( snodes0 > snodes ) THEN
        CALL Info('MarkOuterNodes','Removed number of surface nodes not at BCs: '&
            //I2S(snodes0-snodes),Level=6)
      END IF
    END IF

    SurfaceNodes = snodes
    
  END SUBROUTINE MarkOuterNodes

!------------------------------------------------------------------------------
  SUBROUTINE GaugeTree(Solver,Mesh,TreeEdges,FluxCount,FluxMap,Transient)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Solver_t) :: Solver
    INTEGER :: FluxCount, FluxMap(:)
    LOGICAL :: Transient
    TYPE(Mesh_t) :: Mesh
    LOGICAL, ALLOCATABLE :: TreeEdges(:)

    TYPE(ListMatrixEntry_t), POINTER :: Aentry
    TYPE(ListMatrix_t), POINTER :: Alist(:)
    INTEGER :: i,j,k,l,n,Start
    LOGICAL, ALLOCATABLE :: Done(:), CondReg(:)
    TYPE(ValueList_t), POINTER :: BC
    REAL(KIND=dp) :: Cond1
    TYPE(Element_t), POINTER :: Edge, Boundary, Element

    INTEGER, ALLOCATABLE :: r_e(:), s_e(:,:), iperm(:)
    LOGICAL :: Found
    INTEGER :: ssz, status(MPI_STATUS_SIZE), ierr, ii(ParEnv % PEs)
!------------------------------------------------------------------------------

    IF( ALLOCATED( TreeEdges ) ) THEN
      CALL Info('WhitneyAVSolver','Gauge tree already created',Level=15)
      RETURN
    END IF
      
    ALLOCATE(TreeEdges(Mesh % NumberOfEdges))
    TreeEdges = .FALSE.

    n = Mesh % NumberOfNodes
    ALLOCATE(Done(n)); Done=.FALSE.

    ! Skip Dirichlet BCs in terms of A:
    ! ---------------------------------
    DO i=1,Mesh % NumberOfBoundaryElements
      Boundary => GetBoundaryElement(i)

      SELECT CASE(GetElementFamily())
      CASE(1)
        CYCLE
      CASE(2)
        k = GetBoundaryEdgeIndex(Boundary,1); Element => Mesh % Edges(k)
      CASE(3,4)
        k = GetBoundaryFaceIndex(Boundary)  ; Element => Mesh % Faces(k)
      END SELECT
      IF (.NOT. ActiveBoundaryElement(Element)) CYCLE

      BC => GetBC()
      IF (.NOT.ASSOCIATED(BC)) CYCLE
      IF (.NOT.( ListCheckPresent(BC, 'Mortar BC') .OR. ListCheckPresent( BC, &
                 TRIM(Solver % Variable % Name)//' {e}'))) CYCLE
 
      Done(Element % NodeIndexes) = .TRUE.
    END DO

    IF( Transient ) THEN
      IF ( GetLogical( GetSolverParams(), 'Gauge Tree Skip Conducting Regions', Found) ) THEN
        ! Skip conducting regions:
        ! -------------------------
        ALLOCATE(CondReg(Mesh % NumberOfNodes))
        condReg = .TRUE.
        DO i=1,GetNOFActive()
          Element => GetActiveElement(i)
          Cond1 = GetCReal(GetMaterial(), 'Electric Conductivity',Found)
          IF (cond1==0) condReg(Element % NodeIndexes) = .FALSE.
        END DO

        CALL CommunicateCondReg(Solver,Mesh,CondReg)

        Done = Done.OR.CondReg
        DEALLOCATE(CondReg)
      END IF
    END IF

    ! 
    ! Skip Dirichlet BCs in terms of B:
    ! ---------------------------------
    DO i=1,FluxCount
      j = FluxMap(i)
      IF ( Solver % Variable % Perm(j+n)<=0 ) CYCLE
      Edge => Mesh % Edges(j)
      Done(Edge % NodeIndexes)=.TRUE.
    END DO

    ! 
    ! already set:
    ! ------------

    CALL RecvDoneNodesAndEdges(Solver,Mesh,Done,TreeEdges)

    ! node -> edge list
    ! -----------------
    Alist => NULL()
    n = Mesh % NumberOfNodes
    DO i=1,Mesh % NumberOfEdges
      Edge => Mesh % Edges(i)
      IF ( Solver % Variable % Perm(i+n)<=0 ) CYCLE
      DO j=1,Edge % TYPE % NumberOfNodes
        k=Edge % NodeIndexes(j)
        Aentry=>List_GetMatrixIndex(Alist,k,i)
      END DO
    END DO

    !
    ! generate the tree for all (perhaps disconnected) parts:
    ! -------------------------------------------------------
    DO WHILE(.NOT.ALL(Done))
      DO Start=1,n
        IF (.NOT. Done(Start)) EXIT
      END DO
      CALL DepthFirstSearch(Alist,Done,Start)
    END DO
    CALL List_FreeMatrix(SIZE(Alist),Alist)

    CALL SendDoneNodesAndEdges(Solver,Mesh,Done,TreeEdges)
    DEALLOCATE(Done)

CONTAINS

!------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE DepthFirstSearch(Alist,done,i)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ListMatrix_t) :: Alist(:)
    INTEGER :: i
    LOGICAL :: Done(:)
!------------------------------------------------------------------------------
    TYPE(ListMatrixEntry_t), POINTER :: Aentry
    INTEGER :: j,k,l,n
    TYPE(Element_t), POINTER :: Edge
!------------------------------------------------------------------------------

    ! To give better matrix conditioning some directional heuristics
    ! could be added,e.g. select the order of going through the nodes
    ! edge list here:

    Done(i) = .TRUE.

    Aentry => Alist(i) % Head
    DO WHILE(ASSOCIATED(Aentry))
      k = Aentry % Index
      Aentry => Aentry % Next

      Edge => Mesh % Edges(k)
      IF (ALL(Done(Edge % NodeIndexes))) CYCLE

      IF ( .NOT. TreeEdges(k)) CALL SetDOFToValue(Solver,k,0._dp)
      TreeEdges(k)=.TRUE.
      DO l=1,2
        n = Edge % NodeIndexes(l)
        IF (.NOT. Done(n)) CALL DepthFirstSearch(Alist,done,n)
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE DepthFirstSearch
!------------------------------------------------------------------------------


    
!------------------------------------------------------------------------------
  END SUBROUTINE GaugeTree
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE GaugeTreeFluxBC(Solver,Mesh,TreeEdges,BasicCycles,FluxCount,FluxMap)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ListMatrix_t), POINTER :: BasicCycles(:)
    INTEGER :: FluxCount, FluxMap(:)
    TYPE(Solver_t) :: Solver
    TYPE(Mesh_t) :: Mesh
    LOGICAL, ALLOCATABLE :: TreeEdges(:)
!------------------------------------------------------------------------------
!   TYPE(Mesh_t) :: Mesh
!   LOGICAL, ALLOCATABLE :: TreeEdges(:)
!------------------------------------------------------------------------------
    TYPE(ListMatrixEntry_t), POINTER :: Aentry, Ltmp
    TYPE(ValueList_t), POINTER :: BC
    TYPE(ListMatrix_t), POINTER :: Alist(:)
    INTEGER :: i,j,k,l,n,Start,nCount,fixedge
    LOGICAL, ALLOCATABLE :: Done(:)
    INTEGER, ALLOCATABLE :: NodeList(:)
    TYPE(Element_t), POINTER :: Edge, Boundary, Element
!------------------------------------------------------------------------------

    IF( ALLOCATED( TreeEdges ) ) THEN
      CALL Info('WhitneyAVSolver','Boundary Gauge tree already created',Level=15)
      RETURN
    END IF

    ALLOCATE(TreeEdges(Mesh % NumberOfEdges))
    TreeEdges = .FALSE.

    n = Mesh % NumberOfNodes
    ALLOCATE(Done(n)); Done=.FALSE.

    !
    ! list the candidate nodes:
    ! -------------------------
    DO i=1,FluxCount
      j = FluxMap(i)
      Edge => Mesh % Edges(j)
      Done(Edge % NodeIndexes)=.TRUE.
    END DO

    ALLOCATE(NodeList(COUNT(Done)))
    nCount = 0
    DO i=1,n
      IF ( Done(i) ) THEN
        nCount = nCount+1
        NodeList(nCount)=i
      END IF
    END DO

    Done=.FALSE.
    DO i=1,FluxCount
      IF ( TreeEdges(FluxMap(i)) ) THEN
        Edge => Mesh % Edges(FluxMap(i))
        Done(Edge % NodeIndexes)=.TRUE.
      END IF
    END DO

    ! 
    ! Skip Dirichlet BCs in terms of A:
    ! ---------------------------------
    DO i=1,Mesh % NumberOfBoundaryElements
      Boundary => GetBoundaryElement(i)
      SELECT CASE(GetElementFamily())
      CASE(1)
        CYCLE
      CASE(2)
        k = GetBoundaryEdgeIndex(Boundary,1); Element => Mesh % Edges(k)
      CASE(3,4)
        k = GetBoundaryFaceIndex(Boundary)  ; Element => Mesh % Faces(k)
      END SELECT
      IF (.NOT. ActiveBoundaryElement(Element)) CYCLE
      BC => GetBC()
      IF (.NOT.ASSOCIATED(BC)) CYCLE
      IF (.NOT.ListCheckPresent( BC, &
           TRIM(Solver % Variable % Name)//' {e}')) CYCLE
 
      j=1; k=GetBoundaryEdgeIndex(Boundary,j)
      DO WHILE(k>0)
        Edge => Mesh % Edges(k)
        TreeEdges(k) = .TRUE.
        Done(Edge % NodeIndexes) = .TRUE.
        j=j+1; k=GetBoundaryEdgeIndex(Boundary,j)
      END DO
    END DO

    ! node -> edge list
    ! -----------------
    Alist => NULL()
    DO i=1,FluxCount
      j = FluxMap(i)
      IF ( Solver % Variable % Perm(j+n)<=0 ) CYCLE

      Edge => Mesh % Edges(j)
      DO k=1,Edge % TYPE % NumberOfNodes
        l=Edge % NodeIndexes(k)
        Aentry=>List_GetMatrixIndex(Alist,l,j)
      END DO
    END DO

    ALLOCATE(BasicCycles(FluxCount))
    BasicCycles(:) % Degree = 0
    DO i=1,FluxCount
       BasicCycles(i) % Head => NULL()
    END DO

    ! generate the tree for all (perhaps disconnected) parts:
    ! -------------------------------------------------------
    DO WHILE(.NOT.ALL(Done(NodeList)))
      DO i=1,nCount
        Start = NodeList(i)
        IF ( .NOT. Done(Start) ) EXIT
      END DO
      CALL BreadthFirstSearch(Alist,Done,Start,nCount,NodeList)
    END DO
    DEALLOCATE(Done,NodeList)
    CALL List_FreeMatrix(SIZE(Alist),Alist)

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE BreadthFirstSearch(Alist,done,start,nCount,NodeList)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: start,nCount,NodeList(:)
    LOGICAL :: Done(:)
    TYPE(ListMatrix_t) :: Alist(:)
!------------------------------------------------------------------------------
    TYPE(ListMatrixEntry_t), POINTER :: Aentry, Ltmp, Btmp
    INTEGER :: i,j,k,l,n,m,ll,IF,bcycle
    LOGICAL :: FirstTime= .TRUE.
    TYPE(Element_t), POINTER :: Edge,Edge1,Boundary
    LOGICAL, ALLOCATABLE :: DoneL(:)
    INTEGER, ALLOCATABLE :: Fifo(:), Previous(:), FiFo1(:)
!------------------------------------------------------------------------------

   ALLOCATE(DoneL(Mesh % NumberOfEdges)); DoneL=.FALSE.
   ALLOCATE(Fifo(FluxCount),FiFo1(FluxCount))
   ALLOCATE(Previous(Mesh % NumberOfNodes)); Previous=0;

   bcycle = 0
   DO bcycle=0,FluxCount
     IF(.NOT.ASSOCIATED(BasicCycles(bcycle+1) % Head)) EXIT
   END DO

   IF = 0; m=0

   IF( FirstTime ) THEN
     DO i=1,nCount
       j = NodeList(i)
       IF ( Done(j) ) THEN
         m=m+1; fifo1(m)=j
         IF=IF+1; fifo(IF)=j
       END IF
     END DO
     FirstTime = .FALSE.
   END IF

   IF ( IF==0 ) THEN
     Done(Start)=.TRUE.
     m=m+1; fifo1(m)=start
     IF=1; fifo(IF)=start;
   END IF

   IF ( IF>0 ) THEN
     DO WHILE(m>0)
       j = Fifo1(m); m=m-1

       Aentry => Alist(j) % Head
       DO WHILE(ASSOCIATED(Aentry))
         k = Aentry % Index
         Aentry => Aentry % Next

         Edge => Mesh % Edges(k)
         IF (.NOT. TreeEdges(k) .OR. DoneL(k) ) CYCLE
         DoneL(k)=.TRUE.

         l = Edge % NodeIndexes(1)
         IF (l==j) l=Edge % NodeIndexes(2)

         IF=IF+1; Fifo(IF)=l
         m=m+1; Fifo1(m)=l
         Previous(l)=j
       END DO
     END DO
     Start = l
   END IF


   DO WHILE(IF>0)
     j = Fifo(IF); IF=IF-1

     Aentry => Alist(j) % Head
     DO WHILE(ASSOCIATED(Aentry))
       k = Aentry % Index
       Aentry => Aentry % Next

       Edge => Mesh % Edges(k)
       IF ( DoneL(k) ) CYCLE
       DoneL(k)=.TRUE.

       l = Edge % NodeIndexes(1)
       IF (l==j) l=Edge % NodeIndexes(2)

       IF ( Done(l) ) THEN
         ! Generate fundamental cycle
         bcycle = bcycle+1
         CALL AddToCycle(bcycle,k)

         m = j
         DO WHILE(m/=Previous(l))
           Ltmp => Alist(m) % Head
           DO WHILE(ASSOCIATED(Ltmp))
             Edge1 => Mesh % Edges(Ltmp % Index)
             IF ( ANY(Edge1 % NodeIndexes(1:2)==Previous(m)) ) THEN
               CALL AddToCycle(bcycle,Ltmp % Index); EXIT
             END IF
             Ltmp=>Ltmp % Next
           END DO
           IF (ANY(Edge1 % NodeIndexes(1:2) == l) ) EXIT
           m = Previous(m)
           IF(m==0) EXIT
         END DO

         IF (ALL(Edge1 % NodeIndexes(1:2) /= l) ) THEN
           ltmp => Alist(l) % Head
           DO WHILE(ASSOCIATED(ltmp))
             edge1 => Mesh % Edges(Ltmp % Index)
             IF ( ANY(Edge1 % NodeIndexes(1:2)==Previous(l)) ) THEN
               CALL AddToCycle(bcycle,Ltmp % Index); EXIT
             END IF
             ltmp=>ltmp % Next
           END DO
         END IF
       ELSE
         IF (.NOT.TreeEdges(k)) CALL SetDOFToValue(Solver,k,0._dp)
         IF=IF+1; Fifo(IF)=l
         Previous(l)=j
         Done(l)=.TRUE.
         TreeEdges(k) = .TRUE.
       END IF
     END DO
   END DO
   DEALLOCATE(Fifo, Fifo1, DoneL)
!------------------------------------------------------------------------------
  END SUBROUTINE BreadthFirstSearch
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AddToCycle(bcycle,index)
    IMPLICIT NONE
    INTEGER :: bcycle,index
!------------------------------------------------------------------------------
    TYPE(ListMatrixEntry_t), POINTER :: Btmp

    ALLOCATE(Btmp); Btmp % Next => BasicCycles(bcycle) % Head;
    Btmp % Index = index; BasicCycles(bcycle) % Head => Btmp
    BasicCycles(bcycle) % Degree=BasicCycles(bcycle) % Degree+1
!------------------------------------------------------------------------------
  END SUBROUTINE AddToCycle
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  END SUBROUTINE GaugeTreeFluxBC
!------------------------------------------------------------------------------

END MODULE MagnetoDynamicsUtils
