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
!>  Solve Maxwell equations in vector potential formulation (or the A-V
!>  formulation) and (relatively)low frequency approximation using lowest
!>  order Withney 1-forms (edge elements).
!> \ingroup Solvers
!-------------------------------------------------------------------------------
MODULE MagnetoDynamicsUtils

   USE MGDynMaterialUtils
   IMPLICIT NONE

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
          write (*,*), 'diag', n1, 'is zero'
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
  FUNCTION GetBoundaryEdgeIndex(Boundary,nedge) RESULT(n)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: n,nedge
    TYPE(Element_t) :: Boundary
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,jb1,jb2,je1,je2
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: Parent, Edge, Face
!------------------------------------------------------------------------------
    Mesh => GetMesh()
    n = 0
    SELECT CASE(GetElementFamily(Boundary))
    CASE(1)
      RETURN
    CASE(2)
      IF ( nedge==1 ) THEN
        Parent => Boundary % BoundaryInfo % Left
        IF ( .NOT. ASSOCIATED(Parent) ) &
            Parent => Boundary % BoundaryInfo % Right
 
        jb1 = Boundary % NodeIndexes(1)
        jb2 = Boundary % NodeIndexes(2)
        DO i=1,Parent % TYPE % NumberOfEdges
          Edge => Mesh % Edges(Parent % EdgeIndexes(i))
          je1 = Edge % NodeIndexes(1)
          je2 = Edge % NodeIndexes(2)
          IF ( jb1==je1.AND.jb2==je2 .OR. jb1==je2.AND.jb2==je1) EXIT
        END DO
        n = Parent % EdgeIndexes(i)
      END IF
    CASE(3,4)
      j = GetBoundaryFaceIndex(Boundary)
      Face => Mesh % Faces(j)
      IF ( nedge>0.AND.nedge<=Face % TYPE % NumberOfEdges ) &
        n = Face % EdgeIndexes(nedge) 
    END SELECT
!------------------------------------------------------------------------------
  END FUNCTION GetBoundaryEdgeIndex
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  FUNCTION GetBoundaryFaceIndex(Boundary) RESULT(n)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: n
    TYPE(Element_t) :: Boundary
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,m
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: Parent, Face
!------------------------------------------------------------------------------
    Mesh => GetMesh()
    Parent => Boundary % BoundaryInfo % Left
    IF ( .NOT. ASSOCIATED(Parent) ) &
       Parent => Boundary % BoundaryInfo % Right

    DO i=1,Parent % TYPE % NumberOfFaces
      Face => Mesh % Faces(Parent % FaceIndexes(i))
      m = 0
      DO j=1,Face % TYPE % NumberOfNodes
        DO k=1,Boundary % TYPE % NumberOfNodes
          IF ( Face % NodeIndexes(j)==Boundary % NodeIndexes(k)) m=m+1
        END DO
      END DO
      IF ( m==Face % TYPE % NumberOfNodes) EXIT
    END DO
    n = Parent % FaceIndexes(i)
!------------------------------------------------------------------------------
  END FUNCTION GetBoundaryFaceIndex
!------------------------------------------------------------------------------


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
    INTEGER :: n
    REAL(KIND=dp) :: Acoef(:)
!------------------------------------------------------------------------------
    LOGICAL :: Found, FirstTime = .TRUE., Warned = .FALSE.
    REAL(KIND=dp) :: Avacuum

    SAVE Avacuum 

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
    END IF
    IF( .NOT. Found .AND. .NOT. Warned .AND. &
        .NOT. ListCheckPresent(Material, 'H-B Curve') ) THEN
      CALL Warn('GetReluctivityR','Give > Relative Permeability < or > Reluctivity <  for material!')
      Warned = .TRUE.
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE GetReluctivityR
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
 SUBROUTINE GetReluctivityC(Material,Acoef,n)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: Material
    INTEGER :: n
    COMPLEX(KIND=dp) :: Acoef(:)
!------------------------------------------------------------------------------
    LOGICAL :: L, Found, FirstTime = .TRUE., Warned = .FALSE.
    REAL(KIND=dp) :: Avacuum

    SAVE Avacuum 

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
    IF( .NOT. Found .AND. .NOT. Warned .AND. &
        .NOT. ListCheckPresent(Material, 'H-B Curve') ) THEN
      CALL Warn('GetReluctivityC','Give > Relative Permeability < or > Reluctivity <  for material!')
      Warned = .TRUE.
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE GetReluctivityC
!------------------------------------------------------------------------------

!> Get real tensorial reluctivity
!------------------------------------------------------------------------------
  SUBROUTINE GetReluctivityTensorR(Material, Acoef, n, Found)
!-------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER, INTENT(IN) :: Material
    REAL(KIND=dp), POINTER :: Acoef(:,:,:)
    INTEGER, INTENT(IN) :: n
    LOGICAL , INTENT(OUT) :: Found
!-------------------------------------------------------------------------------
    LOGICAL :: FirstTime = .FALSE.
    INTEGER :: k
    REAL(KIND=dp) :: Avacuum

    SAVE Avacuum

    CALL GetRealArray( Material, Acoef, 'Relative Reluctivity', Found )
!-------------------------------------------------------------------------------
  END SUBROUTINE GetReluctivityTensorR
!-------------------------------------------------------------------------------

!> Get complex tensorial reluctivity
!> Untested
!------------------------------------------------------------------------------
  SUBROUTINE GetReluctivityTensorC(Material, Acoef, n, Found, Cwrk)
!-------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER, INTENT(IN) :: Material
    COMPLEX(KIND=dp), POINTER :: Acoef(:,:,:)
    REAL(KIND=dp), POINTER, OPTIONAL :: Cwrk(:,:,:)
    INTEGER, INTENT(IN) :: n
    LOGICAL , INTENT(OUT) :: Found
!-------------------------------------------------------------------------------
    LOGICAL :: FirstTime = .FALSE.
    LOGICAL :: Found_im
    INTEGER :: k1,k2,k3
    REAL(KIND=dp) :: Avacuum
    REAL(KIND=dp), POINTER :: work(:,:,:)

    SAVE Avacuum

    IF(.NOT. PRESENT(Cwrk)) THEN
      ALLOCATE(work(size(Acoef,1), size(Acoef,2), size(Acoef,3)))
    ELSE
      work => Cwrk
    END IF


    CALL GetRealArray( Material, work, 'Relative Reluctivity', Found )
    Acoef(:,:,:) = work(:,:,:)

    CALL GetRealArray( Material, work, 'Relative Reluctivity im', Found_im )

    Acoef = CMPLX(REAL(Acoef), work)

    Found = Found .OR. Found_im

    IF(.NOT. PRESENT(Cwrk)) THEN
      DEALLOCATE(work)
    END IF
!-------------------------------------------------------------------------------
  END SUBROUTINE GetReluctivityTensorC
!-------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE GetPermittivity(Material,Acoef,n)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: Material
    INTEGER :: n
    REAL(KIND=dp) :: Acoef(:)
!------------------------------------------------------------------------------
    LOGICAL :: Found, FirstTime = .TRUE., Warned = .FALSE.
    REAL(KIND=dp) :: Pvacuum = 0._dp

    IF ( FirstTime ) THEN
      Pvacuum = GetConstReal( CurrentModel % Constants, &
              'Permittivity of Vacuum', Found )
      FirstTime = .FALSE.
    END IF
    

    Acoef(1:n) = GetReal( Material, 'Relative Permittivity', Found )
    IF ( Found ) THEN
      Acoef(1:n) = Pvacuum * Acoef(1:n)
    ELSE
      Acoef(1:n) = GetReal( Material, 'Permittivity', Found )
    END IF

    IF( .NOT. Found ) THEN
      IF(.NOT. Warned ) THEN
        CALL Warn('GetPermittivity','Permittivity not defined in material, defaulting to that of vacuum')
        Warned = .TRUE.
      END IF
      Acoef(1:n) = Pvacuum
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE GetPermittivity
!------------------------------------------------------------------------------

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
      ! edge dof so it must tbe zero. 
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
       IF(.NOT.CondReg(i) .AND. Mesh % ParallelInfo % Interface(i) ) THEN
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
        IF ( Done(i) .AND. Mesh % ParallelInfo % Interface(i) ) THEN
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

END MODULE MagnetoDynamicsUtils
