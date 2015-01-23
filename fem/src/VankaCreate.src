!/*****************************************************************************/
! *
! * Elmer, A Finite Element Software for Multiphysical Problems
! *
! * Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
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

!-------------------------------------------------------------------------------
!> Vanka preconditioning for iterative methods.
!-------------------------------------------------------------------------------
    SUBROUTINE VankaPrec(u,v,ipar)
!-------------------------------------------------------------------------------
      USE DefUtils

      INTEGER :: ipar(*)
      REAL(KIND=dp) u(*), v(*)
!-------------------------------------------------------------------------------
      TYPE(SplittedMatrixT), POINTER :: SP
      TYPE(Matrix_t), POINTER :: A
      INTEGER :: i
      REAL(KIND=dp), POINTER :: SaveValues(:)
      TYPE(BasicMatrix_t), POINTER :: SaveIF(:)
!-------------------------------------------------------------------------------
      A => GlobalMatrix
      SaveValues => A % Values
      A % Values => A % ILUValues

      IF (ParEnv % Pes <= 1 ) THEN
        CALL CRS_MatrixVectorProd(v,u,ipar)
      ELSE
        SP => GlobalData % SplittedMatrix
        ALLOCATE( SaveIF(ParEnv % PEs) )

        DO i=1,ParEnv % PEs
          IF ( SP % IfMatrix(i) % NumberOfRows /= 0 ) THEN
            ALLOCATE(SaveIF(i) % Values(SIZE(SP % IfMatrix(i) % Values)))
            SaveIF(i) % Values = SP % IfMatrix(i) % Values
            SP % IfMatrix(i) % Values = SP % IfMatrix(i) % ILUValues
          END IF
        END DO

        CALL SParMatrixVector(v,u,ipar)

        DO i=1,ParEnv % PEs
          IF ( SP % IfMatrix(i) % NumberOfRows /= 0 ) THEN
            SP % IfMatrix(i) % Values = SaveIF(i) % Values
            DEALLOCATE(SaveIf(i) % Values)
          END IF
        END DO

        DEALLOCATE(SaveIF)
      END IF

      A % Values => SaveValues
!-------------------------------------------------------------------------------
    END SUBROUTINE VankaPrec
!-------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Create the Vanka preconditioning.
!------------------------------------------------------------------------------
  SUBROUTINE VankaCreate(A,Solver)
     USE DefUtils
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     TYPE(Solver_t) :: Solver
!------------------------------------------------------------------------------
     INTEGER, POINTER :: Diag(:), Rows(:), Cols(:), Perm(:), Indexes(:), Ind(:)
     REAL(KIND=dp), POINTER :: ILUValues(:), SValues(:), TotValues(:)
     REAL(KIND=dp), ALLOCATABLE :: al(:,:)
     LOGICAL ::  found
     TYPE(Element_t), POINTER :: Element
     INTEGER :: status(MPI_STATUS_SIZE)
     INTEGER :: i,j,k,l,m,proc,rcnt,nn, dof, dofs, Active
     REAL(KIND=dp), ALLOCATABLE, TARGET :: rval(:)
     INTEGER, ALLOCATABLE :: cnt(:), rrow(:),rcol(:)

     TYPE Buf_t
        REAL(KIND=dp), ALLOCATABLE :: gval(:)
        INTEGER, ALLOCATABLE :: grow(:),gcol(:)
     END TYPE Buf_t
     TYPE(Buf_t), POINTER :: buf(:)

     Diag => A % Diag
     Rows => A % Rows
     Cols => A % Cols

     IF ( .NOT. ASSOCIATED(A % ILUValues) ) THEN
       ALLOCATE( A % ILUvalues(SIZE(A % Values)) )
     END IF
     ILUValues => A % ILUValues

     Dofs  = Solver % Variable % DOFs
     Perm => Solver % Variable % Perm

     Svalues => A % Values

     ALLOCATE(TotValues(SIZE(A % Values)))
     IF(ASSOCIATED(A % PrecValues)) THEN
       TotValues = A % PrecValues
     ELSE
       TotValues = A % Values
     END IF

     IF ( ParEnv  % PEs>1 ) THEN
       ALLOCATE(cnt(0:ParEnv % PEs))
       cnt = 0
       DO i=1,A % NumberOfRows
         DO j=Rows(i),Rows(i+1)-1
           IF ( A % ParallelInfo % Interface(Cols(j)) ) THEN
             DO l=1,SIZE(A % ParallelInfo % NeighbourList(Cols(j)) % Neighbours)
               m = A % ParallelInfo % NeighbourList(Cols(j)) % Neighbours(l)
               IF ( m==ParEnv % myPE ) CYCLE
               cnt(m) = cnt(m)+1
             END DO
           END IF
         END DO
       END DO

       ALLOCATE( buf(0:ParEnv % PEs-1) )
       DO i=0,ParEnv % PEs-1
         IF ( cnt(i) > 0 ) &
           ALLOCATE( Buf(i) % gval(cnt(i)), Buf(i) % grow(cnt(i)), Buf(i) % gcol(cnt(i)) )
       END DO

       cnt = 0
       DO i=1,A % NumberOfRows
         DO j=Rows(i),Rows(i+1)-1
           DO l=1,SIZE(A % ParallelInfo % NeighbourList(Cols(j)) % Neighbours)
             m = A % ParallelInfo % NeighbourList(Cols(j)) % Neighbours(l)
             IF ( m==ParEnv % myPE ) CYCLE
             cnt(m) = cnt(m)+1
             Buf(m) % gcol(cnt(m)) = A % ParallelInfo % GlobalDOFs(Cols(j))
             Buf(m) % gval(cnt(m)) = TotValues(j)
             Buf(m) % grow(cnt(m)) = A % ParallelInfo % GlobalDOFs(i)
           END DO
         END DO
       END DO

       DO i=0,ParEnv % PEs-1
         IF ( ParEnv % IsNeighbour(i+1) ) THEN
           CALL MPI_BSEND( cnt(i), 1, MPI_INTEGER, &
               i, 7001, MPI_COMM_WORLD, status, ierr )
           IF ( cnt(i)>0 ) THEN
             CALL MPI_BSEND( Buf(i) % grow, cnt(i), MPI_INTEGER, &
                 i, 7002, MPI_COMM_WORLD, status, ierr )

             CALL MPI_BSEND( Buf(i) % gcol, cnt(i), MPI_INTEGER, &
                 i, 7003, MPI_COMM_WORLD, status, ierr )

             CALL MPI_BSEND( Buf(i) % gval, cnt(i), MPI_DOUBLE_PRECISION, &
                 i, 7004, MPI_COMM_WORLD, status, ierr )
           END IF
         END IF
       END DO

       DO i=0,ParEnv % PEs-1
         IF ( cnt(i)>0 ) &
           DEALLOCATE( Buf(i) % gval, Buf(i) % grow, Buf(i) % gcol )
       END DO
       DEALLOCATE( cnt,Buf )

       k = SIZE(A % Values)
       ALLOCATE( rrow(k), rcol(k), rval(k) )

       DO i=1,ParEnv % NumOfNeighbours
         CALL MPI_RECV( rcnt, 1, MPI_INTEGER, &
           MPI_ANY_SOURCE, 7001, MPI_COMM_WORLD, status, ierr )

         IF ( rcnt>0 ) THEN
           proc = status(MPI_SOURCE)
           CALL MPI_RECV( rrow, rcnt, MPI_INTEGER, &
              proc, 7002, MPI_COMM_WORLD, status, ierr )

           CALL MPI_RECV( rcol, rcnt, MPI_INTEGER, &
              proc, 7003, MPI_COMM_WORLD, status, ierr )

           CALL MPI_RECV( rval, rcnt, MPI_DOUBLE_PRECISION, &
              proc, 7004, MPI_COMM_WORLD, status, ierr )

           DO j=1,rcnt
             l = SearchNode(A % ParallelInfo,rcol(j),Order=A % Perm)
             IF ( l>0 ) THEN
               k = SearchNode(A % ParallelInfo,rrow(j),Order=A % Perm)
               IF ( k>0 ) THEN
                 IF ( l>=k ) THEN
                   DO m=Diag(k),Rows(k+1)-1
                     IF ( Cols(m)==l ) THEN
                       TotValues(m)=TotValues(m)+rval(j)
                       EXIT
                     END IF
                   END DO
                 ELSE
                   DO m=Rows(k),Diag(k)
                     IF ( Cols(m)==l ) THEN
                       TotValues(m)=TotValues(m)+rval(j)
                       EXIT
                     END IF
                   END DO
                 END IF
               END IF
             END IF
           END DO
         END IF
       END DO
       DEALLOCATE(rrow,rcol,rval)
     END IF

     nn = Solver % Mesh % MaxElementDOFs
     ALLOCATE( Indexes(nn), AL(nn*dofs,nn*dofs), ind(nn*dofs) )

     Active = GetNOFActive(Solver)
     DO i=1,Active
       element => GetActiveElement(i)
       nn = GetElementDOFs(Indexes)

       l = 0
       DO j=1,nn
         k =  Indexes(j)
         DO dof=1,dofs
           l = l + 1
           ind(l) = DOFs*(perm(k)-1)+dof
         END DO
       END DO

       A % Values => TotValues
       DO j=1,nn*dofs
         DO k=1,nn*dofs
           al(j,k) = CRS_GetMatrixElement( A, ind(j), ind(k) )
         END DO
       END DO
       CALL InvertMatrix(al,nn*dofs)
       A % Values => ILUValues
       DO j=1,nn*dofs
         DO k=1,nn*dofs
           CALL CRS_AddToMatrixElement( A,ind(j),ind(k),AL(j,k) )
         END DO
       END DO
     END DO
     A % Values => Svalues
     DEALLOCATE(AL, Indexes, Ind, TotValues)
!------------------------------------------------------------------------------
  END SUBROUTINE VankaCreate
!------------------------------------------------------------------------------

!> \}
