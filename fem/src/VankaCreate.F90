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
      REAL(KIND=dp), POINTER CONTIG :: SaveValues(:)
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
     REAL(KIND=dp), POINTER CONTIG :: ILUValues(:), SValues(:), TotValues(:)
     REAL(KIND=dp), ALLOCATABLE :: al(:,:)
     LOGICAL ::  found
     TYPE(Element_t), POINTER :: Element
     INTEGER :: status(MPI_STATUS_SIZE)
     INTEGER :: i,j,k,l,m,proc,rcnt,nn, dof, dofs, Active, Totcnt
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

       totcnt = SUM(cnt)
       CALL CheckBuffer( ParEnv % PEs*(1+MPI_BSEND_OVERHEAD) + 4*totcnt + &
                  3*COUNT(cnt/=0)*MPI_BSEND_OVERHEAD)

       DO i=0,ParEnv % PEs-1
         IF ( ParEnv % IsNeighbour(i+1) ) THEN
           CALL MPI_BSEND( cnt(i), 1, MPI_INTEGER, &
               i, 7001, ELMER_COMM_WORLD, status, ierr )
           IF ( cnt(i)>0 ) THEN
             CALL MPI_BSEND( Buf(i) % grow, cnt(i), MPI_INTEGER, &
                 i, 7002, ELMER_COMM_WORLD, status, ierr )

             CALL MPI_BSEND( Buf(i) % gcol, cnt(i), MPI_INTEGER, &
                 i, 7003, ELMER_COMM_WORLD, status, ierr )

             CALL MPI_BSEND( Buf(i) % gval, cnt(i), MPI_DOUBLE_PRECISION, &
                 i, 7004, ELMER_COMM_WORLD, status, ierr )
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
           MPI_ANY_SOURCE, 7001, ELMER_COMM_WORLD, status, ierr )

         IF ( rcnt>0 ) THEN
           proc = status(MPI_SOURCE)
           CALL MPI_RECV( rrow, rcnt, MPI_INTEGER, &
              proc, 7002, ELMER_COMM_WORLD, status, ierr )

           CALL MPI_RECV( rcol, rcnt, MPI_INTEGER, &
              proc, 7003, ELMER_COMM_WORLD, status, ierr )

           CALL MPI_RECV( rval, rcnt, MPI_DOUBLE_PRECISION, &
              proc, 7004, ELMER_COMM_WORLD, status, ierr )

           DO j=1,rcnt
             l = SearchNode(A % ParallelInfo,rcol(j),Order=A % ParallelInfo % Gorder )
             IF ( l>0 ) THEN
               k = SearchNode(A % ParallelInfo,rrow(j),Order=A % ParallelInfo % Gorder )
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
     ILUValues = 0._dp
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

!-------------------------------------------------------------------------------
!> Assumes partitionwise smallish invertible "addmatrix"
!-------------------------------------------------------------------------------
    SUBROUTINE CircuitPrec(u,v,ipar)
!-------------------------------------------------------------------------------
      USE DefUtils

      INTEGER :: ipar(*)
      REAL(KIND=dp) u(*), v(*)
!-------------------------------------------------------------------------------
      TYPE(Matrix_t), POINTER :: A
      INTEGER :: i,j
      LOGICAL :: stat

      INTEGER :: ndim, n
      TYPE(Solver_t), POINTER :: sv => Null()
!-------------------------------------------------------------------------------
      A => GlobalMatrix

      ndim = ipar(3)
      u(1:ndim) = v(1:ndim)
      CALL CRS_LUPrecondition( u,v, ipar)

      n = A % CircuitMatrix % NumberOfRows
      IF(n>0) THEN
        IF ( .NOT.ASSOCIATED(sv) ) THEN
          ALLOCATE(sv)
          CALL ListAddString(  sv % Values, 'Linear System Direct Method', 'Umfpack')
          CALL ListAddLogical( sv % Values, 'Linear System Refactorize', .FALSE.)
          CALL ListAddLogical( sv % Values, 'Linear System Free Factorization', .FALSE.)
        END IF
 
        i = ndim - A % ExtraDOFs + 1
        j = ndim - A % ExtraDOFs + n
        CALL Umfpack_SolveSystem( sv, A % CircuitMatrix, u(i:j), v(i:j) )
      END IF
!-------------------------------------------------------------------------------
    END SUBROUTINE CircuitPrec
!-------------------------------------------------------------------------------
 
!-------------------------------------------------------------------------------
    SUBROUTINE CircuitPrecComplex(u,v,ipar)
!-------------------------------------------------------------------------------
      USE DefUtils

      INTEGER :: ipar(*)
      COMPLEX(KIND=dp) u(*), v(*)
!-------------------------------------------------------------------------------
      TYPE(Matrix_t), POINTER :: A
      LOGICAL :: stat
      INTEGER :: i,j,k,l

      REAL(KIND=dp), ALLOCATABLE, SAVE :: ru(:), rv(:)
      INTEGER :: ndim, n
      TYPE(Solver_t), POINTER :: sv => Null()
!-------------------------------------------------------------------------------
      A => GlobalMatrix

      ndim = ipar(3)*2
      u(1:ipar(3)) = v(1:ipar(3))
      CALL CRS_ComplexLUPrecondition( u,v, ipar)

      n = A % CircuitMatrix % NumberOfRows
      IF(n>0) THEN
        IF ( .NOT.ASSOCIATED(sv) ) THEN
          ALLOCATE(sv)
          CALL ListAddString(  sv % Values, 'Linear System Direct Method', 'Umfpack')
          CALL ListAddLogical( sv % Values, 'Linear System Refactorize', .FALSE.)
          CALL ListAddLogical( sv % Values, 'Linear System Free Factorization', .FALSE.)
        END IF
 
        IF(.NOT.ALLOCATED(ru)) THEN
          ALLOCATE(ru(n), rv(n))
        ELSE IF(SIZE(ru)<n) THEN
          DEALLOCATE(ru, rv)
          ALLOCATE(ru(n), rv(n))
        END IF
 
        i = (ndim  - A % ExtraDOFs)/2
        j = 0
        DO k=1,n,2
          j = j + 1
          rv(k)   =  REAL(v(i+j)); rv(k+1) = AIMAG(v(i+j))
        END DO
 
        CALL Umfpack_SolveSystem( sv, A % CircuitMatrix, ru, rv )

        j = 0
        DO k=1,n,2
          j = j + 1
          u(i+j) =  CMPLX( ru(k), ru(k+1), KIND=dp )
        END DO
      END IF
!-------------------------------------------------------------------------------
    END SUBROUTINE CircuitPrecComplex
!-------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Create the Vanka preconditioning.
!------------------------------------------------------------------------------
  SUBROUTINE CircuitPrecCreate(A,Solver)
     USE DefUtils
!------------------------------------------------------------------------------
     TYPE(Matrix_t), TARGET :: A
     TYPE(Solver_t) :: Solver
!------------------------------------------------------------------------------
     INTEGER, POINTER :: Diag(:), Rows(:), Cols(:)
     REAL(KIND=dp), ALLOCATABLE :: TotValues(:)
     LOGICAL ::  found
     INTEGER :: status(MPI_STATUS_SIZE)
     REAL(KIND=dp), ALLOCATABLE, TARGET :: rval(:)
     INTEGER, ALLOCATABLE :: cnt(:), rrow(:),rcol(:), perm(:)
     INTEGER :: i,j,k,l,m,ii,jj,proc,rcnt,nn, dof, dofs, Active, n, nm, ierr,totcnt

     TYPE Buf_t
        REAL(KIND=dp), ALLOCATABLE :: gval(:)
        INTEGER, ALLOCATABLE :: grow(:),gcol(:)
     END TYPE Buf_t
     TYPE(Buf_t), POINTER :: buf(:)

     COMPLEX(KIND=dp) :: c

     TYPE(Matrix_t), POINTER :: tm

     Diag => A % Diag
     Rows => A % Rows
     Cols => A % Cols

     nm = A % NumberOfRows - A % ExtraDOFs
     n  = A % ParallelDOFs
     
     m = SIZE(A % Values)
     ALLOCATE(TotValues(m))

     DO i=1,m
       TotValues(i)=A % Values(i)
     END DO
     IF (ParEnv  % PEs>1 ) THEN
       ALLOCATE(cnt(0:ParEnv % PEs-1))
       cnt = 0
       DO i=nm+1,nm+n
         DO j=Rows(i),Rows(i+1)-1
           IF(Cols(j)<=nm .OR. Cols(j)>nm+n) CYCLE
           IF(TotValues(j)==0) CYCLE

           IF ( A % ParallelInfo % Interface(Cols(j)) ) THEN
             m = A % ParallelInfo % NeighbourList(Cols(j)) % Neighbours(1)
             IF ( m==ParEnv % myPE ) CYCLE
             cnt(m) = cnt(m)+1
           END IF
         END DO
       END DO

       ALLOCATE( buf(0:ParEnv % PEs-1) )
       DO i=0,ParEnv % PEs-1
         IF ( cnt(i) > 0 ) &
           ALLOCATE( Buf(i) % gval(cnt(i)), Buf(i) % grow(cnt(i)), Buf(i) % gcol(cnt(i)) )
       END DO

       cnt = 0
       DO i=nm+1,nm+n
         DO j=Rows(i),Rows(i+1)-1
           IF(Cols(j)<=nm .OR. Cols(j)>nm+n) CYCLE
           IF(TotValues(j)==0) CYCLE

           IF ( A % ParallelInfo % Interface(Cols(j)) ) THEN
             m = A % ParallelInfo % NeighbourList(Cols(j)) % Neighbours(1)
             IF ( m==ParEnv % myPE ) CYCLE
             cnt(m) = cnt(m)+1
             Buf(m) % gcol(cnt(m)) = A % ParallelInfo % GlobalDOFs(Cols(j))
             Buf(m) % gval(cnt(m)) = TotValues(j)
             Buf(m) % grow(cnt(m)) = A % ParallelInfo % GlobalDOFs(i)
           END IF
         END DO
       END DO

       totcnt = SUM(cnt)
       CALL CheckBuffer( ParEnv % PEs*(1+MPI_BSEND_OVERHEAD) + 4*totcnt + &
                  3*COUNT(cnt/=0)*MPI_BSEND_OVERHEAD)

       DO i=0,ParEnv % PEs-1
         IF ( ParEnv % IsNeighbour(i+1) ) THEN
           CALL MPI_BSEND( cnt(i), 1, MPI_INTEGER, i, 7001, ELMER_COMM_WORLD, status, ierr )
           IF ( cnt(i)>0 ) THEN
             CALL MPI_BSEND( Buf(i) % grow, cnt(i), MPI_INTEGER, &
                 i, 7002, ELMER_COMM_WORLD, status, ierr )

             CALL MPI_BSEND( Buf(i) % gcol, cnt(i), MPI_INTEGER, &
                 i, 7003, ELMER_COMM_WORLD, status, ierr )

             CALL MPI_BSEND( Buf(i) % gval, cnt(i), MPI_DOUBLE_PRECISION, &
                 i, 7004, ELMER_COMM_WORLD, status, ierr )
           END IF
         END IF
       END DO

       DO i=0,ParEnv % PEs-1
         IF ( cnt(i)>0 ) &
           DEALLOCATE( Buf(i) % gval, Buf(i) % grow, Buf(i) % gcol )
       END DO
       DEALLOCATE( cnt,Buf )

       DO i=1,ParEnv % NumOfNeighbours
         CALL MPI_RECV( rcnt, 1, MPI_INTEGER, &
           MPI_ANY_SOURCE, 7001, ELMER_COMM_WORLD, status, ierr )

         IF ( rcnt>0 ) THEN
           IF(.NOT.ALLOCATED(rrow)) THEN
             ALLOCATE( rrow(rcnt), rcol(rcnt), rval(rcnt) )
           ELSE IF(SIZE(rrow)<rcnt) THEN
             DEALLOCATE(rrow,rcol,rval)
             ALLOCATE( rrow(rcnt), rcol(rcnt), rval(rcnt) )
           ENDIF

           proc = status(MPI_SOURCE)
           CALL MPI_RECV( rrow, rcnt, MPI_INTEGER, &
              proc, 7002, ELMER_COMM_WORLD, status, ierr )

           CALL MPI_RECV( rcol, rcnt, MPI_INTEGER, &
              proc, 7003, ELMER_COMM_WORLD, status, ierr )

           CALL MPI_RECV( rval, rcnt, MPI_DOUBLE_PRECISION, &
              proc, 7004, ELMER_COMM_WORLD, status, ierr )

           DO j=1,rcnt
             l = SearchNode(A % ParallelInfo,rcol(j),Order=A % ParallelInfo % Gorder )
             IF ( l>0 ) THEN
               k = SearchNode(A % ParallelInfo,rrow(j),Order=A % ParallelInfo % Gorder )
               IF ( k>0 ) THEN
                 IF ( l>=k ) THEN
                   DO m=Diag(k),Rows(k+1)-1
                     IF ( Cols(m) == l ) THEN
                       TotValues(m) = TotValues(m) + rval(j)
                       EXIT
                     ELSE IF( Cols(m)>l) THEN
                       EXIT
                     END IF
                   END DO
                 ELSE
                   DO m=Rows(k),Diag(k)-1
                     IF ( Cols(m)==l ) THEN
                       TotValues(m) = TotValues(m) + rval(j)
                       EXIT
                     ELSE IF( Cols(m)>l) THEN
                       EXIT
                     END IF
                   END DO
                 END IF
               END IF
             END IF
           END DO
         END IF
       END DO
     END IF

     IF(ParEnv % PEs<=1) THEN
       tm => A % CircuitMatrix
     ELSE
       tm => A % ParMatrix % SplittedMatrix % InsideMatrix % CircuitMatrix
     END IF

     IF(ASSOCIATED(tm)) CALL FreeMatrix(tm)

     tm => AllocateMatrix()
     tm % Format = MATRIX_LIST

     IF(ParEnv % PEs<=1) THEN
       A % CircuitMatrix => tm
     ELSE
       A % ParMatrix % SplittedMatrix % InsideMatrix % CircuitMatrix => tm
     END IF
    
     ALLOCATE(Perm(n)); Perm=0

     IF ( A % Complex ) THEN
       j = 0; k = 0
       DO i=nm+1,nm+n,2
         j = j + 1
         IF ( ParEnv % PEs > 1) THEN
           IF ( A % ParallelInfo % NeighbourList(i) % Neighbours(1) /= ParEnv % MyPE ) CYCLE
         END IF
         k = k + 1
         Perm(j) = k
       END DO

       DO i=nm+1,nm+n,2
         IF ( ParEnv % PEs > 1) THEN
           IF ( A % ParallelInfo % NeighbourList(i) % Neighbours(1) /= ParEnv % MyPE ) CYCLE
         END IF
         ii = 2*(Perm((i-nm-1)/2 + 1) - 1)

         DO j=A % Rows(i+1)-2,A % Rows(i),-2
           k = A % Cols(j) - nm
           IF(k <= 0)  EXIT
           IF(k >  n) CYCLE
           jj = 2*(Perm((k-1) / 2 + 1) - 1)
           c = CMPLX( TotValues(j), -TotValues(j+1), KIND=dp )

           IF(ABS(c)>AEPS) THEN
             CALL AddToMatrixElement( tm, ii+1, jj+1,  TotValues(j))
             CALL AddToMatrixElement( tm, ii+1, jj+2, -TotValues(j+1))
             CALL AddToMatrixElement( tm, ii+2, jj+1,  TotValues(j+1))
             CALL AddToMatrixElement( tm, ii+2, jj+2,  TotValues(j))
           END IF
         END DO
       END DO
     ELSE
       ii = 0
       DO i=1,n
         IF ( ParEnv % PEs > 1) THEN
           IF ( A % ParallelInfo % NeighbourList(i+nm) % Neighbours(1) /= ParEnv % MyPE ) CYCLE
         END IF
         ii = ii + 1
         jj = 0
         DO j=A % Rows(i+1)-1,A % Rows(i)
           k = A % Cols(j) - nm
           IF(k <= 0) EXIT
           IF(k  > n) CYCLE
           jj = jj + 1
           IF(ABS(TotValues(j))>AEPS) &
             CALL AddToMatrixElement( tm, ii, jj,  TotValues(j))
         END DO
       END DO
     END IF

     CALL List_ToCRSMatrix(tm)
!------------------------------------------------------------------------------
  END SUBROUTINE CircuitPrecCreate
!------------------------------------------------------------------------------

!> \}

!> \}
