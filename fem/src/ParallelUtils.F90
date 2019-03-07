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
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 02 Apr 2001
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{

!------------------------------------------------------------------------------
!>  Parallel solver utilities for *Solver - routines
!------------------------------------------------------------------------------

MODULE ParallelUtils
     USE SparIterSolve

     IMPLICIT NONE

CONTAINS

#define PARALLEL_FOR_REAL
!-------------------------------------------------------------------------------
  FUNCTION ParallelInit() RESULT ( ParallelEnv )
!-------------------------------------------------------------------------------
    TYPE (ParEnv_t), POINTER :: ParallelEnv
                                                                                                                               
#ifdef PARALLEL_FOR_REAL
    ParallelEnv => ParCommInit( )
#else
    ParEnv % MyPE = 0
    ParEnv % PEs  = 1
    ParallelEnv => ParEnv
#endif
!-------------------------------------------------------------------------------
  END FUNCTION ParallelInit
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  SUBROUTINE ParallelFinalize()
!-------------------------------------------------------------------------------
#ifdef PARALLEL_FOR_REAL
    CALL ParEnvFinalize()
#endif
!-------------------------------------------------------------------------------
  END SUBROUTINE ParallelFinalize
!-------------------------------------------------------------------------------
                                                                                                                               
!--------------------------------'-----------------------------------------------
    SUBROUTINE ParallelInitMatrix( Solver, Matrix, inPerm )
!-------------------------------------------------------------------------------
       INTEGER, POINTER, OPTIONAL :: inPerm(:)
       TYPE(Solver_t) :: Solver
       TYPE(Matrix_t), POINTER :: Matrix
!-------------------------------------------------------------------------------
       TYPE(ParallelInfo_t), POINTER :: MatrixPI, MeshPI
       INTEGER :: i, j, k, l, m, n, DOFs, PDOFs
       LOGICAL :: DGSolver, Found, GB, Global_dof, LocalConstraints, DiscontBC, &
                     OwnersGiven, NeighboursGiven, DGReduced
       TYPE(Mesh_t), POINTER :: Mesh
       TYPE(Element_t), POINTER :: Element
       TYPE(NeighbourList_t), POINTER :: MtrxN, MeshN

       INTEGER :: maxnode, maxedge, maxface, fdofs, maxfdofs, col, &
                  edofs, maxedofs, maxbdofs, l_beg, g_beg, proc, sz

       INTEGER, POINTER :: Perm(:)

       INTEGER, ALLOCATABLE :: Owneddofs(:), Owneddofs2(:)
       INTEGER :: pardofs_bef, pardofs_all
       INTEGER, ALLOCATABLE :: pardofs_in(:), pardofs_out(:), ind(:)
       REAL(KIND=dp) :: realtime,tt

       LOGICAL, ALLOCATABLE :: Done(:)
       TYPE(Matrix_t), POINTER :: G
       INTEGER :: Ierr, status(MPI_STATUS_SIZE)

!-------------------------------------------------------------------------------


!tt = realtime()
#ifdef PARALLEL_FOR_REAL
       IF ( ParEnv % PEs <= 1 .OR. .NOT. ASSOCIATED(Matrix) ) RETURN
       Mesh => Solver % Mesh
       DOFs = Solver % Variable % DOFs


       Perm => Solver % Variable % Perm
       IF(PRESENT(inPerm)) Perm=>InPerm

       n = SIZE(Perm)
       k = n*DOFs + Matrix % ExtraDOFs
       ALLOCATE( Matrix % Perm(k), Matrix % InvPerm(k))

       Matrix % Perm = 0
       DO i=1,n
         IF ( Perm(i) /= 0 )  THEN
            DO j=1,DOFs
               Matrix % Perm((i-1)*DOFs+j) = DOFs * (Perm(i)-1) + j
            END DO
         END IF
       END DO

        DO i=n*DOFs+1,SIZE(Matrix % Perm)
          Matrix % Perm(i) = i
        END DO

       Matrix % INVPerm = 0
       DO i=1,SIZE(Matrix % Perm)
          IF ( Matrix % Perm(i) /= 0 ) THEN
             Matrix % INVPerm(Matrix % Perm(i)) = i
          END IF
       END DO

       IF ( .NOT. Matrix % DGMatrix ) THEN
         n = Matrix % NumberOfRows
         ALLOCATE( Matrix % ParallelInfo )
         ALLOCATE( Matrix % ParallelInfo % NeighbourList(n) )
         CALL AllocateVector( Matrix % ParallelInfo % Interface, n)
         CALL AllocateVector( Matrix % ParallelInfo % GlobalDOFs, n)

         DO i=1,n
           Matrix % ParallelInfo % Interface(i) = .FALSE.
           Matrix % ParallelInfo % GlobalDOFs(i) = 0
           Matrix % ParallelInfo % NeighbourList(i) % Neighbours => NULL()
         END DO

         DO i=1,Mesh % NumberOfNodes
           DO j=1,DOFs
              k = Matrix % Perm((i-1)*DOFs+j)
              IF(k<=0) CYCLE
              Matrix % ParallelInfo % GlobalDOFs(k) = &
                DOFs*(Mesh % ParallelInfo % GlobalDOFs(i)-1)+j
              Matrix % ParallelInfo % Interface(k) = &
                Mesh % ParallelInfo % Interface(i)
              ALLOCATE( Matrix % ParallelInfo % NeighbourList(k) % Neighbours(SIZE( &
                   Mesh % ParallelInfo % NeighbourList(i) % Neighbours)) )
              Matrix % ParallelInfo % NeighbourList(k) % Neighbours = &
                Mesh % ParallelInfo % NeighbourList(i) % Neighbours
           END DO
         END DO

         GB = ListGetLogical( Solver % Values, 'Bubbles in Global System', Found )
         IF (.NOT.Found) GB = .TRUE.

         maxnode = MAXVAL(Mesh % ParallelInfo % GlobalDOFs)
         maxnode = NINT(ParallelReduction(1._dp*maxnode,2))

         edofs = 0; fdofs = 0; maxedofs = 0; maxfdofs = 0
         maxedge = 0; maxface = 0

         IF ( ASSOCIATED(Mesh % Edges) ) THEN
           g_beg = maxnode
           l_beg = Mesh % NumberOfNodes

           n = Mesh % NumberOfEdges

           edofs = Mesh % MaxEdgeDOFS
           maxedofs = NINT(ParallelReduction(edofs*1._dp,2))

           maxedge = 0
           DO i=1,n
             maxedge = MAX(maxedge, Mesh % Edges(i) % GElementindex)
           END DO
           maxedge = NINT(ParallelReduction(1._dp*maxedge,2))

           DO i=1,n
             Element => Mesh % Edges(i)
             DO j=1,Element % BDOFs
               DO m=1,DOFs
                 l = DOFs*(l_beg + edofs*(i-1)+j-1)+m
                 l=Matrix % Perm(l)
                 IF(l==0) CYCLE
                 Matrix % ParallelInfo % GlobalDOFs(l) = &
                     DOFs*(g_beg+maxedofs*(Element % GelementIndex-1)+j-1)+m
                 Matrix % ParallelInfo % Interface(l) = Mesh % ParallelInfo % EdgeInterface(i)
                
                 ALLOCATE(Matrix % Parallelinfo % NeighbourList(l) % Neighbours(SIZE( &
                      Mesh % ParallelInfo % EdgeNeighbourList(i) % Neighbours)))
                 Matrix % Parallelinfo % NeighbourList(l) % Neighbours = &
                      Mesh % ParallelInfo % EdgeNeighbourList(i) % Neighbours
               END DO
             END DO
           END DO
         END IF

         IF ( ASSOCIATED(Mesh % Faces) ) THEN
           g_beg = maxnode + maxedofs*maxedge
           l_beg = Mesh % NumberOfNodes + &
                   Mesh % NumberOfEdges*Mesh % MaxEdgeDOFs

           n = Mesh % NumberOfFaces

           fdofs = Mesh % MaxFaceDOFS
           maxfdofs = NINT(ParallelReduction(fdofs*1._dp,2))

           maxface = 0
           DO i=1,n
             maxface = MAX(maxface, Mesh % Faces(i) % GElementindex)
           END DO
           maxface = NINT(ParallelReduction(1._dp*maxface,2))

           DO i=1,n
             Element => Mesh % Faces(i)
             DO j=1,Element % BDOFs
               DO m=1,DOFs
                 l = Matrix % Perm(DOFs*(l_beg + fdofs*(i-1)+j-1)+m)
                 IF(l==0) CYCLE
                 Matrix % ParallelInfo % GlobalDOFs(l) = &
                     DOFs*(g_beg+maxfdofs*(Element % GelementIndex-1)+j-1)+m
                 Matrix % ParallelInfo % Interface(l) = Mesh % ParallelInfo % FaceInterface(i)

                 ALLOCATE(Matrix % Parallelinfo % NeighbourList(l) % Neighbours(SIZE( &
                      Mesh % ParallelInfo % FaceNeighbourList(i) % Neighbours)))

                 Matrix % Parallelinfo % NeighbourList(l) % Neighbours = &
                      Mesh % ParallelInfo % FaceNeighbourList(i) % Neighbours
               END DO
             END DO
           END DO
         END IF

         IF ( GB ) THEN
           l_beg = Mesh % NumberOfNodes + &
                   Mesh % NumberOfEdges*Mesh % MaxEdgeDOFs + &
                   Mesh % NumberOfFaces*Mesh % MaxFaceDOFs

           g_beg = Maxnode +  maxedge*maxedofs + maxface*maxfdofs
           maxbdofs = NINT(ParallelReduction(1._dp*Mesh % MaxBDOFs,2))

           DO i=1,Mesh % NumberOfBulkElements
             Element=>Mesh % Elements(i)
             DO l=1,Element % BDOFs
               DO j=1,DOFs 
                 k = Matrix % Perm(DOFs*(l_beg+Element % BubbleIndexes(l)-1)+j)
                 IF(k==0) CYCLE
                 Matrix % ParallelInfo % GlobalDOFs(k) = &
                   DOFs*(g_beg+maxbdofs*(Element % GElementIndex-1)+l-1)+j
                 Matrix % ParallelInfo % Interface(k) = .FALSE.
                 ALLOCATE(Matrix % ParallelInfo % NeighbourList(k) % Neighbours(1))
                 Matrix % ParallelInfo % NeighbourList(k) % Neighbours=ParEnv % MyPE
               END DO
             END DO
           END DO
         END IF

         ! Add extra degrees of freedom to parallel structures. The additional
         ! variables are assingned to task zero, and are assumed to be shared by
         ! all tasks (TODO: to be optimized if need be...)
         ! --------------------------------------------------------------------
         g_beg = NINT(ParallelReduction(1._dp*MAXVAL(Matrix % ParallelInfo % GlobalDOFs),2))
         pardofs_all = NINT(ParallelReduction(Matrix % ParallelDOFs*1._dp,2))

         LocalConstraints = ListGetLogical(Solver % Values, 'Partition Local Constraints',Found)
         DiscontBC  = ListGetLogicalAnyBC(CurrentModel,'Discontinuous Boundary' )

         Found = ASSOCIATED(Solver % Matrix % ConstraintMatrix)
         IF(Found) Found = ASSOCIATED(Solver % Matrix % ConstraintMatrix % InvPerm)

         G => Null()
         IF( LocalConstraints  ) THEN
           G => AllocateMatrix()
           G % Format = MATRIX_LIST
         END IF

         OwnersGiven = Matrix % ParallelDOFs>0
         IF(OwnersGiven) THEN
           OwnersGiven = ASSOCIATED(Solver % Matrix % AddMatrix)
           IF(OwnersGiven) OwnersGiven = ALLOCATED(Solver % Matrix % AddMatrix % RowOwner)
         END IF

         NeighboursGiven = Matrix % ParallelDOFs>0
         IF(NeighboursGiven) NeighboursGiven = ASSOCIATED(Solver % Matrix % AddMatrix)
         IF(NeighboursGiven) NeighboursGiven = ASSOCIATED(Solver % Matrix % AddMatrix % ParallelInfo)
         IF(NeighboursGiven) NeighboursGiven = ASSOCIATED(Solver % Matrix % AddMatrix % ParallelInfo % &
                             NeighbourList)

         j=0
         DO i=Matrix % NumberOFRows-Matrix % ExtraDOFs+1,Matrix % NumberOfRows
           j=j+1
           Global_dof = j <= Matrix % ParallelDOFs
           IF (Global_dof) THEN
             Matrix % Parallelinfo % Interface(i) = .TRUE.

             IF(NeighboursGiven) THEN
               ALLOCATE(Matrix % ParallelInfo % NeighbourList(i) % Neighbours( &
                  SIZE(Solver % Matrix % AddMatrix % ParallelInfo % NeighbourList(i) % Neighbours)))
               Matrix % ParallelInfo % NeighbourList(i) % Neighbours = &
                  Solver % Matrix % AddMatrix % ParallelInfo % NeighbourList(i) % Neighbours
             ELSE IF (OwnersGiven) THEN
               ALLOCATE(Matrix % ParallelInfo % NeighbourList(i) % Neighbours(ParEnv % PEs))
               DO k=1,ParEnv % PEs
                 Matrix % ParallelInfo % NeighbourList(i) % Neighbours(k)=k-1
               END DO

               l = j + Solver % Matrix % NumberOfRows
               l = Solver % Matrix % AddMatrix % RowOwner(l)

               Matrix % ParallelInfo % NeighbourList(i) % Neighbours(1)   = l
               Matrix % ParallelInfo % NeighbourList(i) % Neighbours(l+1) = 0
             ELSE
               ALLOCATE(Matrix % ParallelInfo % NeighbourList(i) % Neighbours(ParEnv % PEs))
               l = 0
               DO k=ParEnv % PEs,1,-1
                 l = l + 1
                 Matrix % ParallelInfo % NeighbourList(i) % Neighbours(l)=k-1
               END DO
             END IF

             k = j + g_beg
           ELSE IF(Found) THEN
             n = Solver % Matrix % ConstraintMatrix % InvPerm(j-Matrix % ParallelDOFs)
             IF( n<=0 ) CYCLE
             k = MOD(n-1,Solver % Matrix % NumberOfRows)+1
             n = (n-1) / Solver % Matrix % NumberOfRows

             Matrix % Parallelinfo % Interface(i) = Matrix % ParallelInfo % Interface(k)

             l = SIZE(Matrix % ParallelInfo % NeighbourList(k) % Neighbours)
             ALLOCATE(Matrix % ParallelInfo % NeighbourList(i) % Neighbours(l))

             IF(LocalConstraints) THEN
               IF(.NOT.DiscontBC) THEN
                 IF (Matrix % Parallelinfo % Interface(k)) &
                   CALL Sort(l,Matrix % ParallelInfo % NeighbourList(k) % Neighbours)
               END IF

               DO m=Matrix % Rows(i), Matrix % Diag(i)-1
                 col = Matrix % Cols(m)
                 IF (Matrix % Parallelinfo % Interface(col)) THEN
                   DO l=1,SIZE(Matrix % ParallelInfo % NeighbourList(col) % Neighbours)
                     proc = Matrix % ParallelInfo % NeighbourList(col) % Neighbours(l)
                     IF(proc==ParEnv % mype) THEN
                       CALL List_AddToMatrixElement( G % ListMatrix, proc+1, col, 1._dp )
                     ELSE
                       CALL List_AddToMatrixElement( G % ListMatrix, &
                            proc+1, Matrix % ParallelInfo % GlobalDOFs(col), 1._dp)
                     END IF
                   END DO
                 END IF
               END DO
             END IF

             Matrix % Parallelinfo % NeighbourList(i) % Neighbours = &
                Matrix % ParallelInfo % NeighbourList(k) % Neighbours

             k = Matrix % ParallelInfo % GlobalDOFs(k) + (n+1)*g_beg + pardofs_all
           ELSE
             ALLOCATE(Matrix % ParallelInfo % NeighbourList(i) % Neighbours(1))
             Matrix % Parallelinfo % NeighbourList(i) % Neighbours(1)=ParEnv % myPE
             k = j + g_beg + pardofs_all
           END IF
           Matrix % ParallelInfo % GlobalDofs(i) = k
         END DO

         l=Matrix % ExtraDOFs
         n=SIZE(Perm)
         IF(l>0) THEN
            ALLOCATE(Ind(l))
            k = 0
            DO i=Matrix % NumberOfRows-l+1,Matrix % NumberOfRows
              k = k + 1
              Matrix % Perm(n*DOFs+k) = i
              Ind(k) = Matrix % Parallelinfo % Globaldofs(i)
            END DO
            CALL SortI(l,Ind,Matrix % Perm(n*DOFs+1:))
            DEALLOCATE(Ind)
         END IF

         IF (ASSOCIATED(ParEnv % IsNeighbour).AND. LocalConstraints) THEN
           CALL List_ToCRSMatrix(G)
           ALLOCATE(Done(Matrix % NumberOfROws)); Done=.FALSE.

           DO i=1,ParEnv % PEs
             IF(i-1==ParEnv % Mype) THEN
               IF(G % NumberOfRows>=i) THEN
                 DO j=G % Rows(i), G % Rows(i+1)-1
                   m = G % Cols(j)
                   IF(.NOT.Done(m)) THEN
                     Done(m) = .TRUE.
                     k = SIZE(Matrix % ParallelInfo % Neighbourlist(m) % Neighbours)
                     CALL Sort(k,Matrix % ParallelInfo % NeighbourList(m) % Neighbours)
                   END IF
                 END DO
               END IF
               CYCLE
             END IF
  
             j = 0
             IF(G % NumberOfRows>=i) j = G % Rows(i+1)-G % Rows(i)
  
             IF(ParEnv % IsNeighbour(i)) THEN
               CALL MPI_BSEND( j,1,MPI_INTEGER,i-1,130,Matrix % Comm,ierr )
               IF (j>0) THEN
                 CALL MPI_BSEND(G % Cols(G % Rows(i):),j,MPI_INTEGER,i-1,131,Matrix % Comm,ierr )
               END IF
             ELSE
               IF(j>0) stop 'barf'
             END IF
           END DO

           DO i=1,ParEnv % NumOfNeighbours
             CALL MPI_RECV( j,1,MPI_INTEGER,MPI_ANY_SOURCE,130,Matrix % Comm,status,ierr )
             IF (j>0) THEN
               IF(.NOT.ALLOCATED(owneddofs)) ALLOCATE(owneddofs(j))
               IF(SIZE(owneddofs)<j) THEN
                 DEALLOCATE(owneddofs); ALLOCATE(owneddofs(j))
               END IF

               CALL MPI_RECV(owneddofs,j,MPI_INTEGER,status(MPI_SOURCE),131,Matrix % Comm,status,ierr )
               DO l=1,j
                  m = SearchNode( Matrix % ParallelInfo, Owneddofs(l), Order=Matrix % Perm )
                  IF(m>0) THEN
                    IF(.NOT.Done(m)) THEN
                      Done(m) = .TRUE.
                      k = SIZE(Matrix % ParallelInfo % Neighbourlist(m) % Neighbours)
                      CALL Sort(k,Matrix % ParallelInfo % NeighbourList(m) % Neighbours)
                    END IF
                  END IF
                END DO
             END IF
           END DO
           IF(ALLOCATED(owneddofs)) DEALLOCATE(OwnedDofs)
         END IF
 
         IF(ASSOCIATED(G)) THEN
           IF(ASSOCIATED(G % Cols)) DEALLOCATE(G % Cols)
           IF(ASSOCIATED(G % Rows)) DEALLOCATE(G % Rows)
           IF(ASSOCIATED(G % Diag)) DEALLOCATE(G % Diag)
           IF(ASSOCIATED(G % Values)) DEALLOCATE(G % Values)
           DEALLOCATE(G)
         END IF

       ELSE

         MeshPI => Solver % Mesh % ParallelInfo

         ALLOCATE( Matrix % ParallelInfo )
         MatrixPI => Matrix % ParallelInfo

#if 0
         n = 0
         DO i=1,Mesh % NumberOfBulkElements
           Element => Mesh % Elements(i)
           IF ( .NOT. ASSOCIATED(Element % DGIndexes) ) CYCLE
           n = MAX(n,MAXVAL(Element % DGIndexes))
         END DO
#else
         n = MAXVAL(Matrix % Perm)
#endif

         ALLOCATE( MatrixPI % GlobalDOFs(n) ); MatrixPI % GlobalDOFs=0

         maxnode = MAXVAL(Mesh % ParallelInfo % GlobalDOFs)
         maxnode = NINT(ParallelReduction(1._dp*maxnode,2))

         DGReduced = ListGetLogical(Solver % Values, 'DG Reduced Basis', Found )

         IF( DGReduced ) THEN
           BLOCK 
             INTEGER, POINTER :: DgMap(:), DgMaster(:), DgSlave(:)
             LOGICAL :: GotDgMap, GotMaster, GotSlave
             INTEGER :: group0, group
             LOGICAL, ALLOCATABLE :: Tagged(:)
             
             DgMap => ListGetIntegerArray( Solver % Values,'DG Reduced Basis Mapping',GotDgMap )
             DgMaster => ListGetIntegerArray( Solver % Values,'DG Reduced Basis Master Bodies',GotMaster )
             DgSlave => ListGetIntegerArray( Solver % Values,'DG Reduced Basis Slave Bodies',GotSlave )
             
             IF( GotSlave .AND. GotMaster ) THEN
               DO group0 = 1, 2

                 DO i=1,Mesh % NumberOfBulkElements
                   Element => Mesh % Elements(i)
                   IF ( .NOT. ASSOCIATED(Element % DGIndexes) ) CYCLE                 
                   group = Element % BodyId

                   IF( group0 == 1 ) THEN
                     IF( .NOT. ANY( DgMaster == group ) ) CYCLE
                   ELSE
                     IF( ANY ( DgMaster == group ) ) CYCLE                     
                   END IF
                   group = group0 - 1 

                   DO j=1,Element % TYPE % NumberOfNodes
                     k = Matrix % Perm(Element % DGIndexes(j))
                     IF(k == 0) CYCLE                                                       
                     
                     ! Set the global index for slave dofs only if it not already set for
                     ! local dofs
                     IF( group0 == 2 ) THEN
                       IF( MatrixPI % GlobalDOFs(k) > 0 ) CYCLE
                     END IF
                     
                     MatrixPI % GlobalDOFs(k) = group * maxnode +  &
                         MeshPI % GlobalDOFs(Element % NodeIndexes(j))               
                   END DO
                 END DO
               END DO
             ELSE               
               DO i=1,Mesh % NumberOfBulkElements
                 Element => Mesh % Elements(i)
                 IF ( .NOT. ASSOCIATED(Element % DGIndexes) ) CYCLE

                 group0 = Element % BodyId

                 IF( GotMaster ) THEN
                   IF( ANY( DgMaster == group0 ) ) THEN
                     group = 1
                   ELSE
                     group = 2
                   END IF
                 ELSE IF( GotDgMap ) THEN
                   group = DgMap( group0 )
                 ELSE
                   group = group0                   
                 END IF
                 group = group - 1

                 DO j=1,Element % TYPE % NumberOfNodes
                   k = Matrix % Perm(Element % DGIndexes(j))
                   IF(k == 0) CYCLE                                                       
                   MatrixPI % GlobalDOFs(k) = group * maxnode +  &
                       MeshPI % GlobalDOFs(Element % NodeIndexes(j))               
                 END DO
               END DO
             END IF             
           END BLOCK
         ELSE         
           DO i=1,Mesh % NumberOfBulkElements
             Element => Mesh % Elements(i)
             IF ( .NOT. ASSOCIATED(Element % DGIndexes) ) CYCLE
             DO j=1,Element % TYPE % NumberOfNodes
               k = Matrix % Perm(Element % DGIndexes(j))
               IF(k == 0) CYCLE
               MatrixPI % GlobalDOFs(k) = 8*(Element % GElementIndex-1) + j
             END DO
           END DO
         END IF
         
         ALLOCATE( MatrixPI % INTERFACE(n), MatrixPI % NeighbourList(n) )
         MatrixPI % Interface = .FALSE.
         DO i=1,n
           MtrxN => MatrixPI % NeighbourList(i)
           MtrxN % Neighbours => NULL()
         END DO

         DO i=1,Mesh % NumberOfBulkElements
           Element => Mesh % Elements(i)
           IF ( .NOT. ASSOCIATED(Element % DGIndexes) ) CYCLE
           IF ( ALL(MeshPI % Interface(Element % NodeIndexes)) ) THEN
             DO j=1,Element % Type % NumberOfNodes
                K = Matrix % Perm(Element % DGIndexes(j))
                IF(K==0) CYCLE
                L = Element % Nodeindexes(j)

                MeshN => MeshPI % NeighbourList(L)
                MtrxN => MatrixPI % NeighbourList(K)

                MatrixPI % Interface(k) = .TRUE.

                CALL AllocateVector( MtrxN % Neighbours,  SIZE(MeshN % Neighbours) )
                MtrxN % Neighbours = MeshN % Neighbours
                IF(.NOT.DGReduced) THEN ! ? 
                  DO m=1,SIZE(MeshN % Neighbours)
                   IF ( MeshN % Neighbours(m) == Element % PartIndex ) THEN
                     MtrxN % Neighbours(1) = MeshN % Neighbours(m)
                     MtrxN % Neighbours(m) = MeshN % Neighbours(1)
                     EXIT
                   END IF
                  END DO
                END IF
             END DO
           END IF
         END DO
         DO i=1,n
           MtrxN => MatrixPI % NeighbourList(i)
           IF ( .NOT.ASSOCIATED( MtrxN % Neighbours) ) THEN
             CALL AllocateVector(MtrxN % Neighbours,1)
             MtrxN % Neighbours(1) = ParEnv % myPE
           END IF
         END DO
       END IF

       n = SIZE(Matrix % ParallelInfo % GLobalDOFs)
       ALLOCATE(Matrix % ParallelInfo % Gorder(n), Ind(n))

       Ind = Matrix % ParallelInfo % GlobalDOFs
       DO i=1,n
         Matrix % ParallelInfo % Gorder(i) = i
       END DO
       CALL SortI( n,Ind, Matrix % ParallelInfo % Gorder )

       Matrix % ParMatrix => &
          ParInitMatrix( Matrix, Matrix % ParallelInfo )
!if(parenv%mype==0) print*,'MATRIX INIT TIME: ', realtime()-tt
#endif
!-------------------------------------------------------------------------------
    END SUBROUTINE ParallelInitMatrix
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    SUBROUTINE ParallelInitSolve( Matrix, x, b, r, Update )
!-------------------------------------------------------------------------------
       REAL(KIND=dp) CONTIG :: x(:), b(:), r(:)
       LOGICAL, OPTIONAL :: Update
       TYPE(Matrix_t), POINTER :: Matrix
!-------------------------------------------------------------------------------
       LOGICAL :: Upd
#ifdef PARALLEL_FOR_REAL
       Upd = .TRUE.
       IF ( PRESENT(Update) ) Upd=Update
       CALL SParInitSolve( Matrix, x, b, r, Matrix % ParallelInfo, Upd )
#endif
!-------------------------------------------------------------------------------
    END SUBROUTINE ParallelInitSolve
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
    SUBROUTINE ParallelSumVector( Matrix, x, Op )
!-------------------------------------------------------------------------------
       TYPE(Matrix_t) :: Matrix
       INTEGER, OPTIONAL :: op
       REAL(KIND=dp) CONTIG :: x(:)
!-------------------------------------------------------------------------------
       ParEnv = Matrix % ParMatrix % ParEnv
       ParEnv % ActiveComm = Matrix % Comm

       CALL ExchangeSourceVec( Matrix, Matrix % ParMatrix % SplittedMatrix, &
              Matrix % ParallelInfo, x, op )
!-------------------------------------------------------------------------------
    END SUBROUTINE ParallelSumVector
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
    SUBROUTINE ParallelSumNodalVector( Mesh, x, Perm, Matrix, Op )
!-------------------------------------------------------------------------------
      TYPE(Mesh_t) :: Mesh
      INTEGER, POINTER :: Perm(:)
      REAL(KIND=dp) CONTIG :: x(:)
      TYPE(Matrix_t), OPTIONAL :: Matrix
      INTEGER, OPTIONAL :: op
!-------------------------------------------------------------------------------

      ! We can inherit the ParEnv from the primary matrix even
      ! though the variable is not directly associated to it!
      IF( PRESENT( Matrix ) ) THEN
        ParEnv = Matrix % ParMatrix % ParEnv
        ParEnv % ActiveComm = Matrix % Comm
      END IF

      CALL Info('ParallelSumNodalVector','Summing up parallel nodal vector',Level=12)
      
      CALL ExchangeNodalVec( Mesh % ParallelInfo, Perm, x, op )

      CALL Info('ParallelSumNodalVector','Summing up done',Level=20)
!-------------------------------------------------------------------------------
    END SUBROUTINE ParallelSumNodalVector
!-------------------------------------------------------------------------------

    
!-------------------------------------------------------------------------------
    SUBROUTINE ParallelUpdateSolve( Matrix, x, r )
!-------------------------------------------------------------------------------
       REAL(KIND=dp) CONTIG :: x(:), r(:)
       TYPE(Matrix_t), POINTER :: Matrix
!-------------------------------------------------------------------------------
#ifdef PARALLEL_FOR_REAL
       CALL SParUpdateSolve( Matrix, x, r )
#endif
!-------------------------------------------------------------------------------
    END SUBROUTINE ParallelUpdateSolve
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    SUBROUTINE ParallelMatrixVector( Matrix, x, b, Update, UseMassVals,ZeroNotOwned )
!-------------------------------------------------------------------------------
      REAL(KIND=dp) CONTIG :: x(:), b(:)
      TYPE(Matrix_t), POINTER :: Matrix
      LOGICAL, OPTIONAL :: Update, UseMassVals,ZeroNotOwned
!-------------------------------------------------------------------------------
      INTEGER :: i,ipar(1)
      REAL(KIND=dp), POINTER CONTIG :: Mx(:), Mr(:), Mb(:), r(:)

      TYPE(Matrix_t), POINTER :: SaveMatrix
      TYPE(SplittedMatrixT), POINTER :: SP
      TYPE(Matrix_t), POINTER :: SavePtrIN
      TYPE(BasicMatrix_t), POINTER :: SavePtrIF(:), SavePtrNB(:)
!-------------------------------------------------------------------------------
#ifdef PARALLEL_FOR_REAL
      GlobalData => Matrix % ParMatrix
      SaveMatrix  => GlobalMatrix
      GlobalMatrix => Matrix
      ParEnv = GlobalData % ParEnv
      ParEnv % ActiveComm = Matrix % Comm
      IF ( PRESENT( Update ) ) THEN
         IF ( Update ) THEN
            IF ( PRESENT(UseMassVals) ) THEN
               IF ( UseMassVals ) THEN

                  SP => GlobalData % SplittedMatrix
                  ALLOCATE( SavePtrIF( ParEnv % PEs ) )
                  ALLOCATE( SavePtrNB( ParEnv % PEs ) )
                  ALLOCATE( SavePtrIn )
                  DO i=1,ParEnv % PEs
                    IF ( SP % IfMatrix(i) % NumberOfRows /= 0 ) THEN
                       ALLOCATE(SavePtrIF(i) % Values(SIZE(SP % IfMatrix(i) % Values)))
                       SavePtrIF(i) % Values = SP % IfMatrix(i) % Values
                    END IF
                    IF ( SP % NbsIfMatrix(i) % NumberOfRows /= 0 ) THEN
                       ALLOCATE(SavePtrNB(i) % Values(SIZE(SP % NbsIfMatrix(i) % Values)))
                       SavePtrNB(i) % Values = SP % NbsIfMatrix(i) % Values
                    END IF
                  END DO
                  SavePtrIN % Values => SP % InsideMatrix % Values

                  DO i=1,ParEnv % PEs
                    IF ( SP % IfMatrix(i) % NumberOfRows /= 0 ) &
                      SP % IfMatrix(i) % Values = SP % IfMatrix(i) % MassValues

                    IF ( SP % NbsIfMatrix(i) % NumberOfRows /= 0 ) &
                       SP % NbsIfMatrix(i) % Values = SP % NbsIfMatrix(i) % MassValues 
                  END DO
                  SP % InsideMatrix % Values => SP % InsideMatrix % MassValues
               END IF

               Mx => GlobalData % SplittedMatrix % TmpXVec 
               Mr => GlobalData % SplittedMatrix % TmpRVec 
               CALL SParMatrixVector( Mx, Mr, ipar )
               CALL SParUpdateResult( Matrix, x, b, .FALSE. )

               IF(PRESENT(ZeroNotOwned)) THEN
                 IF(ZeroNotOwned) THEN
                   DO i = 1, Matrix % NumberOFRows
                     IF ( Matrix % ParallelInfo % NeighbourList(i) % Neighbours(1) /= ParEnv % MyPE ) THEN
                       b(i) = 0._dp
                     END IF
                   END DO
                 END IF
               ENDIF
            ELSE
               Mx => GlobalData % SplittedMatrix % TmpXVec 
               Mr => GlobalData % SplittedMatrix % TmpRVec 
               CALL SParMatrixVector( Mx, Mr, ipar )
               CALL SParUpdateResult( Matrix, x, b, .FALSE. )

               IF(PRESENT(ZeroNotOwned)) THEN
                 IF(ZeroNotOwned) THEN
                   DO i = 1, Matrix % NumberOFRows
                     IF ( Matrix % ParallelInfo % NeighbourList(i) % Neighbours(1) /= ParEnv % MyPE ) THEN
                       b(i) = 0._dp
                     END IF
                   END DO
                 END IF
               ENDIF
            END IF

            IF ( PRESENT(UseMassVals) ) THEN
               IF ( UseMassVals ) THEN
                  DO i=1,ParEnv % PEs
                    IF ( SP % IfMatrix(i) % NumberOfRows /= 0 ) THEN
                       ALLOCATE(SP % IfMatrix(i) % Values(SIZE(SavePtrIF(i) % Values)))
                       SP % IfMatrix(i) % Values =  SavePtrIF(i) % Values
                    END IF

                    IF ( SP % NbsIfMatrix(i) % NumberOfRows /= 0 ) THEN
                       ALLOCATE(SP % NbsIfMatrix(i) % Values(SIZE(SavePtrNB(i) % Values)))
                       SP % NbsIfMatrix(i) % Values =  SavePtrNB(i) % Values 
                    END IF
                  END DO
                  SP % InsideMatrix % Values => SavePtrIN % Values 
                  DEALLOCATE( SavePtrIF )
                  DEALLOCATE( SavePtrNB )
                  DEALLOCATE( SavePtrIn )
               END IF
            END IF
         ELSE
            IF ( PRESENT(UseMassVals) ) THEN
               IF ( UseMassVals ) THEN
                  SP => Matrix % ParMatrix % SplittedMatrix
                  ALLOCATE( SavePtrIF( ParEnv % PEs ) )
                  ALLOCATE( SavePtrNB( ParEnv % PEs ) )
                  ALLOCATE( SavePtrIn )
                  DO i=1,ParEnv % PEs
                    IF ( SP % IfMatrix(i) % NumberOfRows /= 0 ) THEN
                       ALLOCATE(SavePtrIF(i) % Values(SIZE(SP % IfMatrix(i) % Values)))
                       SavePtrIF(i) % Values = SP % IfMatrix(i) % Values
                    END IF

                    IF ( SP % NbsIfMatrix(i) % NumberOfRows /= 0 ) THEN
                       ALLOCATE(SavePtrNB(i) % Values(SIZE(SP % NbsIfMatrix(i) % Values)))
                       SavePtrNB(i) % Values = SP % NbsIfMatrix(i) % Values
                    END IF
                  END DO
                  SavePtrIN % Values => SP % InsideMatrix % Values

                  DO i=1,ParEnv % PEs
                    IF ( SP % IfMatrix(i) % NumberOfRows /= 0 ) &
                       SP % IfMatrix(i) % Values = SP % IfMatrix(i) % MassValues

                    IF ( SP % NbsIfMatrix(i) % NumberOfRows /= 0 ) &
                       SP % NbsIfMatrix(i) % Values = SP % NbsIfMatrix(i) % MassValues 
                  END DO
                  SP % InsideMatrix % Values => SP % InsideMatrix % MassValues
               END IF
            END IF

            CALL SParMatrixVector( x, b, ipar )

            IF ( PRESENT(UseMassVals) ) THEN
              IF ( UseMassVals ) THEN
                DO i=1,ParEnv % PEs
                  IF ( SP % IfMatrix(i) % NumberOfRows /= 0 ) THEN
                     SP % IfMatrix(i) % Values = SavePtrIF(i) % Values
                     DEALLOCATE(SavePtrIF(i) % Values)
                  END IF

                  IF ( SP % NbsIfMatrix(i) % NumberOfRows /= 0 ) THEN
                     SP % NbsIfMatrix(i) % Values = SavePtrNB(i) % Values 
                     DEALLOCATE(SavePtrNB(i) % Values)
                  END IF
                END DO
                SP % InsideMatrix % Values => SavePtrIN % Values 
                DEALLOCATE( SavePtrIF )
                DEALLOCATE( SavePtrNB )
                DEALLOCATE( SavePtrIN )
              END IF
            END IF
         END IF
      ELSE
         CALL SParMatrixVector( x, b, ipar )
      END IF
      GlobalMatrix => SaveMatrix
#endif
!-------------------------------------------------------------------------------
    END SUBROUTINE ParallelMatrixVector
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
    SUBROUTINE ParallelMatrixVectorC( Matrix, x, b, Update, UseMassVals,ZeroNotOwned )
!-------------------------------------------------------------------------------
      COMPLEX(KIND=dp) CONTIG :: x(:), b(:)
      TYPE(Matrix_t), POINTER :: Matrix
      LOGICAL, OPTIONAL :: Update, UseMassVals,ZeroNotOwned
!-------------------------------------------------------------------------------
      INTEGER :: i,ipar(1)
      REAL(KIND=dp), POINTER CONTIG :: Mx(:), Mr(:), Mb(:), r(:)

      TYPE(Matrix_t), POINTER :: SaveMatrix
      TYPE(SplittedMatrixT), POINTER :: SP
      TYPE(Matrix_t), POINTER :: SavePtrIN
      TYPE(BasicMatrix_t), POINTER :: SavePtrIF(:), SavePtrNB(:)
!-------------------------------------------------------------------------------
#ifdef PARALLEL_FOR_REAL
      GlobalData => Matrix % ParMatrix
      SaveMatrix  => GlobalMatrix
      GlobalMatrix => Matrix
      ParEnv = GlobalData % ParEnv
      ParEnv % ActiveComm = Matrix % Comm
      IF ( PRESENT( Update ) ) THEN
        CALL Fatal('ParallelMatrixVectorC','Cannot handle parameter > Update <')
      END IF
        
      IF ( PRESENT( UseMassVals ) ) THEN
        CALL Fatal('ParallelMatrixVectorC','Cannot handle parameter > UseMassVals <')
      END IF

      CALL SParCMatrixVector( x, b, ipar )

      GlobalMatrix => SaveMatrix
#endif
!-------------------------------------------------------------------------------
    END SUBROUTINE ParallelMatrixVectorC
!-------------------------------------------------------------------------------

    
!-------------------------------------------------------------------------------
    SUBROUTINE ParallelVectorC(A, vec_out, vec_in)
!-------------------------------------------------------------------------------
      TYPE(Matrix_t), INTENT(in) :: A
      COMPLEX(KIND=dp), INTENT(inout) :: vec_out(:)
      COMPLEX(KIND=dp), INTENT(in), OPTIONAL :: vec_in(:)
!-------------------------------------------------------------------------------
      INTEGER :: i,j,k
!-------------------------------------------------------------------------------
      j = 0

      ! We have a complex valued vector but a real valued matrix.
      ! We use the even (complex) component to check the ownership of the dof.
      ! We could as well use the odd (real) component.
  
      DO i=1,A % NumberOfRows / 2
        IF ( A % ParallelInfo % Neighbourlist(2*i) % &
                   Neighbours(1)==Parenv % Mype ) THEN
          j=j+1
          IF(PRESENT(vec_in)) THEN
            vec_out(j) = vec_in(i)
          ELSE
            vec_out(j) = vec_out(i)
          END IF
        END IF
      END DO
!-------------------------------------------------------------------------------
    END SUBROUTINE ParallelVectorC
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
    SUBROUTINE ParallelVector(A, vec_out, vec_in)
!-------------------------------------------------------------------------------
      TYPE(Matrix_t), INTENT(in) :: A
      REAL(KIND=dp), INTENT(inout) :: vec_out(:)
      REAL(KIND=dp), INTENT(in), OPTIONAL :: vec_in(:)
!-------------------------------------------------------------------------------
      INTEGER :: i,j,k
!-------------------------------------------------------------------------------
      j = 0
      DO i=1,A % NumberOfRows
        IF ( A % ParallelInfo % Neighbourlist(i) % &
                   Neighbours(1)==Parenv % Mype ) THEN
          j=j+1
          IF(PRESENT(vec_in)) THEN
            vec_out(j) = vec_in(i)
          ELSE
            vec_out(j) = vec_out(i)
          END IF
        END IF
      END DO
!-------------------------------------------------------------------------------
    END SUBROUTINE ParallelVector
!-------------------------------------------------------------------------------

    
!-------------------------------------------------------------------------------
    SUBROUTINE PartitionVector(A, vec_out, vec_in)
!-------------------------------------------------------------------------------
      TYPE(Matrix_t), INTENT(in) :: A
      REAL(KIND=dp), INTENT(inout) :: vec_in(:), vec_out(:)
!-------------------------------------------------------------------------------
      INTEGER :: i,j,k
!-------------------------------------------------------------------------------
      j = 0
      vec_out = 0._dp
      DO i=1,A % NumberOfRows
        IF ( A % ParallelInfo % Neighbourlist(i) % &
                   Neighbours(1)==Parenv % Mype ) THEN
          j=j+1
          vec_out(i) = vec_in(j)
        END IF
      END DO
!-------------------------------------------------------------------------------
    END SUBROUTINE PartitionVector
!-------------------------------------------------------------------------------
    

!-------------------------------------------------------------------------------
    SUBROUTINE ParallelUpdateResult( Matrix, x, r )
!-------------------------------------------------------------------------------
       REAL(KIND=dp) :: x(:), r(:)
       TYPE(Matrix_t), POINTER :: Matrix
!-------------------------------------------------------------------------------
#ifdef PARALLEL_FOR_REAL
       CALL SParUpdateResult( Matrix, x, r, .TRUE. )
#endif
!-------------------------------------------------------------------------------
    END SUBROUTINE ParallelUpdateResult
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    SUBROUTINE ParallelUpdateRHS( Matrix, b )
!-------------------------------------------------------------------------------
       REAL(KIND=dp) :: b(:)
       TYPE(Matrix_t), POINTER :: Matrix
!-------------------------------------------------------------------------------
#ifdef PARALLEL_FOR_REAL
       CALL SParUpdateRHS( Matrix, b, Matrix % ParallelInfo )
#endif
!-------------------------------------------------------------------------------
    END SUBROUTINE ParallelUpdateRHS
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------

    FUNCTION ParallelMatrix( A,x,b,r ) RESULT(M)
!-------------------------------------------------------------------------------
       TYPE(Matrix_t), POINTER :: A, M
       REAL(KIND=dp),  POINTER, OPTIONAL :: x(:),b(:),r(:)
!-------------------------------------------------------------------------------
       M => NULL()
#ifdef PARALLEL_FOR_REAL
       M => A % ParMatrix % SplittedMatrix % InsideMatrix
       IF ( PRESENT(x) ) x => A % ParMatrix % SplittedMatrix % TmpXVec
       IF ( PRESENT(b) ) b => M % RHS
       IF ( PRESENT(r) ) r => A % ParMatrix % SplittedMatrix % TmpRVec
#endif
!-------------------------------------------------------------------------------
    END FUNCTION ParallelMatrix
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    FUNCTION ParallelNorm( n, x ) RESULT(s)
!-------------------------------------------------------------------------------
      INTEGER :: n
      REAL(KIND=dp) :: s
      REAL(KIND=dp) CONTIG :: x(:)
!-------------------------------------------------------------------------------
      s = 0.0d0
#ifdef PARALLEL_FOR_REAL
      s = SParNorm( n, x, 1 )
#endif
!-------------------------------------------------------------------------------
    END FUNCTION ParallelNorm
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
    FUNCTION ParallelCNorm( n, x ) RESULT(s)
!-------------------------------------------------------------------------------
      INTEGER :: n
      REAL(KIND=dp) :: s
      COMPLEX(KIND=dp) CONTIG :: x(:)
!-------------------------------------------------------------------------------
      s = 0.0d0
#ifdef PARALLEL_FOR_REAL
      s = SParCNorm( n, x, 1 )
#endif
!-------------------------------------------------------------------------------
    END FUNCTION ParallelCNorm
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
    FUNCTION ParallelDOT( n, x, y ) RESULT(s)
!-------------------------------------------------------------------------------
      INTEGER :: n
      REAL(KIND=dp) :: s
      REAL(KIND=dp) CONTIG :: x(:),y(:)
!-------------------------------------------------------------------------------
      s = 0.0d0
#ifdef PARALLEL_FOR_REAL
      s = SParDotProd( n, x, 1, y, 1 )
#endif
!-------------------------------------------------------------------------------
    END FUNCTION ParallelDOT
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    SUBROUTINE ParallelGlobalNumbering(Mesh,OldMesh,NewNodes,IntCnts,IntArray,Reorder)
!-------------------------------------------------------------------------------
       TYPE(Mesh_t) :: Mesh, OldMesh
       INTEGER :: NewNodes,IntCnts(:),IntArray(:),Reorder(:)
!-------------------------------------------------------------------------------
#ifdef PARALLEL_FOR_REAL
        CALL SparGlobalNumbering( Mesh,OldMesh,NewNodes,IntCnts,IntArray,Reorder )
#endif
!-------------------------------------------------------------------------------
    END SUBROUTINE ParallelGlobalNumbering
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
    SUBROUTINE ParallelIter( SourceMatrix, ParallelInfo, DOFs, XVec, &
              RHSVec, Solver, SParMatrixDesc )
!-------------------------------------------------------------------------------
       TYPE (Matrix_t) :: SourceMatrix
       TYPE (ParallelInfo_t) :: ParallelInfo
       INTEGER :: DOFs
       REAL(KIND=dp), DIMENSION(:) :: XVec, RHSVec
       TYPE (Solver_t) :: Solver
       TYPE (SParIterSolverGlobalD_t), POINTER :: SParMatrixDesc
                                                                                                                               
#ifdef PARALLEL_FOR_REAL
       CALL SParIterSolver( SourceMatrix, ParallelInfo, XVec, &
                 RHSVec, Solver, SParMatrixDesc )
#endif
!-------------------------------------------------------------------------------
    END SUBROUTINE ParallelIter
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    SUBROUTINE ParallelActive( L )
!-------------------------------------------------------------------------------
      LOGICAL :: L
#ifdef PARALLEL_FOR_REAL
       IF ( ParEnv % PEs > 1) THEN
         CALL SParIterBarrier
         CALL SParIterActive(L)
       END IF
#endif
!-------------------------------------------------------------------------------
    END SUBROUTINE ParallelActive
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    SUBROUTINE ParallelAllReduceAnd( L )
!-------------------------------------------------------------------------------
      LOGICAL :: L
#ifdef PARALLEL_FOR_REAL
       IF ( ParEnv % PEs > 1) CALL SParIterAllReduceAnd(L)
#endif
!-------------------------------------------------------------------------------
    END SUBROUTINE ParallelAllReduceAnd
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    FUNCTION ParallelReduction(R,oper_arg) RESULT(rsum)
!-------------------------------------------------------------------------------
      REAL(KIND=dp) :: R, rsum
      INTEGER, OPTIONAL :: oper_arg
!-------------------------------------------------------------------------------
      INTEGER :: oper
!-------------------------------------------------------------------------------
      rsum = r
#ifdef PARALLEL_FOR_REAL
      IF ( ParEnv % PEs>1) THEN
        oper = 0
        IF (PRESENT(oper_arg)) THEN
          oper=oper_arg
        ELSE
          oper = 0
        END IF
        IF(.NOT.ASSOCIATED(ParEnv % Active)) &
          CALL ParallelActive(.TRUE.)
        CALL SparActiveSUM(rsum,oper)
      END IF
#endif
!-------------------------------------------------------------------------------
    END FUNCTION ParallelReduction
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    SUBROUTINE ParallelBarrier
!-------------------------------------------------------------------------------
#ifdef PARALLEL_FOR_REAL
       IF ( ParEnv % PEs > 1) CALL SParIterBarrier
#endif
!-------------------------------------------------------------------------------
    END SUBROUTINE ParallelBarrier
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    SUBROUTINE ParallelActiveBarrier
!-------------------------------------------------------------------------------
#ifdef PARALLEL_FOR_REAL
       IF ( ParEnv % PEs > 1 ) CALL SParIterActiveBarrier
#endif
!-------------------------------------------------------------------------------
    END SUBROUTINE ParallelActiveBarrier
!-------------------------------------------------------------------------------

END MODULE ParallelUtils

!> \}
