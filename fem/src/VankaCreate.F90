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

#include "../config.h"

!-------------------------------------------------------------------------------
!> Vanka preconditioning for iterative methods.
!-------------------------------------------------------------------------------
    SUBROUTINE VankaPrec(u,v,ipar)
!-------------------------------------------------------------------------------
      USE DefUtils
      IMPLICIT NONE
      
      INTEGER :: ipar(*)
      REAL(KIND=dp) u(*), v(*)
!-------------------------------------------------------------------------------
      TYPE(SplittedMatrixT), POINTER :: SP
      TYPE(Matrix_t), POINTER :: A
      INTEGER :: i
      INTEGER, POINTER :: pCols(:), pRows(:)
      REAL(KIND=dp), POINTER CONTIG :: SaveValues(:)
      TYPE(BasicMatrix_t), POINTER :: SaveIF(:)
!-------------------------------------------------------------------------------
      A => GlobalMatrix
      SaveValues => A % Values
      A % Values => A % ILUValues

      IF (ParEnv % Pes <= 1 ) THEN
        IF( ASSOCIATED( A % ILUCols ) ) THEN
          pCols => A % Cols
          pRows => A % Rows
          A % Cols => A % ILUCols
          A % Rows => A % ILURows
        END IF
        CALL CRS_MatrixVectorProd(v,u,ipar)
        IF( ASSOCIATED( A % ILUCols ) ) THEN
          A % Cols => pCols
          A % Rows => pRows
        END IF
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
    IMPLICIT NONE 
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     TYPE(Solver_t), TARGET :: Solver
!------------------------------------------------------------------------------
     INTEGER, POINTER :: Diag(:), Rows(:), Cols(:), Perm(:), Indexes(:), Ind(:)
     REAL(KIND=dp), POINTER CONTIG :: ILUValues(:), SValues(:), TotValues(:)
     REAL(KIND=dp), ALLOCATABLE :: al(:,:)
     LOGICAL ::  found
     TYPE(Element_t), POINTER :: Element
     INTEGER :: status(MPI_STATUS_SIZE)
     INTEGER :: i,j,i2,j2,ierr,k,l,m,proc,rcnt,nn, dof, dofs, Active, Totcnt
     REAL(KIND=dp), ALLOCATABLE, TARGET :: rval(:)
     INTEGER, ALLOCATABLE :: cnt(:), rrow(:),rcol(:)
     REAL(KIND=dp) :: veps
     
     TYPE Buf_t
        REAL(KIND=dp), ALLOCATABLE :: gval(:)
        INTEGER, ALLOCATABLE :: grow(:),gcol(:)
     END TYPE Buf_t
     TYPE(Buf_t), POINTER :: buf(:)

     TYPE(Matrix_t), POINTER :: B
     INTEGER :: VankaMode
     
     Diag => A % Diag
     Rows => A % Rows
     Cols => A % Cols

     Dofs  = Solver % Variable % DOFs
     Perm => Solver % Variable % Perm

     Svalues => A % Values
     
     ! If we have a block solver then we may be coming here with another dof count!
     dofs = dofs * A % NumberOfRows / SIZE(Solver % Variable % Values) 
     
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
           IF ( A % ParallelInfo % GInterface(Cols(j)) ) THEN
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
               i, 7001, ELMER_COMM_WORLD, ierr )
           IF ( cnt(i)>0 ) THEN
             CALL MPI_BSEND( Buf(i) % grow, cnt(i), MPI_INTEGER, &
                 i, 7002, ELMER_COMM_WORLD, ierr )

             CALL MPI_BSEND( Buf(i) % gcol, cnt(i), MPI_INTEGER, &
                 i, 7003, ELMER_COMM_WORLD, ierr )

             CALL MPI_BSEND( Buf(i) % gval, cnt(i), MPI_DOUBLE_PRECISION, &
                 i, 7004, ELMER_COMM_WORLD, ierr )
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

     nn = 10*Solver % Mesh % MaxElementDOFs
     ALLOCATE( Indexes(nn), AL(nn*dofs,nn*dofs), ind(nn*dofs) )

     VankaMode = ListGetInteger(Solver % Values,'Vanka Mode',Found)
     veps = ListGetCReal( Solver % Values,'Vanka epsilon',Found)
     IF(.NOT. Found) veps = 1.0e-6
     
     IF( VankaMode == 0 ) THEN
       ! For the basic mode the filling is exactly the same as for the primary matrix
       IF ( .NOT. ASSOCIATED(A % ILUValues) ) THEN
         ALLOCATE( A % ILUvalues(SIZE(A % Values)) )
       END IF
       ILUValues => A % ILUValues
       ILUValues = 0.0_dp
     ELSE
       ! For the other modes we don't know the fill pattern and hence use a list matrix for assembly
       B => AllocateMatrix()
       B % FORMAT = MATRIX_LIST
       A % Values => TotValues
     END IF
       
     
     SELECT CASE(VankaMode)

     CASE(0)

       ! Pick entries related to a single element and inverse the matrix.
       ! Add the inverse to the preconditioning matrix.
       !-------------------------------------------------------------------
       Active = GetNOFActive(Solver)
       DO i=1,Active
         element => GetActiveElement(i)
         nn = GetElementDOFs(Indexes)

         l = 0
         DO j=1,nn
           k =  Indexes(j)
           IF(perm(k)==0) CYCLE
           DO dof=1,dofs
             l = l + 1
             ind(l) = DOFs*(perm(k)-1)+dof
           END DO
         END DO
         IF(l==0) CYCLE
         
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

     CASE(1)
       CALL Info('VankaCreate','Using block created from connections in matrix row')
       
       ! Add the max index first because list matrix likes this
       i = A % NumberOfRows 
       CALL List_AddToMatrixElement( B % ListMatrix,i,i,0.0_dp )
       
       DO i=1,A % NumberOfRows / dofs
         j = dofs*(i-1)+1
         
         nn = A % Rows(j+1)-A % Rows(j)
         IF(nn > SIZE(ind)) CALL Fatal('VankaCreate','Index too large for "ind" table!')         
         Ind(1:nn) = A % Cols(A % Rows(j):A % Rows(j+1)-1)

         CALL AssembleVankaBlock()
       END DO
         
     CASE(2) 

       BLOCK
         TYPE(Matrix_t), POINTER :: NodeGraph, EdgeGraph
         TYPE(Element_t), POINTER :: Edge
         TYPE(Mesh_t), POINTER :: Mesh
         INTEGER :: n0

         CALL Info('VankaCreate','Using agressive block created around edge')
         
         Mesh => Solver % Mesh
         IF(.NOT. ASSOCIATED(Mesh % Edges)) THEN
           CALL Fatal('VankaCreate','This version requires Edges!')
         END IF
           
         ! Create a graph for node-to-edge connectivity
         !----------------------------------------------
         NodeGraph => AllocateMatrix()
         NodeGraph % FORMAT = MATRIX_LIST         
         DO i = Mesh % NumberOfEdges, 1, -1
           Edge => Mesh % Edges(i)
           DO j=1, Edge % TYPE % NumberOfNodes 
             CALL List_AddToMatrixElement( NodeGraph % ListMatrix,Edge % NodeIndexes(j),i,1.0_dp )
           END DO
         END DO
         CALL List_ToCRSMatrix(NodeGraph)
         PRINT *,'Nonzeros per row NodeGraph:',1.0_dp * SIZE(NodeGraph % Values) / NodeGraph % NumberOfRows

         ! Create a graph for edge-to-edge connectivity
         !----------------------------------------------
         EdgeGraph => AllocateMatrix()
         EdgeGraph % FORMAT = MATRIX_LIST                 
         DO i = Mesh % NumberOfEdges, 1, -1
           Edge => Mesh % Edges(i)
           DO j=1, Edge % TYPE % NumberOfNodes 
             k = Edge % NodeIndexes(j)
             DO l = NodeGraph % Rows(k),NodeGraph % Rows(k+1)-1
               CALL List_AddToMatrixElement( EdgeGraph % ListMatrix,i,NodeGraph % Cols(l),1.0_dp )
             END DO
           END DO
         END DO
         CALL FreeMatrix( NodeGraph)
         
         CALL List_ToCRSMatrix(EdgeGraph)
         PRINT *,'Nonzeros per row EdgeGraph:',1.0_dp * SIZE(EdgeGraph % Values) / EdgeGraph % NumberOfRows

         ! Create the preconditioning matrix
         !-----------------------------------
         i = A % NumberOfRows 
         CALL List_AddToMatrixElement( B % ListMatrix,i,i,0.0_dp )

         n0 = Mesh % NumberOfNodes

         DO i=1, Mesh % NumberOfEdges
           l = 0

           ! First add the nodes
           DO k=EdgeGraph % Rows(i), EdgeGraph % Rows(i+1)-1
             j = EdgeGraph % Cols(k)
             Edge => Mesh % Edges(j)
             DO i2=1,2
               j2 = Edge % NodeIndexes(i2)
               IF( ALL( Indexes(1:l) /= j2 ) )  THEN
                 l = l+1
                 Indexes(l) = j2
               END IF
             END DO
           END DO

           ! Then add the edges
           DO k=EdgeGraph % Rows(i), EdgeGraph % Rows(i+1)-1
             j = EdgeGraph % Cols(k)
             l = l+1
             indexes(l) = j + n0
           END DO

           nn=l
           l=0
           DO j=1,nn
             k = Indexes(j)
             IF(perm(k)==0) CYCLE
             DO dof=1,dofs
               l = l + 1
               IF(l > SIZE(ind)) CALL Fatal('VankaCreate','Index too large for "ind" table!')
               ind(l) = dofs*(perm(k)-1)+dof
             END DO
           END DO
           IF(l==0) CYCLE

           nn = l
           CALL AssembleVankaBlock()
         END DO
         CALL FreeMatrix( EdgeGraph )         
         
       END BLOCK


     CASE(3)
       BLOCK
         TYPE(Mesh_t), POINTER :: Mesh
         TYPE(Element_t), POINTER :: Face
         INTEGER :: nn2, NoElems
         INTEGER, POINTER :: Indexes2(:)
         
         CALL Info('VankaCreate','Using block created around each face')
         
         Mesh => Solver % Mesh
         IF( Mesh % MeshDim == 3 ) THEN
           IF(.NOT. ASSOCIATED(Mesh % Faces)) THEN
             CALL Warn('VankaCreate','This mode requires existence of Faces in 3D!')
             CALL FindMeshFaces3D(Mesh)
           END IF
           NoElems = Mesh % NumberOfFaces
         ELSE
           IF(.NOT. ASSOCIATED(Mesh % Edges)) THEN
             CALL Warn('VankaCreate','This mode requires existence of Edges in 2D!')
             CALL FindMeshEdges2D(Mesh)
           END IF
           
           NoElems = Mesh % NumberOfEdges
         END IF

         ALLOCATE(Indexes2(SIZE(Indexes)))
         
         ! Add the max index first because list matrix likes this
         i = A % NumberOfRows 
         CALL List_AddToMatrixElement( B % ListMatrix,i,i,0.0_dp )
         
         DO i=1, NoElems 
           IF( Mesh % MeshDim == 3 ) THEN
             Face => Mesh % Faces(i)
           ELSE
             Face => Mesh % Edges(i)
           END IF
             
           nn = 0
           nn2 = 0
           DO j=1,2
             IF(j==1) THEN
               Element => Face % BoundaryInfo % Left
             ELSE
               Element => Face % BoundaryInfo % Right
             END IF
             IF(.NOT. ASSOCIATED(Element) ) CYCLE

             IF(j==1) THEN 
               nn = GetElementDOFs(Indexes,Element)
             ELSE
               nn2 = GetElementDOFs(Indexes2,Element)
             END IF
           END DO

           IF(nn2 > 0 .AND. nn > 0 ) THEN
             DO j=1,nn2
               IF( ALL(Indexes(1:nn) /= Indexes2(j) ) ) THEN
                 nn = nn+1
                 Indexes(nn) = Indexes2(j)
               END IF
             END DO
           ELSE IF(nn2 == 0 .AND. nn == 0 ) THEN
             CALL Fatal('VankaCreate','No dofs in parents?')
           ELSE IF( nn2 == 0 ) THEN
             CONTINUE
           ELSE IF( nn == 0 ) THEN
             nn = nn2
             Indexes(1:nn) = Indexes2(1:nn2)             
           END IF
             
           l = 0
           DO j=1,nn
             k = Indexes(j)
             IF(perm(k)==0) CYCLE
             DO dof=1,dofs
               l = l + 1
               ind(l) = dofs*(perm(k)-1)+dof
             END DO
           END DO
           IF(l==0) CYCLE

           nn = l
           CALL AssembleVankaBlock()
         END DO

       END BLOCK

#if 0
     CASE(4)
       ! This does not really work for the problems tested
       BLOCK                  
         TYPE(Mesh_t), POINTER :: Mesh
         INTEGER :: NoLayers = 0
         INTEGER, POINTER :: DownPointer(:), TopPointer(:)
         TYPE(Variable_t), POINTER :: ExtVar
         TYPE(Solver_t), POINTER :: pSolver
         
         SAVE ExtVar, DownPointer, TopPointer, NoLayers
         
         CALL Info('VankaCreate','Using block created by inverse extrusion') 

         pSolver => Solver
         Mesh => Solver % Mesh
         
         ! Find the extruded structure 
         IF( NoLayers == 0 ) THEN
           CALL DetectExtrudedStructure( Mesh, pSolver, ExtVar, &
               TopNodePointer = TopPointer, DownNodePointer = DownPointer, &
               NumberOfLayers = NoLayers )
         END IF

         i = (NoLayers+1)*dofs
         IF( SIZE(Ind) < i ) THEN
           DEALLOCATE( Ind, al)
           ALLOCATE(Ind(i),al(i,i))           
         END IF
           
         DO i=1,Mesh % NumberOfNodes
           IF( TopPointer(i) == i ) THEN
             l = 1
             k = i
             Indexes(1) = k
             DO j = 1,NoLayers
               k = DownPointer(k)
               l = l+1
               Indexes(l) = k
             END DO
             
             nn = l
             l = 0
             DO j=1,nn
               k = Indexes(j)
               IF(perm(k)==0) CYCLE
               DO dof=1,dofs
                 l = l + 1
                 IF(l > SIZE(ind)) THEN
                   PRINT *,'l:',l,size(ind)
                   CALL Fatal('VankaCreate','Index too large for "ind" table!')
                 END IF
                 ind(l) = dofs*(perm(k)-1)+dof
               END DO
             END DO
             IF(l==0) CYCLE

             nn = l
             CALL AssembleVankaBlock()
           END IF
         END DO                      
       END BLOCK
#endif

     CASE DEFAULT

       CALL Fatal('VankaCreate','Unknown vanka mode: '//I2S(VankaMode))
       
     END SELECT

     ! This is common to all other vanka modes except the basic elemental one. 
     IF( VankaMode > 0 ) THEN
       CALL List_ToCRSMatrix(B)
       PRINT *,'Nonzeros per row Vanka:',1.0_dp * SIZE(B % Values) / B % NumberOfRows
       PRINT *,'Fill ratio for Vanka:',1.0_dp * SIZE(B % Values) / SIZE(A % Values)
       
       IF(ASSOCIATED(A % ILUValues)) DEALLOCATE(A % ILUValues)
       IF(ASSOCIATED(A % ILUCols)) DEALLOCATE(A % ILUCols)
       IF(ASSOCIATED(A % ILURows)) DEALLOCATE(A % ILURows)
       
       A % ILUValues => B % Values
       A % ILUCols => B % Cols
       A % ILURows => B % Rows
       
       ! Nullify these so that they won't be destroyed
       NULLIFY( B % Values, B % Cols, B % Rows)
       CALL FreeMatrix( B )               
     END IF
   
     A % Values => Svalues
     DEALLOCATE(AL, Indexes, Ind, TotValues)


   CONTAINS
     
     SUBROUTINE AssembleVankaBlock()
       
       INTEGER :: jj, kk
       REAL(KIND=dp) :: asum, ab
       
       
       al(1:nn,1:nn) = 0.0_dp
       DO j=1,nn
         DO k=1,nn
           IF( CRS_CheckMatrixElement( A,Ind(j), Ind(k) ) ) THEN
             al(j,k) = CRS_GetMatrixElement( A, ind(j), ind(k) )
           END IF
         END DO
       END DO

       CALL InvertMatrix(al,nn)
       asum = SUM(ABS(al(1:nn,1:nn)))

       ! For complex problems we need to have all matrix entries related to same complex number present
       IF( A % COMPLEX ) THEN
         DO j=1,nn/2
           DO k=1,nn/2
             ab = SUM( ABS(AL(2*j-1:2*j,2*k-1:2*k)) )
             IF(ab < veps * asum ) CYCLE
             DO jj=-1,0
               DO kk=-1,0                     
                 CALL List_AddToMatrixElement( B % ListMatrix,ind(2*j+jj),ind(2*k+kk),AL(2*j+jj,2*k+kk) )
               END DO
             END DO
           END DO
         END DO
       ELSE    
         DO j=1,nn
           DO k=1,nn
             ab = ABS(AL(j,k))
             IF(ab < veps * asum ) CYCLE
             CALL List_AddToMatrixElement( B % ListMatrix,ind(j),ind(k),AL(j,k) )
           END DO
         END DO
       END IF
       
     END SUBROUTINE AssembleVankaBlock
              
!------------------------------------------------------------------------------
  END SUBROUTINE VankaCreate
!------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Assumes partitionwise smallish invertible "addmatrix"
!-------------------------------------------------------------------------------
    SUBROUTINE CircuitPrec(u,v,ipar)
!-------------------------------------------------------------------------------
      USE DefUtils
      !USE DirectSolve, ONLY: MumpsLocal_SolveSystem, Umfpack_SolveSystem
      IMPLICIT NONE
      
      INTEGER :: ipar(*)
      REAL(KIND=dp) u(*), v(*)
!-------------------------------------------------------------------------------
      TYPE(Matrix_t), POINTER :: A
      INTEGER :: i,j,k
      LOGICAL :: stat
      INTEGER :: ndim, n
      CHARACTER(:), ALLOCATABLE, SAVE :: str
      TYPE(Solver_t), POINTER, SAVE :: sv => Null()
!-------------------------------------------------------------------------------
      A => GlobalMatrix
      n = A % CircuitMatrix % NumberOfRows

      ndim = ipar(3)

      CALL CRS_LUPrecondition(u,v, ipar)

      IF(n>0) THEN
        IF ( .NOT.ASSOCIATED(sv) ) THEN
          ALLOCATE(sv)
          str = ListGetString( CurrentModel % Solver % Values, &
              'Linear System Direct Method', Stat )
          IF(.NOT. Stat ) str = "umfpack"                    
#if !defined (HAVE_UMFPACK) && defined (HAVE_MUMPS)
          IF( str == "umfpack" ) THEN
            CALL Warn( 'CircuitPrec', 'Umfpack solver not installed, using MUMPS instead!' )
            str = "mumps"
          END IF
#elseif !defined (HAVE_MUMPS) && defined(HAVE_UMFPACK)
          IF( str == "mumps" ) THEN
            CALL Warn( 'CircuitPrec', 'MUMPS solver not installed, using Umfpack instead!' )
            str = "umfpack"
          END IF
#elseif !defined (HAVE_MUMPS) && !defined(HAVE_UMFPACK)
          CALL Fatal( 'CircuitPrec', 'Preconditioner "circuit" needs either Umfpack or MUMPS!')
#endif
          CALL ListAddString( sv % Values, 'Linear System Direct Method', TRIM(str) )
          CALL ListAddLogical( sv % Values, 'Linear System Refactorize', .FALSE.)
          CALL ListAddLogical( sv % Values, 'Linear System Free Factorization', .FALSE.)

          CALL Info('CircuitPrec','Using direct solver '&
              //TRIM(str)//' of size '//I2S(A % ExtraDofs),Level=10)          
        END IF
        i = ndim - A % ExtraDOFs + 1
        j = ndim - A % ExtraDOFs + n
        
        IF(ANY(ABS(A % CircuitMatrix % Values)>0)) THEN
          SELECT CASE( str )
          CASE('umfpack')
            CALL Umfpack_SolveSystem( sv, A % CircuitMatrix, u(i:j), v(i:j) )
          CASE('mumps')
            CALL MumpsLocal_SolveSystem( sv, A % CircuitMatrix, u(i:j), v(i:j) )
          CASE DEFAULT
            CALL Fatal('CircuitPrec','Impossible direct method: '//TRIM(str))
          END SELECT
        END IF
      END IF

!-------------------------------------------------------------------------------
    END SUBROUTINE CircuitPrec
!-------------------------------------------------------------------------------
 
!-------------------------------------------------------------------------------
    SUBROUTINE CircuitPrecComplex(u,v,ipar)
!-------------------------------------------------------------------------------
      USE DefUtils
      IMPLICIT NONE
      
      INTEGER :: ipar(*)
      COMPLEX(KIND=dp) u(*), v(*)
!-------------------------------------------------------------------------------
      TYPE(Matrix_t), POINTER :: A
      LOGICAL :: stat
      INTEGER :: i,j,k,l
      CHARACTER(:), ALLOCATABLE, SAVE :: str
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
          str = ListGetString( CurrentModel % Solver % Values, &
              'Linear System Direct Method', Stat )
          IF(.NOT. Stat ) str = "umfpack"                    
#if !defined (HAVE_UMFPACK) && defined (HAVE_MUMPS)
          IF( str == "umfpack" ) THEN
            CALL Warn( 'CircuitPrecComplex', 'Umfpack solver not installed, using MUMPS instead!' )
            str = "mumps"
          END IF
#elseif !defined (HAVE_MUMPS) && defined(HAVE_UMFPACK)
          IF( str == "mumps" ) THEN
            CALL Warn( 'CircuitPrecComplex', 'MUMPS solver not installed, using Umfpack instead!' )
            str = "umfpack"
          END IF
#elseif !defined (HAVE_MUMPS) && !defined(HAVE_UMFPACK)
          CALL Fatal( 'CircuitPrecComplex', 'Preconditioner "circuit" needs either Umfpack or MUMPS!')
#endif
          CALL ListAddString( sv % Values, 'Linear System Direct Method', TRIM(str) )
          CALL ListAddLogical( sv % Values, 'Linear System Refactorize', .FALSE.)
          CALL ListAddLogical( sv % Values, 'Linear System Free Factorization', .FALSE.)

          CALL Info('CircuitPrecComplex','Using direct solver '&
              //TRIM(str)//' of size '//I2S(A % ExtraDofs/2),Level=10)          
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

        SELECT CASE( str )
        CASE('umfpack')
          CALL Umfpack_SolveSystem( sv, A % CircuitMatrix, ru, rv )
        CASE('mumps')
          CALL MumpsLocal_SolveSystem( sv, A % CircuitMatrix, ru, rv )
        CASE DEFAULT
          CALL Fatal('CircuitPrecComplex','Impossible direct method: '//TRIM(str))
        END SELECT
        
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
     IMPLICIT NONE
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

           IF ( A % ParallelInfo % GInterface(Cols(j)) ) THEN
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

           IF ( A % ParallelInfo % GInterface(Cols(j)) ) THEN
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
           CALL MPI_BSEND( cnt(i), 1, MPI_INTEGER, i, 7001, ELMER_COMM_WORLD, ierr )
           IF ( cnt(i)>0 ) THEN
             CALL MPI_BSEND( Buf(i) % grow, cnt(i), MPI_INTEGER, &
                 i, 7002, ELMER_COMM_WORLD, ierr )

             CALL MPI_BSEND( Buf(i) % gcol, cnt(i), MPI_INTEGER, &
                 i, 7003, ELMER_COMM_WORLD, ierr )

             CALL MPI_BSEND( Buf(i) % gval, cnt(i), MPI_DOUBLE_PRECISION, &
                 i, 7004, ELMER_COMM_WORLD, ierr )
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
       j = 0; k = 0
       DO i=nm+1,nm+n
         j = j + 1
         IF ( ParEnv % PEs > 1) THEN
           IF ( A % ParallelInfo % NeighbourList(i) % Neighbours(1) /= ParEnv % MyPE ) CYCLE
         END IF
         k = k + 1
         Perm(j) = k
       END DO

       DO i=nm+1,nm+n
         IF ( ParEnv % PEs > 1) THEN
           IF ( A % ParallelInfo % NeighbourList(i) % Neighbours(1) /= ParEnv % MyPE ) CYCLE
         END IF
         ii = Perm(i-nm)
         DO j=A % Rows(i+1)-1,A % Rows(i),-1
           k = A % Cols(j) - nm
           IF(k <= 0) EXIT
           IF(k  > n) CYCLE
           jj = Perm(k)
           CALL AddToMatrixElement(tm, ii, jj,  TotValues(j))
         END DO
       END DO
     END IF

     CALL List_ToCRSMatrix(tm)
!------------------------------------------------------------------------------
  END SUBROUTINE CircuitPrecCreate
!------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!> Assumes another solver being used for the preconditioning.
!> Given residual "v" solver Au=v in an approximatite manner.
!-------------------------------------------------------------------------------
  SUBROUTINE SlavePrec(u,v,ipar)
!-------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    INTEGER :: ipar(*)  ! parameters for Hutiter
    REAL(KIND=dp) u(*)  ! new solution
    REAL(KIND=dp) v(*)  ! right-hand-side
!-------------------------------------------------------------------------------
    TYPE(Solver_t), POINTER :: Solver
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Variable_t), POINTER :: pVar
    TYPE(Matrix_t), POINTER :: Amat
    REAL(KIND=dp), POINTER :: b(:), x(:), r(:)
    REAL(KIND=dp) :: rnorm
    LOGICAL :: Found
    CHARACTER(MAX_NAME_LEN) :: str   
    INTEGER :: n
!-------------------------------------------------------------------------------

    Solver => CurrentModel % Solver
    Params => Solver % Values
    Mesh => Solver % Mesh 
    Amat => Solver % Matrix
    
    str = ListGetString( Params,'Slave Prec Residual',UnfoundFatal=.TRUE.)
    pVar => VariableGet( Mesh % Variables, str )
    IF(.NOT. ASSOCIATED(pVar)) CALL Fatal('SlavePrec','Could not find: '//TRIM(str))
    n = SIZE(pVar % Values)
    b => pVar % Values

    b(1:n) = v(1:n)

    CALL DefaultSlaveSolvers( Solver, 'Prec Solvers' )

    str = ListGetString( Params,'Slave Prec Update',UnfoundFatal=.TRUE.)
    pVar => VariableGet( Mesh % Variables, str )
    IF(.NOT. ASSOCIATED(pVar)) CALL Fatal('SlavePrec','Could not find: '//TRIM(str))
    x => pVar % Values
    
    IF( ListCheckPresent( Params,'MG Smoother') ) THEN
      ALLOCATE(r(n))
      
      IF( ListGetLogical( Params,'MG Smoother Normalize Guess',Found) )  THEN
        BLOCK
          REAL(KIND=dp) :: rn, bn    
          CALL MatrixVectorMultiply( Amat, x, r) 
          rn = SUM( r(1:n)**2 )
          bn = SUM( r(1:n) * b(1:n) )
          IF( rn > TINY( rn ) ) THEN
            bn = bn / rn 
            x(1:n) = x(1:n) * bn 
            WRITE( Message,'(A,ES12.3)') 'Preconditioning Normalizing Factor: ',bn
            CALL Info('SlavePrec',Message,Level=6) 
          END IF
        END BLOCK
      END IF

      CALL CRS_MatrixVectorMultiply( Amat, x, r )
      !CALL MGmv( Amat, x, r, .TRUE. )
      r(1:n) = b(1:n) - r(1:n)
      RNorm = MGSmooth( Solver, Amat, Mesh, x, b, r, &
          1, pVar % dofs, PreSmooth = .FALSE.)
    END IF
          
    u(1:n) = x(1:n) 

!-------------------------------------------------------------------------------
  END SUBROUTINE SlavePrec
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
  SUBROUTINE SlavePrecComplex(u,v,ipar)
!-------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    INTEGER :: ipar(*)  ! parameters for Hutiter
    COMPLEX(KIND=dp) u(*)  ! new solution
    COMPLEX(KIND=dp) v(*)  ! right-hand-side
!-------------------------------------------------------------------------------
    TYPE(Solver_t), POINTER :: Solver
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Variable_t), POINTER :: pVar
    TYPE(Matrix_t), POINTER :: Amat
    REAL(KIND=dp), POINTER :: b(:), x(:), r(:)
    REAL(KIND=dp) :: rnorm
    LOGICAL :: Found
    CHARACTER(MAX_NAME_LEN) :: str   
    INTEGER :: n
!-------------------------------------------------------------------------------

    Solver => CurrentModel % Solver
    Params => Solver % Values
    Mesh => Solver % Mesh
    Amat => Solver % Matrix

    str = ListGetString( Params,'Slave Prec Residual',UnfoundFatal=.TRUE.)
    pVar => VariableGet( Mesh % Variables, str )
    IF(.NOT. ASSOCIATED(pVar)) CALL Fatal('SlavePrecComplex','Could not find: '//TRIM(str))
    n = SIZE(pVar % Values)   
    IF(pVar % Dofs /= Solver % Variable % dofs ) THEN
      CALL Fatal('SlavePrecComplex','Residual should have same size as primary variable!')
    END IF
    IF(n /= SIZE(Solver % Variable % Values) ) THEN
      CALL Fatal('SlavePrecComplex','Residual should have same size as primary variable!')
    END IF
    b => pVar % Values

    b(1:n:2) = REAL(v(1:n/2))
    b(2:n:2) = AIMAG(v(1:n/2))
    
    CALL DefaultSlaveSolvers( Solver, 'Prec Solvers' )
    
    str = ListGetString( Params,'Slave Prec Update',UnfoundFatal=.TRUE.)
    pVar => VariableGet( Mesh % Variables, str )    
    IF(.NOT. ASSOCIATED(pVar)) CALL Fatal('SlavePrec','Could not find: '//TRIM(str))
    IF(pVar % Dofs /= Solver % Variable % dofs ) THEN
      CALL Fatal('SlavePrecComplex','Update should have same size as primary variable!')
    END IF
    IF(n /= SIZE(Solver % Variable % Values) ) THEN
      CALL Fatal('SlavePrecComplex','Update should have same size as primary variable!')
    END IF

    x => pVar % Values
    
    IF( ListCheckPresent( Params,'MG Smoother') ) THEN      
      ALLOCATE(r(n))
      
      IF( ListGetLogical( Params,'MG Smoother Normalize Guess',Found) )  THEN
        BLOCK
          REAL(KIND=dp) :: rn, bnre, bnim    
          CALL MatrixVectorMultiply( Amat, x, r) 
          rn = SUM( r(1:n)**2 )
          bnre = SUM( r(1:n) * b(1:n) )          
          bnim = SUM( r(1:n:2) * b(2:n:2) - r(2:n:2) * b(1:n:2) )
          
          IF( rn > TINY( rn ) ) THEN
            bnre = bnre / rn
            bnim = bnim / rn
#if 0
            ! This does not seem to help ...
            b(1:n) = x(1:n)
            x(1:n:2) = bnre * r(1:n:2) - bnim * r(2:n:2)
            x(2:n:2) = bnim * r(1:n:2) + bnre * r(2:n:2)
#else            
            x(1:n) = x(1:n) * bnre
#endif
            WRITE( Message,'(A,2ES12.3)') 'Preconditioning Normalizing Factor: ',bnre,bnim
            CALL Info('SlavePrec',Message,Level=6) 
          END IF
        END BLOCK
      END IF
        
      CALL CRS_MatrixVectorMultiply( Amat, x, r )
      !CALL MGmv( Amat, x, r, .TRUE. )
      r(1:n) = b(1:n) - r(1:n)
      RNorm = MGSmooth( Solver, Amat, Mesh, x, b, r, &
          1, pVar % dofs, PreSmooth = .FALSE.)
      DEALLOCATE(r)
    END IF
      
    u(1:n/2) = CMPLX(x(1:n:2), x(2:n:2) ) 

!-------------------------------------------------------------------------------
  END SUBROUTINE SlavePrecComplex
!-------------------------------------------------------------------------------

  
  
!> \}

!> \}
