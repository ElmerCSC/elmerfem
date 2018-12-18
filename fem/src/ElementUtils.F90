!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
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
! *  Original Date: 01 Oct 1996
! *
! ******************************************************************************/

!> \ingroup ElmerLib
!> \}

!--------------------------------------------------------------------------------
!>  Some basic finite element utilities.
!--------------------------------------------------------------------------------
MODULE ElementUtils

    USE DirectSolve
    USE Integration
    USE ListMatrix
    USE ListMatrixArray
    USE BandMatrix
    USE Lists
    USE CRSMatrix
    USE Interpolation
    USE BandwidthOptimize

CONTAINS

!------------------------------------------------------------------------------
!> Frees structures of the matrix.
!------------------------------------------------------------------------------
   RECURSIVE SUBROUTINE FreeMatrix(Matrix)
!------------------------------------------------------------------------------
     TYPE(Matrix_t), POINTER :: Matrix
!------------------------------------------------------------------------------
     TYPE(Solver_t) :: Solver
     REAL(KIND=dp) :: x(1), b(1)

     TYPE(SplittedMatrixT), POINTER :: s
     TYPE(BasicMatrix_t), POINTER :: m
     TYPE(SParIterSolverGlobalD_t), POINTER :: p

#ifdef HAVE_HYPRE
    INTERFACE
      !! destroy the data structures (should be called when the matrix has
      !! to be updated and SolveHYPRE1 has to be called again).
      SUBROUTINE SolveHYPRE4(hypreContainer)
        USE, INTRINSIC :: iso_c_binding
        INTEGER(KIND=C_INTPTR_T) :: hypreContainer
      END SUBROUTINE SolveHYPRE4

    END INTERFACE
#endif
#ifdef HAVE_TRILINOS
     INTERFACE
      !! destroy the data structures (should be called when the matrix has
      !! to be updated and SolveTrilinos1 has to be called again).
      SUBROUTINE SolveTrilinos4(triliContainer) BIND(C,name='SolveTrilinos4')
        USE, INTRINSIC :: iso_c_binding
        INTEGER(KIND=C_INTPTR_T) :: triliContainer
      END SUBROUTINE SolveTrilinos4

     END INTERFACE
#endif

     IF ( .NOT. ASSOCIATED( Matrix ) ) RETURN

     CALL DirectSolver( Matrix,x,b,Solver,Free_Fact=.TRUE.)

     IF ( ASSOCIATED( Matrix % Perm ) )        DEALLOCATE( Matrix % Perm )
     IF ( ASSOCIATED( Matrix % InvPerm ) )     DEALLOCATE( Matrix % InvPerm )

     IF ( ASSOCIATED( Matrix % Rows ) ) THEN
        IF ( ASSOCIATED(Matrix % Rows, Matrix % ILURows) .OR. SIZE(Matrix % Rows)==0 )  &
           Matrix % ILURows => Null()
        DEALLOCATE( Matrix % Rows )
     END IF

     IF ( ASSOCIATED( Matrix % Cols ) ) THEN
        IF ( ASSOCIATED(Matrix % Cols, Matrix % ILUCols) .OR. SIZE(Matrix % Cols)==0 ) &
           Matrix % ILUCols => Null()
        DEALLOCATE( Matrix % Cols )
     END IF

     IF ( ASSOCIATED( Matrix % Diag ) ) THEN
        IF ( ASSOCIATED(Matrix % Diag, Matrix % ILUDiag) .OR. SIZE(Matrix % Diag)==0 ) &
          Matrix % ILUDiag => Null()
        DEALLOCATE( Matrix % Diag )
     END IF

     IF ( ASSOCIATED( Matrix % RHS   ) )       DEALLOCATE( Matrix % RHS )
     IF ( ASSOCIATED( Matrix % Force ) )       DEALLOCATE( Matrix % Force )
     IF ( ASSOCIATED( Matrix % RHS_im ) )      DEALLOCATE( Matrix % RHS_im )

     IF ( ALLOCATED( Matrix % ExtraVals ) )    DEALLOCATE( Matrix % ExtraVals )
     IF ( ASSOCIATED( Matrix % Values ) )      DEALLOCATE( Matrix % Values )
     IF ( ASSOCIATED( Matrix % MassValues ) )  DEALLOCATE( Matrix % MassValues )
     IF ( ASSOCIATED( Matrix % DampValues ) )  DEALLOCATE( Matrix % DampValues )
     IF ( ASSOCIATED( Matrix % BulkValues ) )  DEALLOCATE( Matrix % BulkValues )
     IF ( ASSOCIATED( Matrix % BulkRHS   ) )   DEALLOCATE( Matrix % BulkRHS )

     IF ( ASSOCIATED( Matrix % ILUValues ) )   DEALLOCATE( Matrix % ILUValues )
     IF ( ASSOCIATED( Matrix % ILURows ) )     DEALLOCATE( Matrix % ILURows )
     IF ( ASSOCIATED( Matrix % ILUCols ) )     DEALLOCATE( Matrix % ILUCols )
     IF ( ASSOCIATED( Matrix % ILUDiag ) )     DEALLOCATE( Matrix % ILUDiag )

     IF ( ASSOCIATED( Matrix % CRHS   ) )      DEALLOCATE( Matrix % CRHS )
     IF ( ASSOCIATED( Matrix % CForce ) )      DEALLOCATE( Matrix % CForce )

     IF ( ASSOCIATED( Matrix % CValues ) )     DEALLOCATE( Matrix % CValues )
     IF ( ASSOCIATED( Matrix % CILUValues ) )  DEALLOCATE( Matrix % CILUValues )

     IF ( ASSOCIATED(Matrix % CMassValues) )  DEALLOCATE( Matrix % CMassValues )
     IF ( ASSOCIATED(Matrix % CDampValues) )  DEALLOCATE( Matrix % CDampValues )

     IF ( ALLOCATED( Matrix % GRows ) )      DEALLOCATE( Matrix % GRows )
     IF ( ALLOCATED( Matrix % RowOwner ) )   DEALLOCATE( Matrix % RowOwner )
     IF ( ASSOCIATED( Matrix % GOrder) )      DEALLOCATE( Matrix % GOrder )

     CALL FreeMatrix( Matrix % EMatrix )
     CALL FreeMatrix( Matrix % ConstraintMatrix )

     IF(ASSOCIATED(Matrix % ParallelInfo)) THEN
       DEALLOCATE(Matrix % ParallelInfo % GlobalDOFs)
       DEALLOCATE(Matrix % ParallelInfo % INTERFACE)
       DO i=1,SIZE(Matrix % ParallelInfo % NeighbourList)
         DEALLOCATE(Matrix % ParallelInfo % NeighbourList(i) % Neighbours)
       END DO
       DEALLOCATE(Matrix % ParallelInfo % NeighbourList)
       DEALLOCATE(Matrix % ParallelInfo)
     END IF

     p=>Matrix % ParMatrix
     IF(ASSOCIATED(p)) THEN
       s => Matrix % ParMatrix % SplittedMatrix

       s % InsideMatrix % EMatrix => NULL()
       CALL FreeMatrix(s % InsideMatrix)

       DO i=1,SIZE(s % IfMatrix)
         m => s % IfMatrix(i)
         IF(ASSOCIATED(m)) THEN
           IF(ALLOCATED(m % Rows)) DEALLOCATE(m % Rows)
           IF(ALLOCATED(m % Cols)) DEALLOCATE(m % Cols)
           IF(ALLOCATED(m % Diag)) DEALLOCATE(m % Diag)
           IF(ALLOCATED(m % Grows)) DEALLOCATE(m % Grows)
           IF(ALLOCATED(m % Values)) DEALLOCATE(m % Values)
           IF(ALLOCATED(m % RowOwner)) DEALLOCATE(m % RowOwner)
           IF(ALLOCATED(m % ILUValues)) DEALLOCATE(m % ILUValues)
           IF(ALLOCATED(m % MassValues)) DEALLOCATE(m % MassValues)
           IF(ALLOCATED(m % DampValues)) DEALLOCATE(m % DampValues)
         END IF
       END DO
       DEALLOCATE(s % IfMatrix)

       DO i=1,SIZE(s % NbsIfMatrix)
         m => s % NbsIfMatrix(i)
         IF(ASSOCIATED(m)) THEN
           IF(ALLOCATED(m % Rows)) DEALLOCATE(m % Rows)
           IF(ALLOCATED(m % Cols)) DEALLOCATE(m % Cols)
           IF(ALLOCATED(m % Diag)) DEALLOCATE(m % Diag)
           IF(ALLOCATED(m % Grows)) DEALLOCATE(m % Grows)
           IF(ALLOCATED(m % Values)) DEALLOCATE(m % Values)
           IF(ALLOCATED(m % RowOwner)) DEALLOCATE(m % RowOwner)
           IF(ALLOCATED(m % ILUValues)) DEALLOCATE(m % ILUValues)
           IF(ALLOCATED(m % MassValues)) DEALLOCATE(m % MassValues)
           IF(ALLOCATED(m % DampValues)) DEALLOCATE(m % DampValues)
         END IF
       END DO
       DEALLOCATE(s % NbsIfMatrix)

       IF(ASSOCIATED(s % VecIndices)) THEN
         DO i=1,SIZE(s % VecIndices)
           IF(ASSOCIATED(s % Vecindices(i) % RevInd)) DEALLOCATE(s % VecIndices(i) % RevInd)
         END DO
         DEALLOCATE(s % VecIndices)
       END IF

       IF(ASSOCIATED(s % IfVecs)) THEN
         DO i=1,SIZE(s % IfVecs)
           IF(ASSOCIATED(s % IfVecs(i) % IfVec)) DEALLOCATE(s % IfVecs(i) % IfVec)
         END DO
         DEALLOCATE(s % ifVecs)
       END IF

       IF(ASSOCIATED(s % IfORows)) THEN
         DO i=1,SIZE(s % IfORows)
           IF(ASSOCIATED(s % IfORows(i) % IfVec)) DEALLOCATE(s % IfORows(i) % IfVec)
         END DO
         DEALLOCATE(s % ifORows)
       END IF

       IF(ASSOCIATED(s % IfLCols)) THEN
         DO i=1,SIZE(s % IfLCols)
           IF(ASSOCIATED(s % IfLCols(i) % IfVec)) DEALLOCATE(s % IfLCols(i) % IfVec)
         END DO
         DEALLOCATE(s % ifLCols)
       END IF

       IF(ASSOCIATED(s % ResBuf)) THEN
         DO i=1,SIZE(s % ResBuf)
           IF(ALLOCATED(s % ResBuf(i) % ResVal)) DEALLOCATE(s % ResBuf(i) % ResVal)
           IF(ALLOCATED(s % ResBuf(i) % ResInd)) DEALLOCATE(s % ResBuf(i) % ResInd)
         END DO
         DEALLOCATE(s % ResBuf)
       END IF

       IF(ASSOCIATED(s % RHS)) THEN
         DO i=1,SIZE(s % RHS)
           IF(ASSOCIATED(s % RHS(i) % RHSVec)) DEALLOCATE(s % RHS(i) % RHSVec)
           IF(ASSOCIATED(s % RHS(i) % RHSInd)) DEALLOCATE(s % RHS(i) % RHSInd)
         END DO
         DEALLOCATE(s % RHS)
       END IF

       IF (ASSOCIATED(s % Work)) DEALLOCATE(s % Work)
       IF (ASSOCIATED(s % TmpXVec)) DEALLOCATE(s % TmpXVec)
       IF (ASSOCIATED(s % TmpRVec)) DEALLOCATE(s % TmpRVec)

       IF(ASSOCIATED(s % GlueTable)) THEN
         IF(ASSOCIATED(s % GlueTable % Rows)) DEALLOCATE(s % GlueTable % Rows)
         IF(ASSOCIATED(s % GlueTable % Cols)) DEALLOCATE(s % GlueTable % Cols)
         IF(ASSOCIATED(s % GlueTable % Inds)) DEALLOCATE(s % GlueTable % Inds)
         IF(ASSOCIATED(s % GlueTable % RowOwner)) DEALLOCATE(s % GlueTable % RowOwner)
         DEALLOCATE(s % GlueTable)
       END IF
       DEALLOCATE(s)
       DEALLOCATE(p)
     END IF

#ifdef HAVE_HYPRE
     IF (Matrix % Hypre /= 0) THEN
       CALL SolveHypre4(Matrix % Hypre)
     END IF
#endif

#ifdef HAVE_TRILINOS
     IF (Matrix % Trilinos /= 0) THEN
       CALL SolveTrilinos4(Matrix % Trilinos)
     END IF
#endif
     DEALLOCATE( Matrix )
!------------------------------------------------------------------------------
   END SUBROUTINE FreeMatrix
!------------------------------------------------------------------------------


   
!------------------------------------------------------------------------------
!> Create a list matrix given the mesh, the active domains and the elementtype 
!> related to the solver. The list matrix is flexible since it can account 
!> for any entries. Also constraints and periodic BCs may give rise to entries
!> in the list matrix topology.
!------------------------------------------------------------------------------
  SUBROUTINE MakeListMatrix( Model,Solver,Mesh,List,Reorder, &
        LocalNodes,Equation, DGSolver, GlobalBubbles, &
        NodalDofsOnly, ProjectorDofs, CalcNonZeros )
!------------------------------------------------------------------------------
    TYPE(Model_t)  :: Model
    TYPE(Mesh_t)   :: Mesh
    TYPE(Matrix_t) :: Matrix
    TYPE(SOlver_t) :: Solver
    TYPE(ListMatrix_t), POINTER :: List(:)
    INTEGER, OPTIONAL ::Reorder(:)
    LOGICAL, OPTIONAL :: DGSolver
    LOGICAL, OPTIONAL :: GlobalBubbles
    LOGICAL, OPTIONAL :: NodalDofsOnly
    LOGICAL, OPTIONAL :: ProjectorDofs
    CHARACTER(LEN=*), OPTIONAL :: Equation
    LOGICAL, OPTIONAL :: CalcNonZeros
!------------------------------------------------------------------------------
    INTEGER :: t,i,j,k,l,m,k1,k2,n,p,q,e1,e2,f1,f2, EDOFs, FDOFs, BDOFs, This, istat
    INTEGER, ALLOCATABLE :: InvPerm(:), IndirectPairs(:)

    LOGICAL :: Flag, FoundDG, GB, DB, Found, Radiation, DoProjectors, &
        DoNonZeros, DgIndirect

    TYPE(Matrix_t), POINTER :: Projector

    INTEGER :: IndexSize, NumberOfFactors
    INTEGER, ALLOCATABLE :: Indexes(:)

    TYPE(ListMatrixEntry_t), POINTER :: CList, Lptr

    TYPE(Matrix_t),POINTER :: PMatrix
    TYPE(Element_t), POINTER :: Element,Elm, Edge1, Edge2, Face1, Face2, Left, Right
!------------------------------------------------------------------------------

    CALL Info('MakeListMatrix','Creating list matrix',Level=14)

    GB = .FALSE.
    IF ( PRESENT(GlobalBubbles) ) GB = GlobalBubbles

    IF( DgSolver ) THEN
      DB = ListGetLogical( Solver % Values,'DG Reduced Basis',Found) 
    ELSE
      DB = .FALSE.
    END IF
    
    List => List_AllocateMatrix(LocalNodes)

    BDOFs = Mesh % MaxBDOFs
    EDOFs = Mesh % MaxEdgeDOFs
    FDOFs = Mesh % MaxFaceDOFs
    IF( PRESENT( NodalDofsOnly ) ) THEN
      IF( NodalDofsOnly ) THEN
        EDOFS = 0
        FDOFS = 0
      END IF
    END IF

    IF( PRESENT( ProjectorDofs ) ) THEN
      DoProjectors = ProjectorDofs 
    ELSE
      DoProjectors = .TRUE.
    END IF

    IF( EDOFS > 0 .AND. .NOT. ASSOCIATED(Mesh % Edges) ) THEN
      CALL Warn('MakeListMatrix','Edge dofs requested but not edges exist in mesh!')
      EDOFS = 0
    END IF

    IF( FDOFS > 0 .AND. .NOT. ASSOCIATED(Mesh % Faces) ) THEN
      CALL Warn('MakeListMatrix','Face dofs requested but not faces exist in mesh!')
      FDOFS = 0
    END IF


    IndexSize = 128
    ALLOCATE( Indexes(IndexSize), STAT=istat )
    IF( istat /= 0 ) THEN
      CALL Fatal('MakeListMatrix','Allocation error for Indexes')
    END IF

    ! Create sparse matrix for the Discontinuous Galerkin solver 
    ! Using either reduced or full basis
    !-------------------------------------------------------------------
    FoundDG = .FALSE.
    IF ( DGSolver ) THEN    
      ! Create the sparse matrix for the Discontinuous Bodies solver
      !-------------------------------------------------------------------
      IF ( DB ) THEN
        DGIndirect = ListGetLogical( Solver % Values,'DG Indirect Connections',Found )

        IF( DGIndirect ) THEN
          CALL Info('MakeListMatrix','Creating also indirect connections!',Level=12)
          ALLOCATE( IndirectPairs( LocalNodes ) )
          IndirectPairs = 0
        END IF
                        

        DO t=1,Mesh % NumberOfBulkElements
          n = 0
          Elm => Mesh % Elements(t)
          IF ( .NOT. CheckElementEquation(Model,Elm,Equation) ) CYCLE

          FoundDG = FoundDG .OR. Elm % DGDOFs > 0
          DO j=1,Elm % DGDOFs
            n = n + 1
            Indexes(n) = Elm % DGIndexes(j)
          END DO

          IF( Elm % DGDofs /= Elm % TYPE % NumberOfNodes ) THEN
            CALL Fatal('MakeListMatrix','Mismatch in sizes in reduced basis DG!')
          END IF

          DO i=1,n
            k1 = Reorder(Indexes(i))
            IF ( k1 <= 0 ) CYCLE
            DO j=1,n
              k2 = Reorder(Indexes(j))
              IF ( k2 <= 0 ) CYCLE
              Lptr => List_GetMatrixIndex( List,k1,k2 )
            END DO
          END DO
        END DO

        IF( Mesh % NumberOfFaces == 0 ) THEN
          DO t=1,Mesh % NumberOfEdges
            n = 0
            Elm => Mesh % Edges(t)
            IF(.NOT. ASSOCIATED( Elm % BoundaryInfo ) ) CYCLE

            Left => Elm % BoundaryInfo % Left
            IF(.NOT. ASSOCIATED( Left ) ) CYCLE
            IF ( .NOT. CheckElementEquation(Model,Left,Equation) ) CYCLE

            Right => Elm % BoundaryInfo %  Right
            IF(.NOT. ASSOCIATED( Right ) ) CYCLE
            IF ( .NOT. CheckElementEquation(Model,Right,Equation) ) CYCLE

            IF( Left % BodyId == Right % BodyId ) CYCLE

            IF( DGIndirect ) THEN
              DO i=1,Left % DGDOFs
                k1 = ReOrder( Left % DgIndexes(i) )
                DO j=1,Right % DGDOFs
                  IF( Left % NodeIndexes(i) == Right % NodeIndexes(j) ) THEN
                    k2 = ReOrder( Right % DgIndexes(j) )
                    IF( k1 /= k2 ) THEN
                      IndirectPairs( k1 ) = k2
                      EXIT
                    END IF
                  END IF
                END DO
              END DO
            END IF

            FoundDG = FoundDG .OR. Left % DGDOFs > 0
            DO j=1,Left % DGDOFs
              n = n + 1
              Indexes(n) = Left % DGIndexes(j)
            END DO

            FoundDG = FoundDG .OR. Right % DGDOFs > 0
            DO j=1,Right % DGDOFs
              n = n + 1
              Indexes(n) = Right % DGIndexes(j)
            END DO

            DO i=1,n
              k1 = Reorder(Indexes(i))
              IF ( k1 <= 0 ) CYCLE
              DO j=1,n
                k2 = Reorder(Indexes(j))
                IF ( k2 <= 0 ) CYCLE
                Lptr => List_GetMatrixIndex( List,k1,k2 )
              END DO
            END DO
          END DO
        END IF

        
        DO t=1,Mesh % NumberOfFaces
          n = 0

          Elm => Mesh % Faces(t)
          IF(.NOT. ASSOCIATED( Elm % BoundaryInfo ) ) CYCLE
          
          Left => Elm  % BoundaryInfo % Left
          IF(.NOT. ASSOCIATED( Left ) ) CYCLE
          IF ( .NOT. CheckElementEquation(Model,Left,Equation) ) CYCLE

          Right => Elm % BoundaryInfo %  Right
          IF(.NOT. ASSOCIATED( Right ) ) CYCLE
          IF ( .NOT. CheckElementEquation(Model,Right,Equation) ) CYCLE

          IF( Left % BodyId == Right % BodyId ) CYCLE

          IF( DGIndirect ) THEN
            DO i=1,Left % DGDOFs
              k1 = ReOrder( Left % DgIndexes(i) )
              DO j=1,Right % DGDOFs
                IF( Left % NodeIndexes(i) == Right % NodeIndexes(j) ) THEN
                  k2 = ReOrder( Right % DgIndexes(j) )
                  IF( k1 /= k2 ) THEN
                    IF( IndirectPairs(k1) > 0 .AND. IndirectPairs(k1) /= k2 ) THEN
                      PRINT *,'Problematic node:',k1,IndirectPairs(k1),k2
                    END IF
                    IndirectPairs( k1 ) = k2
                    EXIT
                  END IF
                END IF
              END DO
            END DO
          END IF            

          FoundDG = FoundDG .OR. Left % DGDOFs > 0
          DO j=1,Left % DGDOFs
            n = n + 1
            Indexes(n) = Left % DGIndexes(j)
          END DO

          FoundDG = FoundDG .OR. Right % DGDOFs > 0
          DO j=1,Right % DGDOFs
            n = n + 1
            Indexes(n) = Right % DGIndexes(j)
          END DO

          DO i=1,n
            k1 = Reorder(Indexes(i))
            IF ( k1 <= 0 ) CYCLE
            DO j=1,n
              k2 = Reorder(Indexes(j))
              IF ( k2 <= 0 ) CYCLE
              Lptr => List_GetMatrixIndex( List,k1,k2 )
            END DO
          END DO
        END DO

        IF( DGIndirect ) THEN
          DO k1 = 1, LocalNodes
            k2 = IndirectPairs(k1)
            IF( k2 == 0 ) CYCLE 
            !PRINT *,'Exchange structure between rows:',k1,k2
            CALL List_ExchangeRowStructure( List, k1, k2 ) 
          END DO
        END IF
          
        
      ELSE
        ! Classical DG solver
        !-------------------------------------        
        DO t=1,Mesh % NumberOfEdges
          n = 0
          Elm => Mesh % Edges(t) % BoundaryInfo % Left
          IF ( ASSOCIATED( Elm ) ) THEN
            IF ( CheckElementEquation(Model,Elm,Equation) ) THEN
              FoundDG = FoundDG .OR. Elm % DGDOFs > 0
              DO j=1,Elm % DGDOFs
                n = n + 1
                Indexes(n) = Elm % DGIndexes(j)
              END DO
            END IF
          END IF

          Elm => Mesh % Edges(t) % BoundaryInfo %  Right
          IF ( ASSOCIATED( Elm ) ) THEN
            IF ( CheckElementEquation(Model,Elm,Equation) ) THEN
              FoundDG = FoundDG .OR. Elm % DGDOFs > 0
              DO j=1,Elm % DGDOFs
                n = n + 1
                Indexes(n) = Elm % DGIndexes(j)
              END DO
            END IF
          END IF

          DO i=1,n
            k1 = Reorder(Indexes(i))
            IF ( k1 <= 0 ) CYCLE
            DO j=1,n
              k2 = Reorder(Indexes(j))
              IF ( k2 <= 0 ) CYCLE
              Lptr => List_GetMatrixIndex( List,k1,k2 )
            END DO
          END DO
        END DO
        DO t=1,Mesh % NumberOfFaces
          n = 0
          Elm => Mesh % Faces(t) % BoundaryInfo % Left
          IF ( ASSOCIATED( Elm ) ) THEN
            IF ( CheckElementEquation(Model,Elm,Equation) ) THEN
              FoundDG = FoundDG .OR. Elm % DGDOFs > 0
              DO j=1,Elm % DGDOFs
                n = n + 1
                Indexes(n) = Elm % DGIndexes(j)
              END DO
            END IF
          END IF

          Elm => Mesh % Faces(t) % BoundaryInfo %  Right
          IF ( ASSOCIATED( Elm ) ) THEN
            IF ( CheckElementEquation(Model,Elm,Equation) ) THEN
              FoundDG = FoundDG .OR. Elm % DGDOFs > 0
              DO j=1,Elm % DGDOFs
                n = n + 1
                Indexes(n) = Elm % DGIndexes(j)
              END DO
            END IF
          END IF

          DO i=1,n
            k1 = Reorder(Indexes(i))
            IF ( k1 <= 0 ) CYCLE
            DO j=1,n
              k2 = Reorder(Indexes(j))
              IF ( k2 <= 0 ) CYCLE
              Lptr => List_GetMatrixIndex( List,k1,k2 )
            END DO
          END DO
        END DO
      END IF
    END IF ! DGSolver

    ! If this is not a GD solver then create permutation considering 
    ! nodal, edge, face and bubble dofs. 
    !-------------------------------------------------------------------
    IF ( .NOT. FoundDG ) THEN
      t = 1
      DO WHILE( t<=Mesh % NumberOfBulkElements+Mesh % NumberOFBoundaryElements )
         Element => Mesh % Elements(t)

         IF ( PRESENT(Equation) ) THEN
           DO WHILE( t<=Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements )
             Element => Mesh % Elements(t)
             IF ( CheckElementEquation(Model,Element,Equation) ) EXIT
             t = t + 1
           END DO
           IF ( t > Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements ) EXIT
         END IF

         n = Element % NDOFs 
         IF( EDOFS > 0 ) n = n + Element % TYPE % NumberOfEdges * EDOFs 
         IF( FDOFS > 0 ) n = n + Element % TYPE % NumberOfFaces * FDOFs
         IF ( GB ) n = n + Element % BDOFs

         IF ( n > IndexSize ) THEN
            IndexSize = n
            IF ( ALLOCATED( Indexes ) ) DEALLOCATE( Indexes )
            ALLOCATE( Indexes(n), STAT=istat )
            IF( istat /= 0 ) THEN
              CALL Fatal('MakeListMatrix','Allocation error for Indexes of size: '//TRIM(I2S(n)))
            END IF
         END IF

         n = 0
         DO i=1,Element % NDOFs
            n = n + 1
            Indexes(n) = Element % NodeIndexes(i)
         END DO

         IF ( EDOFs > 0 ) THEN
            IF ( ASSOCIATED(Element % EdgeIndexes) ) THEN
              DO j=1,Element % TYPE % NumberOFEdges
                 DO i=1, Mesh % Edges(Element % EdgeIndexes(j)) % BDOFs
                   n = n + 1
                   Indexes(n) = EDOFs * (Element % EdgeIndexes(j)-1) + i &
                                + Mesh % NumberOfNodes
                END DO
             END DO
           END IF
         END IF

         IF ( FDOFS > 0 ) THEN
           IF ( ASSOCIATED(Element % FaceIndexes) ) THEN
             DO j=1,Element % TYPE % NumberOFFaces
               DO i=1, Mesh % Faces(Element % FaceIndexes(j)) % BDOFs
                 n = n + 1
                 Indexes(n) = FDOFs*(Element % FaceIndexes(j)-1) + i + &
                     Mesh % NumberOfNodes + EDOFs*Mesh % NumberOfEdges
               END DO
             END DO
           END IF
         END IF

         IF ( GB .AND. ASSOCIATED(Element % BubbleIndexes) ) THEN
            DO i=1,Element % BDOFs
              n = n + 1
              Indexes(n) = FDOFs*Mesh % NumberOfFaces + &
                   Mesh % NumberOfNodes + EDOFs*Mesh % NumberOfEdges + &
                        Element % BubbleIndexes(i)
            END DO
         END IF

         DO i=1,n
            k1 = Reorder(Indexes(i))
            IF ( k1 <= 0 ) CYCLE
            DO j=1,n
               k2 =  Reorder(Indexes(j))
               IF ( k2 <= 0 ) CYCLE
               Lptr => List_GetMatrixIndex( List,k1,k2 )
            END DO
         END DO
         t = t + 1
      END DO

      IF ( ALLOCATED( Indexes ) ) DEALLOCATE( Indexes )
!
!     Diffuse gray radiation condition:
!     ---------------------------------
      Radiation = ListGetLogical( Solver % Values, 'Radiation Solver', Found )
      IF ( .NOT. Found .AND. PRESENT(Equation) ) &
        Radiation = Radiation .OR. (Equation == 'heat equation')

      IF ( Radiation ) THEN
        DO i = Mesh % NumberOfBulkElements+1, &
          Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

          Element => Mesh % Elements(i)
          IF ( ASSOCIATED(Element % BoundaryInfo % GebhardtFactors) ) THEN
             DO j=1,Element % TYPE % NumberOfNodes
                k1 = Reorder(Element % NodeIndexes(j))

                NumberOfFactors = Element % BoundaryInfo % &
                  GebhardtFactors % NumberOfImplicitFactors

                DO n=1,NumberOfFactors

                  Elm => Mesh % Elements( Element % BoundaryInfo % &
                              GebhardtFactors % Elements(n) )

                  DO k=1,Elm % TYPE % NumberOfNodes
                     k2 = Reorder( Elm % NodeIndexes(k) )
                     Lptr => List_GetMatrixIndex( List,k1,k2 )
                  END DO
                END DO
             END DO
          END IF
        END DO
      END IF

      DO i=Mesh % NumberOfBulkElements+1, Mesh % NumberOfBulkElements+ &
                     Mesh % NumberOfBoundaryElements
        IF ( Mesh % Elements(i) % TYPE % ElementCode <  102 .OR. &
             Mesh % Elements(i) % TYPE % ElementCode >= 200 ) CYCLE

        k1 = Reorder( Mesh % Elements(i) % NodeIndexes(1) )
        IF ( k1 > 0 ) THEN
          DO k=1,Mesh % Elements(i) % TYPE % NumberOFNodes
            k2 = Reorder( Mesh % Elements(i) % NodeIndexes(k) )
            IF ( k2 > 0 ) THEN
              Lptr => List_GetMatrixIndex( List,k1,k2 )
              Lptr => List_GetMatrixIndex( List,k2,k1 )
            END IF
          END DO
        END IF

        ! This is a connection element, make a matrix connection for that 
        IF ( Mesh % Elements(i) % TYPE % ElementCode == 102 ) THEN
          k2 = Reorder( Mesh % Elements(i) % NodeIndexes(2) )
          IF ( k2 > 0 ) Lptr => List_GetMatrixIndex( List,k2,k2 )
        END IF
      END DO


      ! Add connection from projectors. These are only needed if the projector is treated
      ! implicely. For explicit projectors or when using Lagrange coefficients the 
      ! connections are not needed. 
      !----------------------------------------------------------------------------------
      IF( DoProjectors ) THEN
        DO This=1,Model % NumberOfBCs
          Projector => Model % BCs(This) % PMatrix
          IF ( .NOT. ASSOCIATED(Projector) ) CYCLE

          IF( ListGetLogical( Model % BCs(This) % Values,&
              'Periodic BC Explicit',Found)) CYCLE
          IF( ListGetLogical( Model % BCs(This) % Values,&
              'Periodic BC Use Lagrange Coefficient',Found)) CYCLE

          CALL Info('MakeListMatrix','Adding matrix topology for BC: '//TRIM(I2S(This)),Level=10)

          DO i=1,Projector % NumberOfRows
            k = Reorder( Projector % InvPerm(i) )
            IF ( k > 0 ) THEN
              DO l=Projector % Rows(i),Projector % Rows(i+1)-1
                IF ( Projector % Cols(l) <= 0 ) CYCLE
                m = Reorder( Projector % Cols(l) )
                IF ( m > 0 ) THEN
                  Lptr => List_GetMatrixIndex( List,k,m )
                  Lptr => List_GetMatrixIndex( List,m,k ) ! keep structure symm.
                  CList => List(k) % Head
                  DO WHILE( ASSOCIATED( CList ) )
                    Lptr => List_GetMatrixIndex( List,m,CList % Index )
                    Lptr => List_GetMatrixIndex( List,CList % Index,m ) ! keep structure symm.
                    CList => CList % Next
                  END DO
                END IF
              END DO
            END IF
          END DO
        END DO
      END IF ! DoProjectors

    END IF
    
    Model % TotalMatrixElements = 0
    DO i=1,LocalNodes
      j = List(i) % Degree
      Model % TotalMatrixElements = Model % TotalMatrixElements + j
    END DO
    
      
    DoNonZeros = .TRUE.
    IF( PRESENT( CalcNonZeros ) ) DoNonZeros = CalcNonZeros
    
    IF( DoNonZeros ) THEN
      ALLOCATE( InvPerm(LocalNodes) )
      InvPerm = 0
      k = 0
      DO i=1,SIZE(Reorder)
        IF (Reorder(i)>0) THEN
          k = k + 1
          InvPerm(Reorder(i)) = k
        END IF
      END DO
      
      Model % Rownonzeros = 0
      DO i=1,LocalNodes
        j = List(i) % Degree
        Model % RowNonzeros(InvPerm(i)) = j
      END DO      
      DEALLOCATE( InvPerm ) 
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE MakeListMatrix
!------------------------------------------------------------------------------
  
!------------------------------------------------------------------------------
!> Create a list matrix array given the mesh, the active domains and the elementtype 
!> related to the solver. The list matrix is flexible since it can account 
!> for any entries. Also constraints and periodic BCs may give rise to entries
!> in the list matrix topology. Multithreaded version using ListMatrixArray for
!> matrix storage.
!------------------------------------------------------------------------------
  SUBROUTINE MakeListMatrixArray( Model,Solver,Mesh,List,Reorder, &
        LocalNodes,Equation, DGSolver, GlobalBubbles, &
        NodalDofsOnly, ProjectorDofs, CalcNonZeros )
!------------------------------------------------------------------------------
    TYPE(Model_t)  :: Model
    TYPE(Mesh_t)   :: Mesh
    TYPE(Matrix_t) :: Matrix
    TYPE(SOlver_t) :: Solver
    TYPE(ListMatrixArray_t) :: List
    INTEGER, OPTIONAL :: Reorder(:)
    LOGICAL, OPTIONAL :: DGSolver
    LOGICAL, OPTIONAL :: GlobalBubbles
    LOGICAL, OPTIONAL :: NodalDofsOnly
    LOGICAL, OPTIONAL :: ProjectorDofs
    CHARACTER(LEN=*), OPTIONAL :: Equation
    LOGICAL, OPTIONAL :: CalcNonZeros
!------------------------------------------------------------------------------
    INTEGER :: t,i,j,k,l,m,k1,k2,n,p,q,e1,e2,f1,f2, EDOFs, FDOFs, BDOFs, This, istat, nthr

    LOGICAL :: Flag, FoundDG, GB, Found, Radiation, DoProjectors, DoNonZeros
    INTEGER, ALLOCATABLE :: InvPerm(:)
    
    TYPE(Matrix_t), POINTER :: Projector

    INTEGER :: IndexSize, NumberOfFactors
    INTEGER, ALLOCATABLE :: Indexes(:), IndexReord(:), IPerm(:)

    TYPE(ListMatrixEntry_t), POINTER :: CList, Lptr

    TYPE(Matrix_t),POINTER :: PMatrix
    TYPE(Element_t), POINTER :: Element,Elm, Edge1, Edge2, Face1, Face2, ElementsList(:)
    TYPE(Graph_t), POINTER :: CurrentColourList
    INTEGER :: CurrentColour, BoundaryColour, CurrentColourStart, &
          CurrentColourEnd, NumberOfMeshColours
    LOGICAL :: NeedLocking
!------------------------------------------------------------------------------

    CALL Info('MakeListMatrixArray','Creating list matrix',Level=14)

    GB = .FALSE.
    IF ( PRESENT(GlobalBubbles) ) GB = GlobalBubbles

    ! No colouring equals a single colour
    NumberOfMeshColours = 1
    IF (ASSOCIATED(Solver % ColourIndexList)) THEN
       NumberOfMeshColours = Solver % ColourIndexList % n
       ! If boundary mesh exists, it has been coloured as well
       IF (ASSOCIATED(Solver % BoundaryColourIndexList)) &
            NumberOfMeshColours = NumberOfMeshColours + Solver % BoundaryColourIndexList % n
    END IF

    nthr=1
    !$ nthr = omp_get_max_threads()
    NeedLocking = (NumberOfMeshColours == 1) .AND. (nthr > 1)
    CALL ListMatrixArray_Allocate(List, LocalNodes, Atomic=NeedLocking)
    
    BDOFs = Mesh % MaxBDOFs
    EDOFs = Mesh % MaxEdgeDOFs
    FDOFs = Mesh % MaxFaceDOFs
    IF( PRESENT( NodalDofsOnly ) ) THEN
      IF( NodalDofsOnly ) THEN
        EDOFS = 0
        FDOFS = 0
      END IF
    END IF

    IF( PRESENT( ProjectorDofs ) ) THEN
      DoProjectors = ProjectorDofs 
    ELSE
      DoProjectors = .TRUE.
    END IF

    IF( EDOFS > 0 .AND. .NOT. ASSOCIATED(Mesh % Edges) ) THEN
      CALL Warn('MakeListMatrixArray','Edge dofs requested but not edges exist in mesh!')
      EDOFS = 0
    END IF

    IF( FDOFS > 0 .AND. .NOT. ASSOCIATED(Mesh % Faces) ) THEN
      CALL Warn('MakeListMatrixArray','Face dofs requested but not faces exist in mesh!')
      FDOFS = 0
    END IF

    ! Create the permutation for the Discontinuous Galerkin solver 
    !-------------------------------------------------------------------
    FoundDG = .FALSE.
    IF ( DGSolver ) THEN
       IndexSize = 128
       ALLOCATE( Indexes(IndexSize), STAT=istat )
       IF( istat /= 0 ) THEN
         CALL Fatal('MakeListMatrixArray','Allocation error for Indexes')
       END IF

       ! TODO: Add multithreading
       DO t=1,Mesh % NumberOfEdges
         n = 0
         Elm => Mesh % Edges(t) % BoundaryInfo % Left
         IF ( ASSOCIATED( Elm ) ) THEN
             IF ( CheckElementEquation(Model,Elm,Equation) ) THEN
                FoundDG = FoundDG .OR. Elm % DGDOFs > 0
                DO j=1,Elm % DGDOFs
                   n = n + 1
                   Indexes(n) = Elm % DGIndexes(j)
                END DO
             END IF
         END IF

         Elm => Mesh % Edges(t) % BoundaryInfo %  Right
         IF ( ASSOCIATED( Elm ) ) THEN
             IF ( CheckElementEquation(Model,Elm,Equation) ) THEN
                FoundDG = FoundDG .OR. Elm % DGDOFs > 0
                DO j=1,Elm % DGDOFs
                   n = n + 1
                   Indexes(n) = Elm % DGIndexes(j)
                END DO
             END IF
         END IF

         DO i=1,n
            k1 = Reorder(Indexes(i))
            IF ( k1 <= 0 ) CYCLE
            DO j=1,n
              k2 = Reorder(Indexes(j))
              IF ( k2 <= 0 ) CYCLE
              
              CALL ListMatrixArray_AddEntry(List, k1, k2)
            END DO
         END DO
      END DO
      DO t=1,Mesh % NumberOfFaces
         n = 0
         Elm => Mesh % Faces(t) % BoundaryInfo % Left
         IF ( ASSOCIATED( Elm ) ) THEN
             IF ( CheckElementEquation(Model,Elm,Equation) ) THEN
                FoundDG = FoundDG .OR. Elm % DGDOFs > 0
                DO j=1,Elm % DGDOFs
                   n = n + 1
                   Indexes(n) = Elm % DGIndexes(j)
                END DO
             END IF
         END IF

         Elm => Mesh % Faces(t) % BoundaryInfo %  Right
         IF ( ASSOCIATED( Elm ) ) THEN
             IF ( CheckElementEquation(Model,Elm,Equation) ) THEN
                FoundDG = FoundDG .OR. Elm % DGDOFs > 0
                DO j=1,Elm % DGDOFs
                   n = n + 1
                   Indexes(n) = Elm % DGIndexes(j)
                END DO
             END IF
         END IF

         DO i=1,n
            k1 = Reorder(Indexes(i))
            IF ( k1 <= 0 ) CYCLE
            DO j=1,n
              k2 = Reorder(Indexes(j))
              IF ( k2 <= 0 ) CYCLE
              
              CALL ListMatrixArray_AddEntry(List, k1, k2)
            END DO
         END DO
      END DO
      DEALLOCATE(Indexes)
    END IF

    ! If this is not a GD solver then create permutation considering 
    ! nodal, edge, face and bubble dofs. 
    !-------------------------------------------------------------------
    IF (.NOT. FoundDG) THEN
      
      !$OMP PARALLEL &
      !$OMP SHARED(LocalNodes, List, ListNew, Equation, EDOFS, FDOFS, GB, &
      !$OMP        Reorder, Model, Mesh, NumberOfMeshColours, CurrentColourStart, NeedLocking, &
      !$OMP        CurrentColourEnd, CurrentColourList, ElementsList, Solver, BoundaryColour) &
      !$OMP PRIVATE(Element, Indexes, istat, IndexSize, IndexReord, &
      !$OMP         IPerm, n, i, j, nReord, k1, k2, Lptr, LptrNew, Errors, &
      !$OMP         CurrentColour) &
      !$OMP DEFAULT(NONE)

      IndexSize = 0

      ! Loop over mesh colours
      DO CurrentColour=1,NumberOfMeshColours

         !$OMP SINGLE
         ! Test if matrix is actually coloured or not
         IF (NumberOfMeshColours == 1) THEN
           CurrentColourList => NULL()
           ElementsList => Mesh % Elements(1:Mesh % NumberOfBulkElements+Mesh % NumberOFBoundaryElements)
           CurrentColourStart = 1
           CurrentColourEnd = Mesh % NumberOfBulkElements+Mesh % NumberOFBoundaryElements
         ELSE IF (CurrentColour <= Solver % ColourIndexList % n) THEN
           CALL Info('MakeListMatrixArray','ListMatrix add colour: '//TRIM(I2S(CurrentColour)),Level=10)
           CurrentColourList => Solver % ColourIndexList
           ElementsList => Mesh % Elements(1:Mesh % NumberOfBulkElements)
           CurrentColourStart = CurrentColourList % Ptr(CurrentColour)
           CurrentColourEnd = CurrentColourList % Ptr(CurrentColour+1)-1
         ELSE
           BoundaryColour = CurrentColour-Solver % ColourIndexList % n
           CALL Info('MakeListMatrixArray','ListMatrix add boundary colour: '//TRIM(I2S(BoundaryColour)),Level=10)

           CurrentColourList => Solver % BoundaryColourIndexList
           ! Boundary elements are stored after bulk elements in Mesh
           ElementsList => Mesh % Elements(Mesh % NumberOfBulkElements+1:&
                 Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements)
           CurrentColourStart = CurrentColourList % Ptr(BoundaryColour)
           CurrentColourEnd = CurrentColourList % Ptr(BoundaryColour+1)-1
         END IF
         !$OMP END SINGLE
         
         !$OMP DO
         DO t=CurrentColourStart,CurrentColourEnd
            IF (ASSOCIATED(CurrentColourList)) THEN
               Element => ElementsList(CurrentColourList % ind(t))
            ELSE
               Element => ElementsList(t)
            END IF


            IF ( PRESENT(Equation) ) THEN
               IF ( .NOT. CheckElementEquation(Model,Element,Equation) ) CYCLE
            END IF

            n = Element % NDOFs 
            IF( EDOFS > 0 ) n = n + Element % TYPE % NumberOfEdges * EDOFs 
            IF( FDOFS > 0 ) n = n + Element % TYPE % NumberOfFaces * FDOFs
            IF ( GB ) n = n + Element % BDOFs
            
            ! Reallocate index permutation if needed
            IF ( n > IndexSize ) THEN
               IF (ALLOCATED(Indexes)) DEALLOCATE(Indexes, IndexReord, IPerm)
               
               IndexSize = MAX(MAX(128, IndexSize*2), n)
               ALLOCATE(Indexes(IndexSize), &
                    IndexReord(IndexSize), &
                    IPerm(IndexSize), STAT=istat )
               IF( istat /= 0 ) THEN
                  CALL Fatal('MakeListMatrixArray','Allocation error for Indexes of size: '//TRIM(I2S(n)))
               END IF
            END IF
            
            n = 0
            DO i=1,Element % NDOFs
               n = n + 1
               Indexes(n) = Element % NodeIndexes(i)
            END DO
            
            IF ( EDOFs > 0 ) THEN
               IF ( ASSOCIATED(Element % EdgeIndexes) ) THEN
                  DO j=1,Element % TYPE % NumberOFEdges
                     DO i=1, Mesh % Edges(Element % EdgeIndexes(j)) % BDOFs
                        n = n + 1
                        Indexes(n) = EDOFs * (Element % EdgeIndexes(j)-1) + i &
                             + Mesh % NumberOfNodes
                     END DO
                  END DO
               END IF
            END IF
            
            IF ( FDOFS > 0 ) THEN
               IF ( ASSOCIATED(Element % FaceIndexes) ) THEN
                  DO j=1,Element % TYPE % NumberOFFaces
                     DO i=1, Mesh % Faces(Element % FaceIndexes(j)) % BDOFs
                        n = n + 1
                        Indexes(n) = FDOFs*(Element % FaceIndexes(j)-1) + i + &
                             Mesh % NumberOfNodes + EDOFs*Mesh % NumberOfEdges
                     END DO
                  END DO
               END IF
            END IF
            
            IF ( GB .AND. ASSOCIATED(Element % BubbleIndexes) ) THEN
               DO i=1,Element % BDOFs
                  n = n + 1
                  Indexes(n) = FDOFs*Mesh % NumberOfFaces + &
                       Mesh % NumberOfNodes + EDOFs*Mesh % NumberOfEdges + &
                       Element % BubbleIndexes(i)
               END DO
            END IF

            nReord=0
            DO i=1,n
               IF (Reorder(Indexes(i)) > 0) THEN
                  nReord=nReord+1
                  IndexReord(nReord)=Reorder(Indexes(i))
               END IF
            END DO
            CALL InsertionSort(nReord, IndexReord, IPerm)

            DO i=1,nReord
               k1 = IndexReord(IPerm(i))
             
               ! Bulk add all sorted nonzero indices to a matrix row in one go
               CALL ListMatrixArray_AddEntries(List, k1, nReord, IndexReord, IPerm, NeedLocking)
            END DO
         END DO
         !$OMP END DO
      END DO ! Loop over mesh colours

      IF (ALLOCATED(Indexes)) DEALLOCATE(Indexes, IndexReord, IPerm)
      !$OMP END PARALLEL

!
!     Diffuse gray radiation condition:
!     ---------------------------------
      Radiation = ListGetLogical( Solver % Values, 'Radiation Solver', Found )
      IF ( .NOT. Found .AND. PRESENT(Equation) ) &
        Radiation = Radiation .OR. (Equation == 'heat equation')

      IF ( Radiation ) THEN
        DO i = Mesh % NumberOfBulkElements+1, &
          Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

          Element => Mesh % Elements(i)
          IF ( ASSOCIATED(Element % BoundaryInfo % GebhardtFactors) ) THEN
             DO j=1,Element % TYPE % NumberOfNodes
                k1 = Reorder(Element % NodeIndexes(j))

                NumberOfFactors = Element % BoundaryInfo % &
                  GebhardtFactors % NumberOfImplicitFactors

                DO n=1,NumberOfFactors

                  Elm => Mesh % Elements( Element % BoundaryInfo % &
                              GebhardtFactors % Elements(n) )

                  DO k=1,Elm % TYPE % NumberOfNodes
                     k2 = Reorder( Elm % NodeIndexes(k) )
                     CALL ListMatrixArray_AddEntry(List, k1, k2)
                  END DO
                END DO
             END DO
          END IF
        END DO
      END IF

      ! TODO: Add multithreading
      DO i=Mesh % NumberOfBulkElements+1, Mesh % NumberOfBulkElements+ &
                     Mesh % NumberOfBoundaryElements
        IF ( Mesh % Elements(i) % TYPE % ElementCode <  102 .OR. &
             Mesh % Elements(i) % TYPE % ElementCode >= 200 ) CYCLE

        k1 = Reorder( Mesh % Elements(i) % NodeIndexes(1) )
        IF ( k1 > 0 ) THEN
          DO k=1,Mesh % Elements(i) % TYPE % NumberOFNodes
            k2 = Reorder( Mesh % Elements(i) % NodeIndexes(k) )
            IF ( k2 > 0 ) THEN
              CALL ListMatrixArray_AddEntry(List, k1, k2)
              CALL ListMatrixArray_AddEntry(List, k2, k1)
            END IF
          END DO
        END IF

        ! This is a connection element, make a matrix connection for that 
        IF ( Mesh % Elements(i) % TYPE % ElementCode == 102 ) THEN
          k2 = Reorder( Mesh % Elements(i) % NodeIndexes(2) )
          IF ( k2 > 0 ) CALL ListMatrixArray_AddEntry(List, k2, k2)
        END IF
      END DO


      ! Add connection from projectors. These are only needed if the projector is treated
      ! implicely. For explicit projectors or when using Lagrange coefficients the 
      ! connections are not needed. 
      !----------------------------------------------------------------------------------
      IF( DoProjectors ) THEN
        DO This=1,Model % NumberOfBCs
          Projector => Model % BCs(This) % PMatrix
          IF ( .NOT. ASSOCIATED(Projector) ) CYCLE

          IF( ListGetLogical( Model % BCs(This) % Values,&
              'Periodic BC Explicit',Found)) CYCLE
          IF( ListGetLogical( Model % BCs(This) % Values,&
              'Periodic BC Use Lagrange Coefficient',Found)) CYCLE

          CALL Info('MakeListMatrixArray','Adding matrix topology for BC: '//TRIM(I2S(This)),Level=10)

          ! TODO: Add multithreading
          DO i=1,Projector % NumberOfRows
            k = Reorder( Projector % InvPerm(i) )
            IF ( k > 0 ) THEN
              DO l=Projector % Rows(i),Projector % Rows(i+1)-1
                IF ( Projector % Cols(l) <= 0 ) CYCLE
                m = Reorder( Projector % Cols(l) )
                IF ( m > 0 ) THEN
                  CALL ListMatrixArray_AddEntry(List, k, m)
                  CALL ListMatrixArray_AddEntry(List, m, k)
                  CList => List % Rows(k) % Head
                  DO WHILE( ASSOCIATED( CList ) )
                    CALL ListMatrixArray_AddEntry(List, m, CList % Index)
                    CALL ListMatrixArray_AddEntry(List, CList % Index, m)
                    CList => CList % Next
                  END DO
                END IF
              END DO
            END IF
          END DO
        END DO
      END IF ! DoProjectors

    END IF ! (.NOT. FoundDG)



    Model % TotalMatrixElements = 0
    DO i=1,LocalNodes
      j = List % Rows(i) % Degree
      Model % TotalMatrixElements = Model % TotalMatrixElements + j
    END DO
    
      
    DoNonZeros = .TRUE.
    IF( PRESENT( CalcNonZeros ) ) DoNonZeros = CalcNonZeros
    
    IF( DoNonZeros ) THEN
      ALLOCATE( InvPerm(LocalNodes) )
      InvPerm = 0
      k = 0
      DO i=1,SIZE(Reorder)
        IF (Reorder(i)>0) THEN
          k = k + 1
          InvPerm(Reorder(i)) = k
        END IF
      END DO
      
      Model % Rownonzeros = 0
      DO i=1,LocalNodes
        j = List % Rows(i) % Degree
        Model % RowNonzeros(InvPerm(i)) = j
      END DO      
      DEALLOCATE( InvPerm ) 
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE MakeListMatrixArray
!------------------------------------------------------------------------------

  

!------------------------------------------------------------------------------
!>    Initialize a CRS format matrix to the effect that it will be ready to
!>    accept values when CRS_GlueLocalMatrix is called (build up the index
!>    tables of a CRS format matrix)....
!------------------------------------------------------------------------------
  SUBROUTINE InitializeMatrix( Matrix, n, List, DOFs, Reorder, InvInitialReorder )
!------------------------------------------------------------------------------
    INTEGER :: DOFs, n
    TYPE(Matrix_t),POINTER :: Matrix
    TYPE(ListMatrix_t) :: List(:)
    INTEGER, OPTIONAL :: Reorder(:), InvInitialReorder(:)
!------------------------------------------------------------------------------
    TYPE(ListMatrixEntry_t), POINTER :: Clist
    INTEGER :: i,j,k,l,m,k1,k2
    INTEGER, POINTER CONTIG :: Rows(:), Cols(:)
    LOGICAL :: DoReorder 

    
    Rows => Matrix % Rows
    Cols => Matrix % Cols
    
    DoReorder = .FALSE.
    IF( PRESENT( Reorder ) ) THEN
      IF( .NOT. PRESENT( InvInitialReorder ) ) THEN
        CALL Fatal('InitializeMatrix','Need both old and new numbering!')
      END IF
      DoReorder = .TRUE.
    END IF
    
    
    IF( DoReorder ) THEN
      ! In case of reordering we need to compute the row offset in advance
      Rows(1) = 1
      DO i=1,n       
        DO l=1,DOFs
          j = Reorder( InvInitialReorder(i) )
          k1 = DOFs * (j-1) + l
          
          Rows(k1+1) = Dofs * List(i) % Degree
        END DO
      END DO
      DO i=1,Dofs*n
        Rows(i+1) = Rows(i) + Rows(i+1)
      END DO
      
      !$OMP PARALLEL DO SHARED(Rows, Cols, List, n, DOFs, Reorder, InvInitialReorder) &
      !$OMP PRIVATE(CList, l, j, k1, k2, k, m) &
      !$OMP DEFAULT(NONE)
      DO i=1,n
        DO l=1,DOFs
          CList => List(i) % Head
          j = Reorder( InvInitialReorder(i) )
          k1 = DOFs * (j-1) + l
          k2 = Rows(k1)-1
          
          DO WHILE( ASSOCIATED( CList ) )
            k = Reorder( InvInitialReorder(Clist % INDEX))
            k = DOFs*(k-1)
            DO m=k+1,k+DOFs
              k2 = k2+1
              Cols(k2) = m
            END DO
            CList => Clist % Next
          END DO
        END DO
      END DO
      !$OMP END PARALLEL DO      
    ELSE
      ! If there is no renumbering then the reordering is one-to-one mapping
      !$OMP PARALLEL DO SHARED(Rows, Cols, List, n, DOFs) &
      !$OMP PRIVATE(CList, l, j, k1, k2, k, m) &
      !$OMP DEFAULT(NONE)

      DO i=1,n
        DO l=1,DOFs
          CList => List(i) % Head
          j = i
          k1 = DOFs * (j-1) + l
          k2 = Rows(k1)-1
          
          DO WHILE( ASSOCIATED( CList ) )
            k = Clist % index 
            k = DOFs*(k-1)
            DO m=k+1,k+DOFs
              k2 = k2+1
              Cols(k2) = m
            END DO
            CList => Clist % Next
          END DO
          
          Rows(k1+1) = k2+1
        END DO
      END DO
      !$OMP END PARALLEL DO      
    END IF
    
    IF ( Matrix % FORMAT == MATRIX_CRS ) CALL CRS_SortMatrix( Matrix )
!------------------------------------------------------------------------------
  END SUBROUTINE InitializeMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   FUNCTION CreateMatrix( Model, Solver, Mesh, Perm, DOFs, MatrixFormat, &
          OptimizeBW, Equation, DGSolver, GlobalBubbles, &
          NodalDofsOnly, ProjectorDofs, ThreadedStartup ) RESULT(Matrix)
!------------------------------------------------------------------------------
     IMPLICIT NONE
     TYPE(Model_t) :: Model
     TYPE(Mesh_t) :: Mesh
     TYPE(Solver_t), TARGET :: Solver
     INTEGER :: DOFs, MatrixFormat
     INTEGER, TARGET :: Perm(:)
     LOGICAL :: OptimizeBW
     LOGICAL, OPTIONAL :: DGSolver, GlobalBubbles
     LOGICAL, OPTIONAL :: NodalDofsOnly, ProjectorDofs
     LOGICAL, OPTIONAL :: ThreadedStartup

     CHARACTER(LEN=*), OPTIONAL :: Equation

     TYPE(Matrix_t),POINTER :: Matrix
!------------------------------------------------------------------------------
     TYPE(ListMatrix_t), POINTER :: ListMatrix(:)
     TYPE(ListMatrixArray_t) :: ListMatrixArray
     TYPE(Matrix_t), POINTER :: A
     TYPE(Element_t), POINTER :: Element
     TYPE(ListMatrixEntry_t), POINTER :: CList
     CHARACTER(LEN=MAX_NAME_LEN) :: Eq, str
     LOGICAL :: GotIt, DG, GB, UseOptimized, Found
     INTEGER i,j,k,l,k1,t,n, p,m, minEdgeDOFs, maxEdgeDOFs, &
           minFaceDOFs, maxFaceDOFs, BDOFs, cols, istat
     INTEGER, POINTER :: Ivals(:)
     INTEGER, ALLOCATABLE :: InvInitialReorder(:)
     INTEGER :: nthr
     LOGICAL :: UseThreads
     LOGICAL, ALLOCATABLE :: ConstrainedNode(:)
     
!------------------------------------------------------------------------------

     NULLIFY( Matrix )

     DG = .FALSE.
     IF ( PRESENT(DGSolver) )  DG = DGSolver

     GB = .FALSE.
     IF ( PRESENT(GlobalBubbles) ) GB = GlobalBubbles
       
     IF( OptimizeBW ) THEN
       IF( ListGetLogical( Solver % Values,'DG Reduced Basis',Found ) ) THEN
         CALL Info('CreateMatrix','Suppressing bandwidth optimization for discontinuous bodies',Level=8)
         OptimizeBW = .FALSE.
       END IF
     END IF

     
     UseThreads = .FALSE.
     nthr = 1
     IF ( PRESENT(ThreadedStartup) ) THEN
#ifdef _OPENMP
       IF (ThreadedStartup) THEN
          UseThreads = .TRUE.
          nthr = omp_get_max_threads()
       END IF
#endif
     END IF

     minEdgeDOFs = HUGE(minEdgeDOFs)
     maxEdgeDOFs = 0
     minFaceDOFs = HUGE(minFaceDOFs)
     maxFaceDOFs = 0
     BDOFs = 0

     !$OMP PARALLEL SHARED(Mesh) &
     !$OMP          REDUCTION(min:minEdgeDOFs) REDUCTION(max:maxEdgeDOFs) &
     !$OMP          REDUCTION(min:minFaceDOFs) REDUCTION(max:maxFaceDOFs) &
     !$OMP          REDUCTION(max:BDOFs) &
     !$OMP          DEFAULT(NONE) NUM_THREADS(nthr)
     
     !$OMP DO 
     DO i=1,Mesh % NumberOfEdges
        minEdgeDOFs = MIN( minEdgeDOFs, Mesh % Edges(i) % BDOFs )
        maxEdgeDOFs = MAX( maxEdgeDOFs, Mesh % Edges(i) % BDOFs )
     END DO
     !$OMP END DO NOWAIT
     !$OMP DO
     DO i=1,Mesh % NumberOfFaces
        minFaceDOFs = MIN( minFaceDOFs, Mesh % Faces(i) % BDOFs )
        maxFaceDOFs = MAX( maxFaceDOFs, Mesh % Faces(i) % BDOFs )
     END DO
     !$OMP END DO NOWAIT
     !$OMP DO
     DO i=1,Mesh % NumberOfBulkElements
        BDOFs = MAX( BDOFs, Mesh % Elements(i) % BDOFs )
     END DO
     !$OMP END DO NOWAIT
     !$OMP END PARALLEL
     
     Mesh % MaxEdgeDOFs = maxEdgeDOFs
     IF(minEdgeDOFs <= maxEdgeDOFs ) THEN
       Mesh % MinEdgeDOFs = minEdgeDOFs
     ELSE
       Mesh % MinEdgeDOFs = maxEdgeDOFs
     END IF

     Mesh % MaxFaceDOFs = maxFaceDOFs
     IF(minFaceDOFs <= maxFaceDOFs ) THEN
       Mesh % MinFaceDOFs = minFaceDOFs
     ELSE
       Mesh % MinFaceDOFs = maxFaceDOFs
     END IF
     
     Mesh % MaxBDOFs = BDOFs
     
     IF ( PRESENT( Equation ) ) n = StringToLowerCase( Eq,Equation )

     Perm = 0
     IF ( PRESENT(Equation) ) THEN
       CALL Info('CreateMatrix','creating initial permutation',Level=14)
       k = InitialPermutation( Perm,Model,Solver,Mesh,Eq,DG,GB )
       IF ( k <= 0 ) THEN
         RETURN
       END IF
     ELSE
       k = SIZE( Perm )
     END IF
     
     IF ( k == SIZE(Perm) ) THEN
       IF(PRESENT(NodalDofsOnly)) THEN
         IF(NodalDofsOnly) k=Mesh % NumberOfNodes
       END IF       
       DO i=1,k 
         Perm(i) = i
       END DO
     END IF

     IF( ParEnv % PEs > 1 .AND. &
         ListGetLogical( Solver % Values,'Skip Pure Halo Nodes',Found ) ) THEN
       CALL Info('CreateMatrix','Skipping pure halo nodes',Level=14)
       j = 0
       DO i=1,Mesh % NumberOfNodes 
         ! These are pure halo nodes that need not be communicated. They are created only 
         ! for sufficient geometric information on the boundaries.
         IF( .NOT. ANY( ParEnv % Mype == Mesh % ParallelInfo % NeighbourList(i) % Neighbours ) ) THEN
           Perm(i) = 0
         ELSE IF( Perm(i) > 0 ) THEN
           j = j + 1
           Perm(i) = j
         END IF
       END DO
       PRINT *,'Eliminating '//TRIM(I2S(k-j))//' halo nodes out of '&
           //TRIM(I2S(k))//' in partition '//TRIM(I2S(ParEnv % MyPe))
       k = j
     END IF

     
     IF( OptimizeBW ) THEN
       CALL Info('CreateMatrix','Creating inverse of initial order of size: '//TRIM(I2S(k)),Level=14)
       ALLOCATE( InvInitialReorder(k), STAT=istat )
       IF( istat /= 0 ) THEN
         CALL Fatal('CreateMatrix','Allocation error for InvInitialReorder of size: '//TRIM(I2S(k)))
       END IF

       ! We need to keep the initial numbering only in case we optimize the bandwidth!
       InvInitialReorder = 0
       DO i=1,SIZE(Perm)
         IF (Perm(i)>0) InvInitialReorder(Perm(i)) = i
       END DO
     END IF
     
     UseOptimized = ListGetLogical( Solver % Values, &
       'Optimize Bandwidth Use Always', GotIt )
          
     Matrix => NULL()

     ! check if matrix structures really need to be created:
     ! -----------------------------------------------------
     IF ( ListGetLogical( Solver % Values, 'No matrix',GotIt)) RETURN

     !------------------------------------------------------------------------------
     ! Note that Model % RowNonZeros is not used anymore!!!!
     ! For this to be ok the "SetRows" flag must be set .FALSE.
     ! The reason behind the change is that using that was not flexible enough with
     ! DB fields. 
     !------------------------------------------------------------------------------
     !     ALLOCATE( Model % RowNonZeros(k), STAT=istat )
     !     IF( istat /= 0 ) THEN
     !       CALL Fatal('CreateMatrix','Allocation error for RowNonZeros of size: '//TRIM(I2S(k)))
     !     END IF
     !     Model % RowNonzeros=0

     IF (UseThreads) THEN
       CALL Info('CreateMatrix','Creating threaded list matrix array for equation',Level=14)
       IF ( PRESENT(Equation) ) THEN
         CALL MakeListMatrixArray( Model, Solver, Mesh, ListMatrixArray, Perm, k, Eq, DG, GB,&
               NodalDofsOnly, ProjectorDofs, CalcNonZeros = .FALSE. )
         n = OptimizeBandwidth( ListMatrixArray % Rows, Perm, InvInitialReorder, &
               k, OptimizeBW, UseOptimized, Eq )
       ELSE
         CALL MakeListMatrixArray( Model, Solver, Mesh, ListMatrixArray, Perm, k, &
               DGSolver=DG, GlobalBubbles=GB, NodalDofsOnly=NodalDofsOnly, &
               ProjectorDofs=ProjectorDofs, CalcNonZeros = .FALSE. )
         n = OptimizeBandwidth( ListMatrixArray % Rows, Perm, InvInitialReorder, &
               k, OptimizeBW,UseOptimized, ' ' )
       END IF
       
       !------------------------------------------------------------------------------
       ! Initialize the matrix. Multithreading only supports CRS.
       !------------------------------------------------------------------------------
       CALL Info('CreateMatrix','Initializing list matrix array for equation',Level=14)
       IF ( MatrixFormat == MATRIX_CRS) THEN
         Matrix => CRS_CreateMatrix( DOFs*k, Model % TotalMatrixElements, Ndeg=DOFs, &
             Reorder=Perm, AllocValues=.TRUE., SetRows = .FALSE.)
         Matrix % FORMAT = MatrixFormat
         IF( OptimizeBW ) THEN
           CALL InitializeMatrix( Matrix, k, ListMatrixArray % Rows, &
               DOFs, Perm, InvInitialReorder )
         ELSE
           CALL InitializeMatrix( Matrix, k, ListMatrixArray % Rows, DOFs )          
         END IF
       ELSE
         CALL Fatal('CreateMatrix','Multithreaded startup only supports CRS matrix format')
       END IF
       
       CALL Info('CreateMatrix','Matrix created',Level=14)

       CALL ListMatrixArray_Free( ListMatrixArray )       
     ELSE
       NULLIFY( ListMatrix )
       CALL Info('CreateMatrix','Creating list matrix for equation',Level=14)
       IF ( PRESENT(Equation) ) THEN
         CALL MakeListMatrix( Model, Solver, Mesh, ListMatrix, Perm, k, Eq, DG, GB,&
               NodalDofsOnly, ProjectorDofs, CalcNonZeros = .FALSE.)
         n = OptimizeBandwidth( ListMatrix, Perm, InvInitialReorder, &
               k, OptimizeBW, UseOptimized, Eq )
       ELSE
         CALL MakeListMatrix( Model, Solver, Mesh, ListMatrix, Perm, k, &
               DGSolver=DG, GlobalBubbles=GB, NodalDofsOnly=NodalDofsOnly, &
               ProjectorDofs=ProjectorDofs, CalcNonZeros = .FALSE.)
         n = OptimizeBandwidth( ListMatrix, Perm, InvInitialReorder, &
               k, OptimizeBW,UseOptimized, ' ' )
       END IF
       
       !------------------------------------------------------------------------------
       ! Initialize the matrix. 
       !------------------------------------------------------------------------------
       CALL Info('CreateMatrix','Initializing list matrix for equation',Level=14)
       SELECT CASE( MatrixFormat )
       CASE( MATRIX_CRS )
         Matrix => CRS_CreateMatrix( DOFs*k, Model % TotalMatrixElements, Ndeg=DOFs, &
             Reorder=Perm, AllocValues=.TRUE., SetRows = .FALSE.)
         Matrix % FORMAT = MatrixFormat
         IF( OptimizeBW ) THEN
           CALL InitializeMatrix( Matrix, k, ListMatrix, &
               DOFs, Perm, InvInitialReorder )
         ELSE
           CALL InitializeMatrix( Matrix, k, ListMatrix, DOFs )
         END IF
            
      CASE( MATRIX_BAND )
        Matrix => Band_CreateMatrix( DOFs*k, DOFs*n,.FALSE.,.TRUE. )
         
       CASE( MATRIX_SBAND )
         Matrix => Band_CreateMatrix( DOFs*k, DOFs*n,.TRUE.,.TRUE. )
       END SELECT
       CALL Info('CreateMatrix','Matrix created',Level=14)

       CALL List_FreeMatrix( k, ListMatrix )
     END IF
     
     NULLIFY( Matrix % MassValues, Matrix % DampValues, Matrix % Force, Matrix % RHS_im )
!------------------------------------------------------------------------------
     Matrix % Solver => Solver
     Matrix % DGMatrix = DG
     Matrix % Subband = DOFs * n
     Matrix % COMPLEX = .FALSE.
     Matrix % FORMAT  = MatrixFormat
!------------------------------------------------------------------------------

     n = ListGetInteger( Solver % Values, 'Constraint DOFs', Found )
     IF ( n>0 ) THEN       
       Matrix % ConstraintMatrix => AllocateMatrix()
       A => Matrix % ConstraintMatrix
       A % NumberOfRows = n
       ALLOCATE( A % Rows(n+1), A % Diag(n), A % RHS(n), &
           ConstrainedNode(Mesh % NumberOfNodes), STAT=istat )
       IF( istat /= 0 ) THEN
         CALL Fatal('CreateMatrix','Allocation error for CRS matrix topology: '//TRIM(I2S(n)))
       END IF

       DO i=1,n
         A % RHS(i:i) = ListGetConstReal( Solver % Values,  &
           'Constraint DOF ' // TRIM(i2s(i)) // ' Value' )
       END DO

       Cols = 0
       A % Rows(1) = 1
       DO i=1,n
         WRITE( str, '(a)' ) 'Constraint DOF ' // TRIM(i2s(i)) // ' Body' 
         ivals => ListGetIntegerArray( Solver % Values, str, Found )
         IF ( ASSOCIATED(ivals) ) THEN
           ConstrainedNode = .FALSE.
           DO k=1,Solver % Mesh % NumberOfBulkElements
             Element => Solver % Mesh % Elements(k)
             IF ( ALL(ivals /= Element % Bodyid) ) CYCLE
             IF ( ALL( Perm(Element % NodeIndexes) > 0 ) ) &
               ConstrainedNode(Element % NodeIndexes) = .TRUE.
           END DO
           Cols = Cols+DOFs*COUNT(ConstrainedNode)
         END IF

         WRITE( str, '(a)' ) 'Constraint DOF ' // TRIM(i2s(i)) // ' BC'
         Ivals => ListGetIntegerArray( Solver % Values, str, Found )
         IF ( ASSOCIATED(Ivals) ) THEN
           ConstrainedNode = .FALSE.
           DO k=Solver % Mesh % NumberOfBulkElements+1, &
                Solver % Mesh % NumberOfBulkElements+Solver % Mesh % NumberOfBoundaryElements
             Element => Solver % Mesh % Elements(k)
             IF ( ALL(Element % Boundaryinfo % Constraint /= ivals) ) CYCLE
             IF ( ALL( Perm(Element % NodeIndexes) > 0 ) ) &
                 ConstrainedNode(Element % NodeIndexes) = .TRUE.
           END DO
           Cols = Cols+DOFs*COUNT(ConstrainedNode)
         END IF
         A % Rows(i+1) = A % Rows(i)+Cols
       END DO

       ALLOCATE( A % Cols(cols), A % Values(cols), STAT=istat )
       IF( istat /= 0 ) THEN
         CALL Fatal('CreateMatrix','Allocation error for CRS cols and values: '//TRIM(I2S(cols)))
       END IF
       A % Cols = 0
       A % Values = 0

       DO i=1,n
         WRITE( str, '(a)' ) 'Constraint DOF ' // TRIM(i2s(i)) // ' Body' 
         ivals => ListGetIntegerArray( Solver % Values, str, Found )
         IF ( ASSOCIATED(ivals) ) THEN
           DO k=1,Solver % Mesh % NumberOfBulkElements
             Element => Solver % Mesh % Elements(k)
             IF ( ALL(ivals /= Element % Bodyid) ) CYCLE

             IF ( ALL( Perm(Element % NodeIndexes) > 0 ) ) THEN
               DO p=1,Element % TYPE % NumberOfNodes
                 l = Perm(Element % NodeIndexes(p))
                 DO m=1,DOFs 
                   k1 = DOFs*(l-1)+m
                   CALL CRS_MakeMatrixIndex( A,i,k1 )
                 END DO
               END DO
             END IF
           END DO
         END IF

         WRITE( str, '(a)' ) 'Constraint DOF ' // TRIM(i2s(i)) // ' BC'
         ivals => ListGetIntegerArray( Solver % Values, str, Found )
         IF ( ASSOCIATED(ivals) ) THEN
           DO k=Solver % Mesh % NumberOfBulkElements+1, &
                Solver % Mesh % NumberOfBulkElements+Solver % Mesh % NumberOfBoundaryElements
             Element => Solver % Mesh % Elements(k)
             IF ( ALL(ivals /= Element % Boundaryinfo % Constraint) ) CYCLE
             IF ( ALL(Perm(Element % NodeIndexes) > 0) ) THEN
               DO p=1,Element % TYPE % NumberOfNodes
                 l = Perm(Element % NodeIndexes(p))
                 DO m=1,DOFs 
                   k1 = DOFs*(l-1)+m
                   CALL CRS_MakeMatrixIndex( A,i,k1 )
                 END DO
               END DO
             END IF
           END DO
         END IF
       END DO
       CALL CRS_SortMatrix(A)
     END IF

!     DEALLOCATE( Model % RowNonZeros )
     IF( OptimizeBW ) DEALLOCATE( InvInitialReorder )
!------------------------------------------------------------------------------
   END FUNCTION CreateMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   FUNCTION CreateOdeMatrix( Model, Solver, Dofs ) RESULT(Matrix)
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     TYPE(Solver_t), TARGET :: Solver
     INTEGER :: DOFs
     TYPE(Matrix_t), POINTER :: Matrix
!------------------------------------------------------------------------------
     LOGICAL :: Found
     INTEGER i,j,k
!------------------------------------------------------------------------------

     Matrix => NULL()

     IF ( ListGetLogical( Solver % Values, 'No matrix',Found)) RETURN
     
     ! Create a list matrix that allows for unspecified entries in the matrix 
     ! structure to be introduced.
     Matrix => AllocateMatrix()
     Matrix % FORMAT = MATRIX_LIST
     
     ! Initialize matrix indices
     DO i = 1, Dofs
       DO j = 1, Dofs
          CALL List_AddMatrixIndex(Matrix % ListMatrix, i, j) 
       END DO
     END DO

     CALL List_ToCRSMatrix(Matrix)
     CALL CRS_SortMatrix(Matrix,.TRUE.)
     
     CALL Info('CreateOdeMatrix','Number of rows in ode matrix: '//&
         TRIM(I2S(Matrix % NumberOfRows)), Level=9)
     CALL Info('CreateOdeMatrix','Number of entries in ode matrix: '//&
         TRIM(I2S(SIZE(Matrix % Cols)) ), Level=9)
     
     Matrix % Solver => Solver
     Matrix % DGMatrix = .FALSE.
     Matrix % Subband = DOFs
     Matrix % COMPLEX = .FALSE.
     ! Matrix % FORMAT  = MatrixFormat

!------------------------------------------------------------------------------
   END FUNCTION CreateOdeMatrix
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE RotateMatrix( Matrix,Vector,n,DIM,DOFs,NodeIndexes,  &
                   Normals,Tangent1,Tangent2 )
!------------------------------------------------------------------------------

    REAL(KIND=dp) :: Matrix(:,:),Vector(:)
    REAL(KIND=dp), POINTER :: Normals(:,:), Tangent1(:,:),Tangent2(:,:)
    INTEGER :: n,DIM,DOFs,NodeIndexes(:)
!------------------------------------------------------------------------------

    INTEGER :: i,j,k,l
    REAL(KIND=dp) :: s,R(n*DOFs,n*DOFs),Q(n*DOFs,n*DOFs),N1(3),T1(3),T2(3)
!------------------------------------------------------------------------------

    DO i=1,MIN(n,SIZE(NodeIndexes))
      IF ( NodeIndexes(i)<=0 .OR. NodeIndexes(i)>SIZE(Normals,1) ) CYCLE

      R = 0.0d0
      DO j=1,n*DOFs
        R(j,j) = 1.0d0
      END DO

      N1 = Normals( NodeIndexes(i),: )

      SELECT CASE(DIM)
      CASE (2)
        R(DOFs*(i-1)+1,DOFs*(i-1)+1) =  N1(1)
        R(DOFs*(i-1)+1,DOFs*(i-1)+2) =  N1(2)

        R(DOFs*(i-1)+2,DOFs*(i-1)+1) = -N1(2)
        R(DOFs*(i-1)+2,DOFs*(i-1)+2) =  N1(1)
      CASE (3)
        T1 = Tangent1( NodeIndexes(i),: )
        T2 = Tangent2( NodeIndexes(i),: )

        R(DOFs*(i-1)+1,DOFs*(i-1)+1) = N1(1)
        R(DOFs*(i-1)+1,DOFs*(i-1)+2) = N1(2)
        R(DOFs*(i-1)+1,DOFs*(i-1)+3) = N1(3)

        R(DOFs*(i-1)+2,DOFs*(i-1)+1) = T1(1)
        R(DOFs*(i-1)+2,DOFs*(i-1)+2) = T1(2)
        R(DOFs*(i-1)+2,DOFs*(i-1)+3) = T1(3)

        R(DOFs*(i-1)+3,DOFs*(i-1)+1) = T2(1)
        R(DOFs*(i-1)+3,DOFs*(i-1)+2) = T2(2)
        R(DOFs*(i-1)+3,DOFs*(i-1)+3) = T2(3)
      END SELECT

      DO j=1,n*DOFs
        DO k=1,n*DOFs
          s = 0.0D0
          DO l=1,n*DOFs
            s = s + R(j,l) * Matrix(l,k)
          END DO
          Q(j,k) = s
        END DO
      END DO

      DO j=1,n*DOFs
        DO k=1,n*DOFs
          s = 0.0D0
          DO l=1,n*DOFs
            s = s + Q(j,l) * R(k,l)
          END DO
          Matrix(j,k) = s
        END DO
      END DO

      DO j=1,n*DOFs
        s = 0.0D0
        DO k=1,n*DOFs
          s = s + R(j,k) * Vector(k)
        END DO
        Q(j,1) = s
      END DO
      Vector(1:n*DOFs) = Q(:,1)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE RotateMatrix
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Given the normal return the tangent directions. The 
!> First tangent direction will always be on the xy-plane if
!> also the normal is in the xy-plane.
!------------------------------------------------------------------------------
  SUBROUTINE TangentDirections( Normal,Tangent1,Tangent2 )
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: Normal(3),Tangent1(3)
   REAL(KIND=dp), OPTIONAL :: Tangent2(3)
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: n1,n2,n3
!------------------------------------------------------------------------------
   n1 = ABS(Normal(1))
   n2 = ABS(Normal(2))
   n3 = ABS(Normal(3))

   IF( PRESENT( Tangent2 ) ) THEN   
     IF ( n1 <= n3 .AND. n2 <= n3 ) THEN
       Tangent1(1) =  0.0_dp
       Tangent1(2) = -Normal(3)
       Tangent1(3) =  Normal(2)
     ELSE
       Tangent1(1) = -Normal(2)
       Tangent1(2) =  Normal(1)
       Tangent1(3) =  0.0_dp
     END IF

     Tangent1 = Tangent1 / SQRT(SUM(Tangent1**2))
     Tangent2(1) = Normal(2)*Tangent1(3) - Normal(3)*Tangent1(2)
     Tangent2(2) = Normal(3)*Tangent1(1) - Normal(1)*Tangent1(3)
     Tangent2(3) = Normal(1)*Tangent1(2) - Normal(2)*Tangent1(1)
     Tangent2 = Tangent2 / SQRT(SUM(Tangent2**2))
   ELSE
     ! This is a 2D tangent only
     Tangent1(1) = Normal(2)
     Tangent1(2) = -Normal(1)
     Tangent1(3) = 0.0_dp
     Tangent1 = Tangent1 / SQRT(SUM(Tangent1**2))
   END IF
     
!------------------------------------------------------------------------------
 END SUBROUTINE TangentDirections
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>    Integrates a user-defined function over the specified bulk elements.
!> \deprecated Is this used for something?
!------------------------------------------------------------------------------
   FUNCTION VolumeIntegrate( Model, ElementList, IntegrandFunctionName ) &
       RESULT(Integral)
!------------------------------------------------------------------------------
   TYPE(Model_t) :: Model                     !< All model information (mesh, materials, BCs, etc...)
   INTEGER, DIMENSION(:) :: ElementList       !< List of elements that belong to the integration volume
   CHARACTER(LEN=*) :: IntegrandFunctionName  !< Name the function has in the .sif file
   REAL(KIND=dp) :: Integral                  !< The value of the volume integral
!------------------------------------------------------------------------------
     INTEGER :: n
     TYPE(Element_t), POINTER :: CurrentElement
     INTEGER, POINTER :: NodeIndexes(:)
     TYPE(Nodes_t)   :: ElementNodes

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     REAL(KIND=dp), DIMENSION(:), POINTER :: &
         U_Integ,V_Integ,W_Integ,S_Integ

     REAL(KIND=dp), DIMENSION(Model % MaxElementNodes) :: IntegrandFunction
     REAL(KIND=dp) :: s,ug,vg,wg
     REAL(KIND=dp) :: Basis(Model % MaxElementNodes)
     REAL(KIND=dp) :: dBasisdx(Model % MaxElementNodes,3),SqrtElementMetric
     REAL(KIND=dp) :: IntegrandAtGPt, dV
     INTEGER :: N_Integ, t, tg, i, istat
     LOGICAL :: stat

! Need MaxElementNodes only in allocation
     n = Model % MaxElementNodes
     ALLOCATE( ElementNodes % x( n ),   &
               ElementNodes % y( n ),   &
               ElementNodes % z( n ), STAT=istat )
     IF( istat /= 0 ) THEN
       CALL Fatal('VolumeIntegrate','Allocation error for ElementNodes')
     END IF

     Integral = 0.0d0

! Loop over all elements in the list
     DO i=1,SIZE(ElementList)

       t = ElementList(i)

       IF ( t < 1 .OR. t > Model % NumberOfBulkElements ) THEN
! do something
       END IF

       CurrentElement => Model % Elements(t)
       n = CurrentElement % TYPE % NumberOfNodes
       NodeIndexes => CurrentElement % NodeIndexes

!------------------------------------------------------------------------------
! Get element nodal coordinates
!------------------------------------------------------------------------------
       ElementNodes % x(1:n) = Model % Nodes % x(NodeIndexes)
       ElementNodes % y(1:n) = Model % Nodes % y(NodeIndexes)
       ElementNodes % z(1:n) = Model % Nodes % z(NodeIndexes)

! Read this from Simulation block in the .sif file
       IntegrandFunction(1:n) = ListGetReal( Model % Simulation, &
           IntegrandFunctionName, n, NodeIndexes )

!------------------------------------------------------------------------------
!    Gauss integration stuff
!------------------------------------------------------------------------------
       IntegStuff = GaussPoints( CurrentElement )
       U_Integ => IntegStuff % u
       V_Integ => IntegStuff % v
       W_Integ => IntegStuff % w
       S_Integ => IntegStuff % s
       N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
! Loop over Gauss integration points
!------------------------------------------------------------------------------
       DO tg=1,N_Integ

         ug = U_Integ(tg)
         vg = V_Integ(tg)
         wg = W_Integ(tg)

!------------------------------------------------------------------------------
! Need SqrtElementMetric and Basis at the integration point
!------------------------------------------------------------------------------
         stat = ElementInfo( CurrentElement,ElementNodes,ug,vg,wg, &
             SqrtElementMetric,Basis,dBasisdx )

         s = SqrtElementMetric * S_Integ(tg)

! Calculate the function to be integrated at the Gauss point
         IntegrandAtGPt = SUM( IntegrandFunction(1:n) * Basis )

! Use general coordinate system for dV
         dV = CoordinateSqrtMetric( SUM( ElementNodes % x(1:n) * Basis), &
             SUM( ElementNodes % y(1:n) * Basis), &
             SUM( ElementNodes % z(1:n) * Basis) )

         Integral = Integral + s*IntegrandAtGPt*dV

       END DO! of the Gauss integration points

     END DO! of the bulk elements

     DEALLOCATE( ElementNodes % x, &
         ElementNodes % y, &
         ElementNodes % z )
!------------------------------------------------------------------------------
   END FUNCTION VolumeIntegrate
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>    Integrates the normal component of a user-defined vector function
!>    over the specified boundary elements. 
!>  \deprecated Is this used for something?
!------------------------------------------------------------------------------
   FUNCTION FluxIntegrate( Model, ElementList, IntegrandFunctionName ) &
!------------------------------------------------------------------------------
       RESULT(Integral)
!------------------------------------------------------------------------------
   TYPE(Model_t) :: Model                     !< All model information (mesh, materials, BCs, etc...)
   INTEGER, DIMENSION(:) :: ElementList       !< List of elements that belong to the integration boundary
   CHARACTER(LEN=*) :: IntegrandFunctionName  !< Name the function has in the .sif file
   REAL(KIND=dp) :: Integral                  !< The value of the boundary integral
!------------------------------------------------------------------------------
     INTEGER :: n
     TYPE(Element_t), POINTER :: CurrentElement
     INTEGER, POINTER :: NodeIndexes(:)
     TYPE(Nodes_t)    :: ElementNodes

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     REAL(KIND=dp), DIMENSION(:), POINTER :: &
         U_Integ,V_Integ,W_Integ,S_Integ

     REAL(KIND=dp), DIMENSION(Model % MaxElementNodes,3) :: IntegrandFunction
!     REAL(KIND=dp), POINTER :: IntegrandFunction(:,:)
     CHARACTER(LEN=2) :: Component
     CHARACTER(LEN=MAX_NAME_LEN) :: IntegrandFunctionComponent
     REAL(KIND=dp) :: s,ug,vg,wg
     REAL(KIND=dp) :: Basis(Model % MaxElementNodes)
     REAL(KIND=dp) :: dBasisdx(Model % MaxElementNodes,3),SqrtElementMetric
     REAL(KIND=dp) :: Normal(3)
     REAL(KIND=dp) :: IntegrandAtGPt(3), FluxAtGPt, dS
     INTEGER :: N_Integ, t, tg, i, j, DIM, istat
     LOGICAL :: stat

     DIM = CoordinateSystemDimension()

! Need MaxElementNodes only in allocation
     n = Model % MaxElementNodes
     ALLOCATE( ElementNodes % x( n ),   &
               ElementNodes % y( n ),   &
               ElementNodes % z( n ), STAT=istat )
     IF( istat /= 0 ) THEN
       CALL Fatal('FluxIntegrate','Allocation error for ElementNodes')
     END IF

     Integral = 0.0d0

! Loop over all elements in the list
     DO i=1,SIZE(ElementList)

       t = ElementList(i)

       IF ( t < 1 .OR. t > Model % NumberOfBulkElements ) THEN
! do something
       END IF

       CurrentElement => Model % Elements(t)
       n = CurrentElement % TYPE % NumberOfNodes
       NodeIndexes => CurrentElement % NodeIndexes

       Model % CurrentElement => Model % Elements(t)
!------------------------------------------------------------------------------
! Get element nodal coordinates
!------------------------------------------------------------------------------
       ElementNodes % x(1:n) = Model % Nodes % x(NodeIndexes)
       ElementNodes % y(1:n) = Model % Nodes % y(NodeIndexes)
       ElementNodes % z(1:n) = Model % Nodes % z(NodeIndexes)

! Read the integrand from Simulation block in the .sif file
! It is assumed to be a contravariant vector, but
! ListGetRealArray doesn t exist, so we READ it component by component
! naming them with suffixes " 1" etc.
       DO j=1,DIM
         WRITE (Component, '(" ",I1.1)') j
         IntegrandFunctionComponent = IntegrandFunctionName(1: &
             LEN_TRIM(IntegrandFunctionName))
         IntegrandFunctionComponent(LEN_TRIM(IntegrandFunctionName)+1: &
             LEN_TRIM(IntegrandFunctionName)+2) = Component
         IntegrandFunction(1:n,j) = ListGetReal( Model % Simulation, &
          IntegrandFunctionComponent(1:LEN_TRIM(IntegrandFunctionComponent)), &
          n, NodeIndexes )
       END DO

!------------------------------------------------------------------------------
!    Gauss integration stuff
!------------------------------------------------------------------------------
       IntegStuff = GaussPoints( CurrentElement )
       U_Integ => IntegStuff % u
       V_Integ => IntegStuff % v
       W_Integ => IntegStuff % w
       S_Integ => IntegStuff % s
       N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
! Loop over Gauss integration points
!------------------------------------------------------------------------------
       DO tg=1,N_Integ

         ug = U_Integ(tg)
         vg = V_Integ(tg)
         wg = W_Integ(tg)

!------------------------------------------------------------------------------
! Need SqrtElementMetric and Basis at the integration point
!------------------------------------------------------------------------------
         stat = ElementInfo( CurrentElement,ElementNodes,ug,vg,wg, &
             SqrtElementMetric,Basis,dBasisdx )

! If we want to allow covariant integrand vectors given
! we would also need Metric, and read it as follows...
!     IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
!       x = SUM( nodes % x(1:n) * Basis )
!       y = SUM( nodes % y(1:n) * Basis )
!       z = SUM( nodes % z(1:n) * Basis )
!     END IF
!
!     CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )
!     s = SqrtMetric * SqrtElementMetric * S_Integ(t)
! And also FluxAtGpt...
! Now we won t need Metric, so get SqrtMetric later...

         s = SqrtElementMetric * S_Integ(tg)

! Get normal to the element boundary at the integration point
! N.B. NormalVector returns covariant normal vector
         Normal = NormalVector( CurrentElement,ElementNodes,ug,vg,.TRUE. )

! Calculate the contravariant vector function to be integrated
! at the Gauss point
         DO j=1,DIM
           IntegrandAtGPt(j) = SUM( IntegrandFunction(1:n,j) * Basis )
         END DO
! Calculate the normal component of the vector function
         FluxAtGPt = SUM( IntegrandAtGPt * Normal )

! Use general coordinate system for dS
! Would be included in s by SqrtMetric
         dS = CoordinateSqrtMetric( SUM( ElementNodes % x(1:n) * Basis), &
             SUM( ElementNodes % y(1:n) * Basis), &
             SUM( ElementNodes % z(1:n) * Basis) )

         Integral = Integral + s*FluxAtGPt*dS

       END DO! of the Gauss integration points

     END DO! of the boundary elements

     DEALLOCATE( ElementNodes % x, &
         ElementNodes % y, &
         ElementNodes % z )
!------------------------------------------------------------------------------
   END FUNCTION FluxIntegrate
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>    Integrates A user-defined vector function 
!>    over the specified boundary elements.
!> \deprecated Is this used for something?
!------------------------------------------------------------------------------
   FUNCTION SurfaceIntegrate( Model, ElementList, IntegrandFunctionName ) &
       RESULT(Integral)
!------------------------------------------------------------------------------
   TYPE(Model_t) :: Model                     !< All model information (mesh, materials, BCs, etc...)
   INTEGER, DIMENSION(:) :: ElementList       !< List of elements that belong to the integration boundary
   CHARACTER(LEN=*) :: IntegrandFunctionName  !< Name the function has in the .sif file
   REAL(KIND=dp), DIMENSION(3) :: Integral    !< The vector value of the integral
!------------------------------------------------------------------------------
     INTEGER :: n
     TYPE(Element_t), POINTER :: CurrentElement
     INTEGER, POINTER :: NodeIndexes(:)
     TYPE(Nodes_t)    :: ElementNodes

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     REAL(KIND=dp), DIMENSION(:), POINTER :: &
         U_Integ,V_Integ,W_Integ,S_Integ

     REAL(KIND=dp), DIMENSION(Model % MaxElementNodes,3) :: IntegrandFunction
!     REAL(KIND=dp), POINTER :: IntegrandFunction(:,:)
     CHARACTER(LEN=2) :: Component
     CHARACTER(LEN=MAX_NAME_LEN) :: IntegrandFunctionComponent
     REAL(KIND=dp) :: s,ug,vg,wg
     REAL(KIND=dp) :: Basis(Model % MaxElementNodes)
     REAL(KIND=dp) :: dBasisdx(Model % MaxElementNodes,3),SqrtElementMetric
!     REAL(KIND=dp) :: Normal(3)
     REAL(KIND=dp) :: IntegrandAtGPt(3), dS
     INTEGER :: N_Integ, t, tg, i, j, DIM, istat
     LOGICAL :: stat

     DIM = CoordinateSystemDimension()

! Need MaxElementNodes only in allocation
     n = Model % MaxElementNodes
     ALLOCATE( ElementNodes % x( n ),   &
               ElementNodes % y( n ),   &
               ElementNodes % z( n ), STAT=istat )
     IF( istat /= 0 ) THEN
       CALL Fatal('SurfaceIntegrate','Allocation error for ElementNodes')
     END IF

     Integral = 0.0d0

! Loop over all elements in the list
     DO i=1,SIZE(ElementList)

       t = ElementList(i)

       IF ( t < 1 .OR. t > Model % NumberOfBulkElements ) THEN
! do something
       END IF

       CurrentElement => Model % Elements(t)
       Model % CurrentElement => CurrentElement
       n = CurrentElement % TYPE % NumberOfNodes
       NodeIndexes => CurrentElement % NodeIndexes

!------------------------------------------------------------------------------
! Get element nodal coordinates
!------------------------------------------------------------------------------
       ElementNodes % x(1:n) = Model % Nodes % x(NodeIndexes)
       ElementNodes % y(1:n) = Model % Nodes % y(NodeIndexes)
       ElementNodes % z(1:n) = Model % Nodes % z(NodeIndexes)

! Read the integrand from Simulation block in the .sif file
! It is assumed to be a contravariant vector, but
! ListGetRealArray doesn t exist, so we READ it component by component
! naming them with suffixes " 1" etc.
       DO j=1,DIM
         WRITE (Component, '(" ",I1.1)') j
         IntegrandFunctionComponent = IntegrandFunctionName(1: &
             LEN_TRIM(IntegrandFunctionName))
         IntegrandFunctionComponent(LEN_TRIM(IntegrandFunctionName)+1: &
             LEN_TRIM(IntegrandFunctionName)+2) = Component
         IntegrandFunction(1:n,j) = ListGetReal( Model % Simulation, &
          IntegrandFunctionComponent(1:LEN_TRIM(IntegrandFunctionComponent)), &
          n, NodeIndexes )
       END DO

!------------------------------------------------------------------------------
!    Gauss integration stuff
!------------------------------------------------------------------------------
       IntegStuff = GaussPoints( CurrentElement )
       U_Integ => IntegStuff % u
       V_Integ => IntegStuff % v
       W_Integ => IntegStuff % w
       S_Integ => IntegStuff % s
       N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
! Loop over Gauss integration points
!------------------------------------------------------------------------------
       DO tg=1,N_Integ

         ug = U_Integ(tg)
         vg = V_Integ(tg)
         wg = W_Integ(tg)

!------------------------------------------------------------------------------
! Need SqrtElementMetric and Basis at the integration point
!------------------------------------------------------------------------------
         stat = ElementInfo( CurrentElement,ElementNodes,ug,vg,wg, &
             SqrtElementMetric,Basis,dBasisdx )

! If we want to allow covariant integrand vectors given
! we would also need Metric, and read it as follows...
!     IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
!       x = SUM( nodes % x(1:n) * Basis )
!       y = SUM( nodes % y(1:n) * Basis )
!       z = SUM( nodes % z(1:n) * Basis )
!     END IF
!
!     CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )
!     s = SqrtMetric * SqrtElementMetric * S_Integ(t)
! Now we won t need Metric, so get SqrtMetric later...

         s = SqrtElementMetric * S_Integ(tg)

! If you need normal directly at the integration point
!         Normal = NormalVector( CurrentElement,ElementNodes,ug,vg,.TRUE. )

! Calculate the contravariant vector function to be integrated
! at the Gauss point
         DO j=1,DIM
           IntegrandAtGPt(j) = SUM( IntegrandFunction(1:n,j) * Basis )
         END DO

! Use general coordinate system for dS
! Would be included in s by SqrtMetric
         dS = CoordinateSqrtMetric( SUM( ElementNodes % x(1:n) * Basis), &
             SUM( ElementNodes % y(1:n) * Basis), &
             SUM( ElementNodes % z(1:n) * Basis) )

         DO j=1,DIM
           Integral(j) = Integral(j) + s*IntegrandAtGPt(j)*dS
         END DO

       END DO! of the Gauss integration points

     END DO! of the boundary elements

     DEALLOCATE( ElementNodes % x, &
         ElementNodes % y, &
         ElementNodes % z )
!------------------------------------------------------------------------------
   END FUNCTION SurfaceIntegrate
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>    Integrates the normal component of a user-defined vector function
!>    over a specified line element.
!> \deprecated Is this used for something?
!------------------------------------------------------------------------------
   FUNCTION LineIntegrate( Model, LineElement, LineElementNodes, &
       IntegrandFunctionName, QuadrantTreeExists, RootQuadrant ) &
       RESULT(Integral)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!
!  TYPE(Model_t), POINTER :: Model,  
!     INPUT: All model information (mesh, materials, BCs, etc...)
!
!  TYPE(Element_t) :: LineElement
!     INPUT: Line element that belongs to the line of integration
!
!  REAL(KIND=dp), DIMENSION(LineElement % Type % NumberOfNodes,3) ::
!     LineElementNodes
!     INPUT: List of nodal point coordinates
!
!  CHARACTER(LEN=MAX_NAME_LEN) :: IntegrandFunctionName
!     INPUT: Name the function has in the .sif file or somewhere else
!
!  LOGICAL :: QuadrantTreeExists
!     INPUT: QuadrantTree has been built, use it in element search
!
!  TYPE(Quadrant_t), POINTER :: RootQuadrant
!     OUTPUT: Quadrant tree structure root
!
!  FUNCTION RETURN VALUE:
!    REAL(KIND=dp) :: Integral
!     The value of the flux integral
!      
!------------------------------------------------------------------------------
   TYPE(Model_t) :: Model
! TARGET only for CurrentElement
   TYPE(Element_t), TARGET :: LineElement
   REAL(KIND=dp), DIMENSION(:,:), TARGET CONTIG :: LineElementNodes
   CHARACTER(LEN=*) :: IntegrandFunctionName
   REAL(KIND=dp) :: Integral
   LOGICAL :: QuadrantTreeExists
   TYPE(Quadrant_t), POINTER :: RootQuadrant
!------------------------------------------------------------------------------
     INTEGER :: n
! Only one element at a time, need CurrentElement only for NormalVector!
     TYPE(Element_t), POINTER :: CurrentElement
! LineElement nodes don t belong to global node structure
! Need the structure for the function calls
     TYPE(Nodes_t) :: ElementNodes

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     REAL(KIND=dp), DIMENSION(:), POINTER :: &
         U_Integ,V_Integ,W_Integ,S_Integ

! IntegrandFunction at the bulk element nodal points
     REAL(KIND=dp), DIMENSION(Model % MaxElementNodes,3) :: IntegrandAtNodes
! IntegrandFunction at the Gauss points
     REAL(KIND=dp), DIMENSION(LineElement % TYPE % GaussPoints,3) :: IntegrandFunction
     CHARACTER(LEN=2) :: Component
     CHARACTER(LEN=MAX_NAME_LEN) :: IntegrandFunctionComponent
     REAL(KIND=dp) :: s,ug,vg,wg
     REAL(KIND=dp) :: Basis(LineElement % TYPE % NumberOfNodes)
     REAL(KIND=dp) :: dBasisdx(LineElement % TYPE % NumberOfNodes,3),SqrtElementMetric
     REAL(KIND=dp) :: Normal(3)
! IntegrandFunction already at GPts
     REAL(KIND=dp) :: FluxAtGPt, dS
     INTEGER :: N_Integ, t, tg, i, j, DIM
     LOGICAL :: stat
! Search for the bulk element each Gauss point belongs to
     TYPE(Element_t), POINTER :: BulkElement
     TYPE(Nodes_t) :: BulkElementNodes
     INTEGER, POINTER :: NodeIndexes(:)
     REAL(KIND=dp), DIMENSION(3) :: LocalCoordsInBulkElem
     REAL(KIND=dp), DIMENSION(3) :: Point
     INTEGER :: nBulk, maxlevel=10, k, Quadrant
     TYPE(Quadrant_t), POINTER :: LeafQuadrant

     DIM = CoordinateSystemDimension()

     n = LineElement % TYPE % NumberOfNodes

     Integral = 0.0d0

!------------------------------------------------------------------------------
! Get element nodal coordinates
!------------------------------------------------------------------------------
! Move from LineElementNodes to Nodes_t structure
     ElementNodes % x => LineElementNodes(1:n,1)
     ElementNodes % y => LineElementNodes(1:n,2)
     ElementNodes % z => LineElementNodes(1:n,3)

!------------------------------------------------------------------------------
!    Gauss integration stuff
!------------------------------------------------------------------------------
     IntegStuff = GaussPoints( LineElement )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
! Loop over Gauss integration points
!------------------------------------------------------------------------------
     DO tg=1,N_Integ

       ug = U_Integ(tg)
       vg = V_Integ(tg)
       wg = W_Integ(tg)

!------------------------------------------------------------------------------
! Need SqrtElementMetric and Basis at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( LineElement,ElementNodes,ug,vg,wg, &
           SqrtElementMetric,Basis,dBasisdx )

! If we want to allow covariant integrand vectors given... (see FluxIntegrate)

       s = SqrtElementMetric * S_Integ(tg)

! Find in which bulk element the Gauss point belongs to
       IF ( QuadrantTreeExists ) THEN
! Find the last existing quadrant that the point belongs to
         Point = [ SUM( ElementNodes % x(1:n) * Basis), &
             SUM( ElementNodes % y(1:n) * Basis), &
             SUM( ElementNodes % z(1:n) * Basis) ]
!         PRINT*,'Point:', Point
         CALL FindLeafElements(Point, DIM, RootQuadrant, LeafQuadrant)         
!         PRINT*,'Elems in LeafQuadrant',LeafQuadrant % NElemsInQuadrant
! Go through the bulk elements in the last ChildQuadrant only
         nBulk = Model % MaxElementNodes
         ALLOCATE( BulkElementNodes % x( nBulk ),   &
             BulkElementNodes % y( nBulk ),   &
             BulkElementNodes % z( nBulk ) )
!         PRINT*,'Elements:', LeafQuadrant % Elements
         DO k=1, LeafQuadrant % NElemsInQuadrant
           BulkElement => Model % Elements( &
               LeafQuadrant % Elements(k) )
           nBulk = BulkElement % TYPE % NumberOfNodes
           NodeIndexes => BulkElement % NodeIndexes
           BulkElementNodes % x(1:nBulk) = Model % Nodes % x(NodeIndexes)
           BulkElementNodes % y(1:nBulk) = Model % Nodes % y(NodeIndexes)
           BulkElementNodes % z(1:nBulk) = Model % Nodes % z(NodeIndexes)
           IF ( PointInElement( BulkElement,BulkElementNodes, &
               [ SUM( ElementNodes % x(1:n) * Basis), &
               SUM( ElementNodes % y(1:n) * Basis), &
               SUM( ElementNodes % z(1:n) * Basis) ], &
               LocalCoordsInBulkElem) ) EXIT
         END DO
!         PRINT*,'Point in Element: ', LeafQuadrant % Elements(k)
        ELSE
! Go through all BulkElements
! Need MaxElementNodes only in allocation
         nBulk = Model % MaxElementNodes
         ALLOCATE( BulkElementNodes % x( nBulk ),   &
             BulkElementNodes % y( nBulk ),   &
             BulkElementNodes % z( nBulk ) )
         DO k=1,Model % NumberOfBulkElements
           BulkElement => Model % Elements(k)
           nBulk = BulkElement % TYPE % NumberOfNodes
           NodeIndexes => BulkElement % NodeIndexes
           BulkElementNodes % x(1:nBulk) = Model % Nodes % x(NodeIndexes)
           BulkElementNodes % y(1:nBulk) = Model % Nodes % y(NodeIndexes)
           BulkElementNodes % z(1:nBulk) = Model % Nodes % z(NodeIndexes)
           IF ( PointInElement(BulkElement,BulkElementNodes, &
               [ SUM(ElementNodes % x(1:n) * Basis), &
                 SUM(ElementNodes % y(1:n) * Basis), &
                 SUM(ElementNodes % z(1:n) * Basis) ], &
               LocalCoordsInBulkElem) ) EXIT
         END DO
!         PRINT*,'Point in Element: ', k
       END IF
! Calculate value of the function in the bulk element
! Read the integrand from Simulation block in the .sif file
! It is assumed to be a contravariant vector, but
! ListGetRealArray doesn t exist, so we read it component by component
! naming them with suffixes " 1" etc.

      DO j=1,DIM
        WRITE (Component, '(" ",I1.1)') j
        IntegrandFunctionComponent = IntegrandFunctionName(1: &
            LEN_TRIM(IntegrandFunctionName))
        IntegrandFunctionComponent(LEN_TRIM(IntegrandFunctionName)+1: &
            LEN_TRIM(IntegrandFunctionName)+2) = Component
        IntegrandAtNodes(1:nBulk,j) = ListGetReal( Model % Simulation, &
            IntegrandFunctionComponent(1:LEN_TRIM(IntegrandFunctionComponent)), &
            nBulk, NodeIndexes )
      END DO

      DO j=1,DIM
        IntegrandFunction(tg,j) = InterpolateInElement( BulkElement, &
            IntegrandAtNodes(1:nBulk,j),LocalCoordsInBulkElem(1), &
            LocalCoordsInBulkElem(2),LocalCoordsInBulkElem(3) )
      END DO

      DEALLOCATE( BulkElementNodes % x, &
          BulkElementNodes % y, &
          BulkElementNodes % z )

! Get normal to the element boundary at the integration point
! N.B. NormalVector returns covariant normal vector
! NormalVector defined weirdly, doesn t accept LineElement as an argument,
! but wants a pointer
         CurrentElement => LineElement
         Normal = NormalVector( CurrentElement,ElementNodes,ug,vg,.FALSE. )
! Might be consistently in wrong direction, since no check
!         Normal = NormalVector( CurrentElement,ElementNodes,ug,vg,.TRUE. )

! Contravariant vector function to be integrated is already
! at the Gauss point

! Calculate the normal component of the vector function
         FluxAtGPt = SUM( IntegrandFunction(tg,1:DIM) * Normal(1:DIM) )

! Use general coordinate system for dS
! Would be included in s by SqrtMetric
         dS = CoordinateSqrtMetric( SUM( ElementNodes % x(1:n) * Basis), &
             SUM( ElementNodes % y(1:n) * Basis), &
             SUM( ElementNodes % z(1:n) * Basis) )

         Integral = Integral + s*FluxAtGPt*dS

       END DO! of the Gauss integration points

!------------------------------------------------------------------------------
   END FUNCTION LineIntegrate
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   FUNCTION ElementArea( Mesh,Element,N ) RESULT(A)
!------------------------------------------------------------------------------
     TYPE(Mesh_t), POINTER :: Mesh
     INTEGER :: N
     TYPE(Element_t) :: Element
!------------------------------------------------------------------------------

     REAL(KIND=dp), TARGET :: NX(N),NY(N),NZ(N)

     REAL(KIND=dp) :: A,R1,R2,Z1,Z2,S,U,V,W,X,Y,Z

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     INTEGER :: N_Integ,t

     REAL(KIND=dp) :: Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3), &
              SqrtMetric,SqrtElementMetric

     TYPE(Nodes_t) :: Nodes

     LOGICAL :: stat

     REAL(KIND=dp) :: Basis(n)
     REAL(KIND=dp) :: dBasisdx(n,3)

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
!------------------------------------------------------------------------------
 
#if 0
     IF ( ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
            CurrentCoordinateSystem() == CylindricSymmetric ) .AND. Element % TYPE % ELementCode / 100 == 2 ) THEN
       R1 = Mesh % Nodes % x(Element % NodeIndexes(1))
       R2 = Mesh % Nodes % x(Element % NodeIndexes(2))

       Z1 = Mesh % Nodes % y(Element % NodeIndexes(1))
       Z2 = Mesh % Nodes % y(Element % NodeIndexes(2))

       A = PI*ABS(R1+R2)*SQRT((Z1-Z2)*(Z1-Z2)+(R1-R2)*(R1-R2))
     ELSE 
#endif
       Nodes % x => NX
       Nodes % y => NY
       Nodes % z => NZ

       Nodes % x = Mesh % Nodes % x(Element % NodeIndexes)
       Nodes % y = Mesh % Nodes % y(Element % NodeIndexes)
       Nodes % z = Mesh % Nodes % z(Element % NodeIndexes)

       IntegStuff = GaussPoints( element )
       U_Integ => IntegStuff % u
       V_Integ => IntegStuff % v
       W_Integ => IntegStuff % w
       S_Integ => IntegStuff % s
       N_Integ  = IntegStuff % n
!
!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------
!
       A = 0.0
       DO t=1,N_Integ
!
!        Integration stuff
!
         u = U_Integ(t)
         v = V_Integ(t)
         w = W_Integ(t)
!
!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
         stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
                    Basis,dBasisdx )
!------------------------------------------------------------------------------
!        Coordinatesystem dependent info
!------------------------------------------------------------------------------
         IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
           X = SUM( Nodes % x(1:n)*Basis )
           Y = SUM( Nodes % y(1:n)*Basis )
           Z = SUM( Nodes % z(1:n)*Basis )

           SqrtMetric = CoordinateSqrtMetric( x,y,z )
           A =  A + SqrtMetric * SqrtElementMetric * S_Integ(t)
         ELSE
           A =  A + SqrtElementMetric * S_Integ(t)
         END IF
       END DO
#if 0
     END IF
#endif

   END FUNCTION ElementArea
!------------------------------------------------------------------------------


   !------------------------------------------------------------------------------
   !> If element has two of the same indexes regard the element as degenerate.
   !------------------------------------------------------------------------------
   FUNCTION DegenerateElement( Element ) RESULT ( Stat ) 
     TYPE(Element_t), POINTER :: Element
     LOGICAL Stat

     INTEGER :: i,n
     INTEGER, POINTER :: Indexes(:)
     
     Stat = .FALSE.

     n = Element % TYPE % NumberOfNodes
     Indexes => Element % NodeIndexes
     
     DO i = 1, n
       IF( ANY( Indexes(i+1:n) == Indexes(i) ) ) THEN
         Stat = .TRUE.           
         EXIT
       END IF
     END DO
     
   END FUNCTION DegenerateElement

   
END MODULE ElementUtils

!> \} ElmerLib
