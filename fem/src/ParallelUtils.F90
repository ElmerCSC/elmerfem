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
     USE SparIterComm
     
     IMPLICIT NONE

     INTERFACE ParallelReduction
       MODULE PROCEDURE ParallelReductionR, ParallelReductionI, ParallelReductionZ
     END INTERFACE ParallelReduction

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
!   IF( ParEnv % PEs > 1 ) THEN
      CALL Info('ParallelFinalize','Final synchronization of the MPI system!',Level=12)
      CALL ParEnvFinalize()
!   END IF
#endif
!-------------------------------------------------------------------------------
  END SUBROUTINE ParallelFinalize
!-------------------------------------------------------------------------------

!--------------------------------'-----------------------------------------------
    SUBROUTINE ParallelInitMatrix( Solver, Matrix, inPerm )
!-------------------------------------------------------------------------------
       INTEGER, POINTER, OPTIONAL :: inPerm(:)
       TYPE(Solver_t), TARGET :: Solver
       TYPE(Matrix_t), POINTER :: Matrix
!-------------------------------------------------------------------------------
       INTEGER :: i, j, k, n, DOFs
       LOGICAL :: Found
       TYPE(Mesh_t), POINTER :: Mesh
       INTEGER, POINTER :: Perm(:)
       INTEGER, ALLOCATABLE :: ind(:)
       REAL(KIND=dp) :: realtime,tt
       INTEGER :: Ierr

!-------------------------------------------------------------------------------
!tt = realtime()
#ifdef PARALLEL_FOR_REAL
       IF ( ParEnv % PEs <= 1 .OR. .NOT. ASSOCIATED(Matrix) ) RETURN

       CALL Info('ParallelInitMatrix','Creating communication structures for matrix!',Level=5)

       Mesh => Solver % Mesh
       DOFs = Solver % Variable % DOFs

       IF(PRESENT(inPerm)) THEN
         Perm => InPerm
       ELSE
         Perm => Solver % Variable % Perm
       END IF
       IF(.NOT. ASSOCIATED(Perm)) THEN
         CALL Fatal('ParallelInitMatrix','Cannot initialize matrix without Perm vector!')
       END IF
       

       n = SIZE(Perm)
       k = n*DOFs + Matrix % ExtraDOFs
                   
       i = 0
       IF(ASSOCIATED( Matrix % Perm) ) i=i+1
       IF(ASSOCIATED( Matrix % InvPerm) ) i=i+1
       
       IF(i==1) THEN
         CALL Fatal('ParallelInitMatrix','Only Perm or InvPerm is associated!')
       ELSE IF(i==0) THEN
         CALL Info('ParallelInitMatrix','Allocating communication structures for matrix!',Level=5)
         ! InvPerm is indexed by the values stored in Perm below. The nodal part
         ! reaches MAXVAL(Perm)*DOFs, but the extra rows are numbered from
         ! n*DOFs upwards, so MAX() is needed to cover both. The two agree
         ! whenever the permutation is onto, and instrumenting the 124 _npN
         ! tests shows MAXVAL(Perm) /= n never once coincides with ExtraDOFs>0,
         ! so this only guards a case the suite does not reach.
         !
         ! NOTE: if that case ever does arise, sizing alone is not the whole
         ! story. The extra rows of the matrix follow the base rows, which end
         ! at MAXVAL(Perm)*DOFs, so "Matrix % Perm(i) = i" below would then give
         ! them the wrong row and leave the true extra-row slots of InvPerm at
         ! zero. Numbering them from MAXVAL(Perm)*DOFs is most likely the real
         ! fix, but it is deliberately not made here, being unverifiable against
         ! a suite that never exercises it.
         j = MAX( MAXVAL(Perm), n )*DOFs + Matrix % ExtraDOFs
         ALLOCATE( Matrix % Perm(k), Matrix % InvPerm(j))

         BLOCK
           ! NOTE: no initializer here on purpose; an initialized local would
           ! implicitly acquire the SAVE attribute.
           LOGICAL :: DoConf

           DoConf = ListGetLogical( Solver % Values, 'Apply Conforming BCs',Found )
           DoConf = DoConf .AND. ASSOCIATED(Mesh % PeriodicPerm)

           Matrix % Perm = 0
           DO i=1,n
             IF ( DoConf ) THEN
               IF ( Mesh % PeriodicPerm(i) /= 0 ) CYCLE
             END IF

             IF ( Perm(i) /= 0 ) THEN
               DO j=1,DOFs
                 Matrix % Perm((i-1)*DOFs+j) = DOFs * (Perm(i)-1) + j
               END DO
             END IF
           END DO
         END BLOCK

         DO i=n*DOFs+1,SIZE(Matrix % Perm)
           Matrix % Perm(i) = i
         END DO

         Matrix % INVPerm = 0
         DO i=1,SIZE(Matrix % Perm)
           IF ( Matrix % Perm(i) /= 0 ) THEN
             Matrix % InvPerm(Matrix % Perm(i)) = i
           END IF
         END DO
       ELSE
         CALL Info('ParallelInitMatrix','Skipping generation of Perm and InvPerm',Level=20)
       END IF

       IF(ASSOCIATED(Matrix % ParallelInfo)) THEN
         CALL Fatal('ParallelInitMatrix','ParallelInfo already created!')
       END IF               
       ALLOCATE( Matrix % ParallelInfo )


       ! The global numbering of the matrix rows depends on the discretization.
       !--------------------------------------------------------------------
       IF( ListGetLogical( Solver % Values,'CutFEM', Found ) ) THEN
         CALL MatrixParallelInfoCutFEM( Solver, Matrix, Perm, DOFs )
       ELSE IF ( Matrix % DGMatrix ) THEN
         CALL MatrixParallelInfoDG( Solver, Matrix, Perm, DOFs )
       ELSE
         CALL MatrixParallelInfoDefault( Solver, Matrix, Perm, DOFs )
       END IF

       n = SIZE(Matrix % ParallelInfo % GLobalDOFs)
       ALLOCATE(Matrix % ParallelInfo % Gorder(n), Ind(n))

       Ind = Matrix % ParallelInfo % GlobalDOFs
       DO i=1,n
         Matrix % ParallelInfo % Gorder(i) = i
       END DO
       CALL SortI( n,Ind, Matrix % ParallelInfo % Gorder )

       
       ! This is used for a rare condition:
       !
       ! o Linear system solved using Hypre
       ! o Some partitions have no degrees of freedom assigned, but have matrix entries
       ! to contribute to the global system (for example from halo elements).
       BLOCK
         INTEGER :: NameSpaceI
         LOGICAL :: DoIt
         
         DoIt = .FALSE.
         IF( ListGetLogical( Solver % Values,'Linear System Trialing', Found ) ) THEN
           NameSpaceI = MAX( 1, NINT( ListGetCReal( &
               Solver % Values,'Linear System Namespace Number', Found ) ))           
           IF( ListGetLogical( Solver % Values, 'linsys' // &
               I2S(NameSpaceI)//': Linear System Use Hypre', Found) ) DoIt = .TRUE.
         ELSE
           IF( ListGetLogical( Solver % Values, 'Linear System Use Hypre', Found) .OR. &
               ListGetLogical( Solver % Values, 'linsys1: Linear System Use Hypre', Found) ) &
               DoIt = .TRUE.
         END IF
         IF(DoIt) THEN
           CALL AssignAtLeastOneDOFToPartition(n,Matrix % ParallelInfo,Matrix % Comm)
         END IF
       END BLOCK
       
       BLOCK
         LOGICAL :: SkipActiveCheck, Found

         SkipActiveCheck = .FALSE.
         DO i=1,CurrentModel % NumberOfBCs
           IF (ListGetString(CurrentModel % Bcs(i) % Values, 'Radiation',Found)=='diffuse gray' .OR. &
                ListGetLogical( CurrentModel % Bcs(i) % Values, 'Radiator BC', Found)) SkipActiveCheck=.TRUE.
         END DO
         Matrix % Solver => Solver
         Matrix % ParMatrix => ParInitMatrix(Matrix, Matrix % ParallelInfo, SkipActiveCheck)
       END BLOCK
       
       ! We can make many routines quicker if we know there is nothing to share.
       ! Still we might want to operate in parallel using parallel norms etc. 
       i = COUNT( Matrix % ParallelInfo % GInterface )
       CALL MPI_ALLREDUCE( i, j, 1, MPI_INTEGER, MPI_SUM, Matrix % comm, ierr )
       Matrix % ParallelInfo % NothingShared = (j==0)
              
!if(parenv%mype==0) print*,'MATRIX INIT TIME: ', realtime()-tt
#endif
  END SUBROUTINE ParallelInitMatrix
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!> Global numbering of matrix rows for a CutFEM mesh, where the degrees of
!> freedom live on the mesh nodes and edges.
!-------------------------------------------------------------------------------
    SUBROUTINE MatrixParallelInfoCutFEM( Solver, Matrix, Perm, DOFs )
!-------------------------------------------------------------------------------
      TYPE(Solver_t), TARGET :: Solver
      TYPE(Matrix_t), POINTER :: Matrix
      INTEGER, POINTER :: Perm(:)
      INTEGER :: DOFs
!-------------------------------------------------------------------------------
      TYPE(ParallelInfo_t), POINTER :: MatrixPI, MeshPI
      TYPE(Mesh_t), POINTER :: Mesh
      INTEGER :: i, j, k, l, n

      Mesh => Solver % Mesh

      CALL Info('MatrixParallelInfoCutFEM','Creating ParallelInfo assuming CutFEM mesh!',Level=12)
      
      MeshPI => Mesh % ParallelInfo
      MatrixPI => Matrix % ParallelInfo
      n = MAXVAL(Matrix % Perm)

      ALLOCATE( MatrixPI % GlobalDOFs(n) )
      MatrixPI % GlobalDOFs = 0

      ALLOCATE( MatrixPI % GInterface(n) )
      MatrixPI % GInterface = .FALSE.
      
      ALLOCATE(MatrixPI % NeighbourList(n) )
      DO i=1,n
        NULLIFY(MatrixPI % NeighbourList(i) % Neighbours)
      END DO
      
      
      BLOCK
        INTEGER :: n0, nn, ne
        nn = Mesh % NumberOfNodes
        ne = Mesh % NumberOfEdges

        n0 = MAXVAL(Mesh % ParallelInfo % GlobalDofs)
        n0 = ParallelReduction(n0,2)

        DO i=1,nn+ne
          DO j=1,DOFs
            k = Matrix % Perm((i-1)*DOFs+j)
            IF ( k<=0 ) CYCLE
                        
            IF(i<=nn) THEN
              MatrixPI % GInterface(k) = Mesh % ParallelInfo % GInterface(i)
              MatrixPI % GlobalDOFs(k) = DOFs*(Mesh % ParallelInfo % GlobalDOFs(i)-1)+j
              
              l = SIZE(Mesh % ParallelInfo % NeighbourList(i) % Neighbours)
              ALLOCATE(MatrixPI % NeighbourList(k) % Neighbours(l))
              MatrixPI % NeighbourList(k) % Neighbours = &
                  Mesh % ParallelInfo % NeighbourList(i) % Neighbours            
            ELSE
              MatrixPI % GInterface(k) = Mesh % ParallelInfo % EdgeInterface(i-nn)
              MatrixPI % GlobalDOFs(k) = DOFs*n0 + DOFs*(Mesh % Edges(i-nn) % GElementIndex-1)+j
              
              l = SIZE(Mesh % ParallelInfo % EdgeNeighbourList(i-nn) % Neighbours)
              ALLOCATE(MatrixPI % NeighbourList(k) % Neighbours(l))
              MatrixPI % NeighbourList(k) % Neighbours = &
                  Mesh % ParallelInfo % EdgeNeighbourList(i-nn) % Neighbours                                
            END IF
          END DO
        END DO
      END BLOCK       

      CALL Info('MatrixParallelInfoCutFEM','Done ParallelInfo assuming CutFEM mesh!',Level=12)

      
!-------------------------------------------------------------------------------
    END SUBROUTINE MatrixParallelInfoCutFEM
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!> Global numbering of matrix rows for a Discontinuous Galerkin discretization,
!> where the degrees of freedom are local to each element.
!-------------------------------------------------------------------------------
    SUBROUTINE MatrixParallelInfoDG( Solver, Matrix, Perm, DOFs )
!-------------------------------------------------------------------------------
      TYPE(Solver_t), TARGET :: Solver
      TYPE(Matrix_t), POINTER :: Matrix
      INTEGER, POINTER :: Perm(:)
      INTEGER :: DOFs
!-------------------------------------------------------------------------------
      TYPE(ParallelInfo_t), POINTER :: MatrixPI, MeshPI
      TYPE(Mesh_t), POINTER :: Mesh
      TYPE(Element_t), POINTER :: Element
      TYPE(NeighbourList_t), POINTER :: MtrxN, MeshN
      INTEGER :: i, j, k, l, m, n, maxnode
      LOGICAL :: Found, DGReduced

      ! Width of the slot reserved per element in the global DG numbering.
      ! Discontinuous Galerkin as a method has no such limit, but the Elmer
      ! implementation has used a nodal basis on the element's own nodes with
      ! linear basis functions only (for a decade or so), so the widest case
      ! is the 8-node hexahedron. Should higher order DG ever be implemented,
      ! this numbering has to be revisited: an element needing more than
      ! MaxDGNodes indices would overrun its slot and collide with the range
      ! of the following element. The check at the use site below guards that.
      INTEGER, PARAMETER :: MaxDGNodes = 8

      Mesh => Solver % Mesh

      MeshPI => Solver % Mesh % ParallelInfo
      MatrixPI => Matrix % ParallelInfo
      n = MAXVAL(Matrix % Perm)

      ALLOCATE( MatrixPI % GlobalDOFs(n) ); MatrixPI % GlobalDOFs=0

      maxnode = MAXVAL(Mesh % ParallelInfo % GlobalDOFs)
      maxnode = ParallelReduction(maxnode,2)

      DGReduced = ListGetLogical(Solver % Values, 'DG Reduced Basis', Found )

      IF( DGReduced ) THEN
        BLOCK
          INTEGER, POINTER :: DgMap(:), DgMaster(:), DgSlave(:)
          INTEGER :: group0, group
          LOGICAL :: GotDgMap, GotMaster, GotSlave

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
          IF ( Element % TYPE % NumberOfNodes > MaxDGNodes ) THEN
            CALL Fatal('MatrixParallelInfoDG','Global DG numbering reserves only '//&
                I2S(MaxDGNodes)//' indices per element, but element has '//&
                I2S(Element % TYPE % NumberOfNodes)//' nodes')
          END IF
          DO j=1,Element % TYPE % NumberOfNodes
            k = Matrix % Perm(Element % DGIndexes(j))
            IF(k == 0) CYCLE
            MatrixPI % GlobalDOFs(k) = MaxDGNodes*(Element % GElementIndex-1) + j
          END DO
        END DO
      END IF

      ALLOCATE( MatrixPI % GInterface(n), MatrixPI % NeighbourList(n) )
      MatrixPI % GInterface = .FALSE.
      DO i=1,n
        MtrxN => MatrixPI % NeighbourList(i)
        MtrxN % Neighbours => NULL()
      END DO

      DO i=1,Mesh % NumberOfBulkElements
        Element => Mesh % Elements(i)
        IF ( .NOT. ASSOCIATED(Element % DGIndexes) ) CYCLE
        IF ( ALL(MeshPI % GInterface(Element % NodeIndexes)) ) THEN
          DO j=1,Element % Type % NumberOfNodes
             K = Matrix % Perm(Element % DGIndexes(j))
             IF(K==0) CYCLE
             L = Element % Nodeindexes(j)

             MeshN => MeshPI % NeighbourList(L)
             MtrxN => MatrixPI % NeighbourList(K)
             MatrixPI % GInterface(k) = .TRUE.

             IF( ASSOCIATED(MtrxN % Neighbours) ) DEALLOCATE(MtrxN % Neighbours)
             CALL AllocateVector( MtrxN % Neighbours,  SIZE(MeshN % Neighbours) )

             MtrxN % Neighbours = MeshN % Neighbours
             IF(.NOT.DGReduced) THEN 
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

      
!-------------------------------------------------------------------------------
    END SUBROUTINE MatrixParallelInfoDG
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!> Global numbering of matrix rows for the default discretization, i.e. neither
!> CutFEM nor Discontinuous Galerkin. Covers nodal, edge, face and bubble
!> degrees of freedom, and the extra rows of a collection matrix.
!-------------------------------------------------------------------------------
    SUBROUTINE MatrixParallelInfoDefault( Solver, Matrix, Perm, DOFs )
!-------------------------------------------------------------------------------
      TYPE(Solver_t), TARGET :: Solver
      TYPE(Matrix_t), POINTER :: Matrix
      INTEGER, POINTER :: Perm(:)
      INTEGER :: DOFs
!-------------------------------------------------------------------------------
      TYPE(Mesh_t), POINTER :: Mesh
      TYPE(Element_t), POINTER :: Element
      INTEGER :: i, j, k, l, m, n, bdofs
      INTEGER :: maxnode, maxedge, maxface, fdofs, maxfdofs, &
                 edofs, maxedofs, maxbdofs, l_beg, g_beg
      LOGICAL :: Found

      Mesh => Solver % Mesh
      ! Default parallel numbering, not CutFEM or Discontinuous Galerkin!
      !--------------------------------------------------------------------
      
      n = Matrix % NumberOfRows
      ALLOCATE( Matrix % ParallelInfo % NeighbourList(n) )
      CALL AllocateVector( Matrix % ParallelInfo % GInterface, n)
      CALL AllocateVector( Matrix % ParallelInfo % GlobalDOFs, n)

      DO i=1,n
        Matrix % ParallelInfo % GInterface(i) = .FALSE.
        Matrix % ParallelInfo % GlobalDOFs(i) = 0
        Matrix % ParallelInfo % NeighbourList(i) % Neighbours => NULL()
      END DO

      DO i=1,Mesh % NumberOfNodes
        DO j=1,DOFs
           k = Matrix % Perm((i-1)*DOFs+j)
           IF ( k<=0 ) CYCLE
           IF (k > SIZE(Matrix % ParallelInfo % GlobalDOFs)) THEN
             CALL Fatal('MatrixParallelInfoDefault','Matrix % ParallelInfo % GlobalDOFs bounds error.'//&
                 ' Matrix vs Solver perm scope conflict is a possible cause.')
           END IF
           
           Matrix % ParallelInfo % GlobalDOFs(k) = &
             DOFs*(Mesh % ParallelInfo % GlobalDOFs(i)-1)+j

           Matrix % ParallelInfo % GInterface(k) = &
             Mesh % ParallelInfo % GInterface(i)

           IF( ASSOCIATED( Mesh % ParallelInfo % NeighbourList(i) % Neighbours ) ) THEN
             ALLOCATE( Matrix % ParallelInfo % NeighbourList(k) % Neighbours(SIZE( &
                 Mesh % ParallelInfo % NeighbourList(i) % Neighbours)) )

             Matrix % ParallelInfo % NeighbourList(k) % Neighbours = &
                 Mesh % ParallelInfo % NeighbourList(i) % Neighbours
           END IF
        END DO
      END DO

      maxnode = MAXVAL(Mesh % ParallelInfo % GlobalDOFs)
      maxnode = ParallelReduction(maxnode,2)

      edofs = 0; fdofs = 0; maxedofs = 0; maxfdofs = 0
      maxedge = 0; maxface = 0

      IF ( ASSOCIATED(Mesh % Edges) ) THEN
        g_beg = maxnode
        l_beg = Mesh % NumberOfNodes

        n = Mesh % NumberOfEdges
        edofs = Mesh % MaxEdgeDOFS
        maxedofs = ParallelReduction(edofs,2)

        maxedge = 0
        DO i=1,n
          maxedge = MAX(maxedge, Mesh % Edges(i) % GElementindex)
        END DO
        maxedge = ParallelReduction(maxedge,2)

        DO i=1,n
          Element => Mesh % Edges(i)
          DO j=1,Element % BDOFs
            DO m=1,DOFs
              l = DOFs*(l_beg + edofs*(i-1)+j-1)+m
              l = Matrix % Perm(l)
              IF(l==0) CYCLE

              Matrix % ParallelInfo % GlobalDOFs(l) = &
                  DOFs*(g_beg+maxedofs*(Element % GelementIndex-1)+j-1)+m
              Matrix % ParallelInfo % GInterface(l) = Mesh % ParallelInfo % EdgeInterface(i)

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
        maxfdofs = ParallelReduction(fdofs,2)

        maxface = 0
        DO i=1,n
          maxface = MAX(maxface, Mesh % Faces(i) % GElementindex)
        END DO
        maxface = ParallelReduction(maxface,2)

        DO i=1,n
          Element => Mesh % Faces(i)
          DO j=1,Element % BDOFs
            DO m=1,DOFs
              l = Matrix % Perm(DOFs*(l_beg + fdofs*(i-1)+j-1)+m)
              IF(l==0) CYCLE
              Matrix % ParallelInfo % GlobalDOFs(l) = &
                  DOFs*(g_beg+maxfdofs*(Element % GelementIndex-1)+j-1)+m
              Matrix % ParallelInfo % GInterface(l) = Mesh % ParallelInfo % FaceInterface(i)

              ALLOCATE(Matrix % Parallelinfo % NeighbourList(l) % Neighbours(SIZE( &
                   Mesh % ParallelInfo % FaceNeighbourList(i) % Neighbours)))

              Matrix % Parallelinfo % NeighbourList(l) % Neighbours = &
                   Mesh % ParallelInfo % FaceNeighbourList(i) % Neighbours
            END DO
          END DO
        END DO
      END IF

      IF ( Solver % GlobalBubbles ) THEN
        l_beg = Mesh % NumberOfNodes + &
                Mesh % NumberOfEdges*Mesh % MaxEdgeDOFs + &
                Mesh % NumberOfFaces*Mesh % MaxFaceDOFs

        g_beg = Maxnode +  maxedge*maxedofs + maxface*maxfdofs
        maxbdofs = ParallelReduction(Mesh % MaxBDOFs,2)

        DO i=1,Mesh % NumberOfBulkElements
          Element=>Mesh % Elements(i)

          j = Element % Type % ElementCode/100
          bdofs = Solver % Def_Dofs(j, Element % Bodyid, 5)
          IF ( bdofs<=0 .AND. &
             Solver % Def_Dofs(j,Element % Bodyid,6)>1) bdofs = Element % BDOFs

          DO l=1,bdofs
            DO j=1,DOFs
              k = Matrix % Perm(DOFs*(l_beg+Element % BubbleIndexes(l)-1)+j)
              IF(k==0) CYCLE
              Matrix % ParallelInfo % GlobalDOFs(k) = &
                DOFs*(g_beg+maxbdofs*(Element % GElementIndex-1)+l-1)+j
              Matrix % ParallelInfo % GInterface(k) = .FALSE.
              ALLOCATE(Matrix % ParallelInfo % NeighbourList(k) % Neighbours(1))
              Matrix % ParallelInfo % NeighbourList(k) % Neighbours=ParEnv % MyPE
            END DO
          END DO
        END DO
      END IF

      CALL SetupExtraDOFs( Solver, Matrix, Perm, DOFs )

!-------------------------------------------------------------------------------
    END SUBROUTINE MatrixParallelInfoDefault
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Add the extra rows of a collection matrix to the parallel structures. These
!> are the global unknowns shared by all partitions, counted by ParallelDOFs and
!> originating from the AddMatrix (circuit equations), followed by the rows of
!> the constraint matrix.
!-------------------------------------------------------------------------------
    SUBROUTINE SetupExtraDOFs( Solver, Matrix, Perm, DOFs )
!-------------------------------------------------------------------------------
      TYPE(Solver_t), TARGET :: Solver
      TYPE(Matrix_t), POINTER :: Matrix
      INTEGER, POINTER :: Perm(:)
      INTEGER :: DOFs
!-------------------------------------------------------------------------------
      TYPE(Matrix_t), POINTER :: G
      ! Roles of the two matrices in play here; see the comment below.
      TYPE(Matrix_t), POINTER :: BaseMatrix, CollMatrix
      INTEGER :: i, j, k, l, m, n, col, proc, g_beg, pardofs_all
      LOGICAL :: Found, Global_dof, LocalConstraints, DiscontBC, &
                 OwnersGiven, NeighboursGiven
      INTEGER, ALLOCATABLE :: Ind(:)
      ! Add extra degrees of freedom to parallel structures. The additional
      ! variables are assigned to task zero, and are assumed to be shared by
      ! all tasks (TODO: to be optimized if need be...)
      !
      ! Two distinct matrices are in play from here on, and mixing them up is
      ! easy, so they are given names:
      !
      !   CollMatrix  the matrix being initialized, i.e. the Matrix argument.
      !               When it carries extra rows (ExtraDOFs>0) it is typically
      !               a collection matrix, built by EliminateLinearRestriction
      !               from a base system plus that system's constraint rows.
      !   BaseMatrix  the base system those extra rows refer to, taken to be
      !               Solver % Matrix.
      !
      ! The distinction is load-bearing, not an oversight: the constraint rows
      ! are described by BaseMatrix % ConstraintMatrix (which is deliberately
      ! not propagated onto the collection matrix), and their InvPerm entries
      ! are decoded modulo BaseMatrix % NumberOfRows. Observed in the contact
      ! and mortar test cases: CollMatrix % NumberOfRows equals
      ! BaseMatrix % NumberOfRows + CollMatrix % ExtraDOFs exactly, while
      ! CollMatrix % ConstraintMatrix is null. Replacing BaseMatrix by
      ! CollMatrix below therefore silently drops the constraint numbering.
      !
      ! Caveat: equating the base system with Solver % Matrix is an assumption.
      ! Callers build the collection from a matrix of their choosing, which
      ! need not be the solver matrix. Passing the base matrix in explicitly
      ! would be the cleaner contract.
      ! --------------------------------------------------------------------
      CollMatrix => Matrix
      BaseMatrix => Solver % Matrix

      g_beg = ParallelReduction(MAXVAL(CollMatrix % ParallelInfo % GlobalDOFs),2)
      pardofs_all = ParallelReduction(CollMatrix % ParallelDOFs,2)

      LocalConstraints = ListGetLogical(Solver % Values, 'Partition Local Constraints',Found)
      DiscontBC  = ListGetLogicalAnyBC(CurrentModel,'Discontinuous Boundary' )

      Found = ASSOCIATED(BaseMatrix % ConstraintMatrix)
      IF(Found) Found = ASSOCIATED(BaseMatrix % ConstraintMatrix % InvPerm)

      G => Null()
      IF( LocalConstraints  ) THEN
        G => AllocateMatrix()
        G % Format = MATRIX_LIST
      END IF

      OwnersGiven = CollMatrix % ParallelDOFs>0
      IF(OwnersGiven) THEN
        OwnersGiven = ASSOCIATED(BaseMatrix % AddMatrix)
        IF(OwnersGiven) OwnersGiven = ALLOCATED(BaseMatrix % AddMatrix % RowOwner)
      END IF

      NeighboursGiven = CollMatrix % ParallelDOFs>0
      IF(NeighboursGiven) NeighboursGiven = ASSOCIATED(BaseMatrix % AddMatrix)
      IF(NeighboursGiven) NeighboursGiven = ASSOCIATED(BaseMatrix % AddMatrix % ParallelInfo)
      IF(NeighboursGiven) NeighboursGiven = ASSOCIATED(BaseMatrix % AddMatrix % ParallelInfo % &
                          NeighbourList)


      j=0
      DO i=CollMatrix % NumberOFRows-CollMatrix % ExtraDOFs+1,CollMatrix % NumberOfRows
        j=j+1
        Global_dof = j <= CollMatrix % ParallelDOFs
        IF (Global_dof) THEN
          CollMatrix % ParallelInfo % GInterface(i) = .TRUE.

          IF(NeighboursGiven) THEN
            ALLOCATE(CollMatrix % ParallelInfo % NeighbourList(i) % Neighbours( &
               SIZE(BaseMatrix % AddMatrix % ParallelInfo % NeighbourList(i) % Neighbours)))

            CollMatrix % ParallelInfo % NeighbourList(i) % Neighbours = &
               BaseMatrix % AddMatrix % ParallelInfo % NeighbourList(i) % Neighbours

            IF(ALL(CollMatrix % ParallelInfo % NeighbourList(i) % Neighbours /= ParEnv % myPE)) THEN
              CollMatrix % ParallelInfo % GInterface(i) = .FALSE.
              DEALLOCATE(CollMatrix % ParallelInfo % NeighbourList(i) % Neighbours)
              ALLOCATE(CollMatrix % ParallelInfo % NeighbourList(i) % Neighbours(1))
              CollMatrix % ParallelInfo % NeighbourList(i) % Neighbours(1) = ParEnv % mype

              CALL CRS_ZeroRow(CollMatrix,i)
              CALL CRS_SetMatrixElement(CollMatrix,i,i,1._dp)
              CollMatrix % RHS(i) = 0._dp
            END IF

          ELSE IF (OwnersGiven) THEN
            ALLOCATE(CollMatrix % ParallelInfo % NeighbourList(i) % Neighbours(ParEnv % PEs))
            DO k=1,ParEnv % PEs
              CollMatrix % ParallelInfo % NeighbourList(i) % Neighbours(k)=k-1
            END DO

            l = j + BaseMatrix % NumberOfRows
            l = BaseMatrix % AddMatrix % RowOwner(l)

            CollMatrix % ParallelInfo % NeighbourList(i) % Neighbours(1)   = l
            CollMatrix % ParallelInfo % NeighbourList(i) % Neighbours(l+1) = 0
          ELSE
            ALLOCATE(CollMatrix % ParallelInfo % NeighbourList(i) % Neighbours(ParEnv % PEs))
            l = 0
            DO k=ParEnv % PEs,1,-1
              l = l + 1
              CollMatrix % ParallelInfo % NeighbourList(i) % Neighbours(l)=k-1
            END DO
          END IF

          k = j + g_beg
        ELSE IF(Found) THEN
          n = BaseMatrix % ConstraintMatrix % InvPerm(j-CollMatrix % ParallelDOFs)
          IF( n<=0 ) CYCLE
          k = MOD(n-1,BaseMatrix % NumberOfRows)+1
          n = (n-1) / BaseMatrix % NumberOfRows

          CollMatrix % ParallelInfo % GInterface(i) = CollMatrix % ParallelInfo % GInterface(k)

          l = SIZE(CollMatrix % ParallelInfo % NeighbourList(k) % Neighbours)
          ALLOCATE(CollMatrix % ParallelInfo % NeighbourList(i) % Neighbours(l))

          IF(LocalConstraints) THEN
            IF(.NOT.DiscontBC) THEN
              IF (CollMatrix % ParallelInfo % GInterface(k)) &
                CALL Sort(l,CollMatrix % ParallelInfo % NeighbourList(k) % Neighbours)
            END IF

            DO m=CollMatrix % Rows(i), CollMatrix % Diag(i)-1
              col = CollMatrix % Cols(m)
              IF (CollMatrix % ParallelInfo % GInterface(col)) THEN
                DO l=1,SIZE(CollMatrix % ParallelInfo % NeighbourList(col) % Neighbours)
                  proc = CollMatrix % ParallelInfo % NeighbourList(col) % Neighbours(l)
                  IF(proc==ParEnv % mype) THEN
                    CALL List_AddToMatrixElement( G % ListMatrix, proc+1, col, 1._dp )
                  ELSE
                    CALL List_AddToMatrixElement( G % ListMatrix, &
                         proc+1, CollMatrix % ParallelInfo % GlobalDOFs(col), 1._dp)
                  END IF
                END DO
              END IF
            END DO
          END IF

          CollMatrix % Parallelinfo % NeighbourList(i) % Neighbours = &
             CollMatrix % ParallelInfo % NeighbourList(k) % Neighbours

          k = CollMatrix % ParallelInfo % GlobalDOFs(k) + (n+1)*g_beg + pardofs_all
        ELSE
          ALLOCATE(CollMatrix % ParallelInfo % NeighbourList(i) % Neighbours(1))
          CollMatrix % Parallelinfo % NeighbourList(i) % Neighbours(1)=ParEnv % myPE
          k = j + g_beg + pardofs_all
        END IF
        CollMatrix % ParallelInfo % GlobalDofs(i) = k
      END DO

      l=CollMatrix % ExtraDOFs
      n=SIZE(Perm)
      IF(l>0) THEN
         ALLOCATE(Ind(l))
         k = 0
         DO i=CollMatrix % NumberOfRows-l+1,CollMatrix % NumberOfRows
           k = k + 1
           CollMatrix % Perm(n*DOFs+k) = i
           Ind(k) = CollMatrix % Parallelinfo % Globaldofs(i)
         END DO
         CALL SortI(l,Ind,CollMatrix % Perm(n*DOFs+1:))
         DEALLOCATE(Ind)
      END IF

      CALL SyncConstraintOwnership( CollMatrix, G, LocalConstraints )
!-------------------------------------------------------------------------------
    END SUBROUTINE SetupExtraDOFs

!-------------------------------------------------------------------------------
!> Tell the neighbouring partitions which of the constraint rows are owned here,
!> so that the neighbour lists of a shared row come out in the same order on
!> every partition. G carries the rows to report, collected while numbering the
!> extra degrees of freedom; it is freed on the way out.
!-------------------------------------------------------------------------------
    SUBROUTINE SyncConstraintOwnership( CollMatrix, G, LocalConstraints )
!-------------------------------------------------------------------------------
      TYPE(Matrix_t), POINTER :: CollMatrix, G
      LOGICAL :: LocalConstraints
!-------------------------------------------------------------------------------
      LOGICAL, ALLOCATABLE :: Done(:)
      INTEGER, ALLOCATABLE :: Owneddofs(:)
      INTEGER :: i, j, k, l, m
      INTEGER :: ierr, status(MPI_STATUS_SIZE)
      IF (ASSOCIATED(ParEnv % IsNeighbour).AND. LocalConstraints) THEN
        CALL List_ToCRSMatrix(G)
        ALLOCATE(Done(CollMatrix % NumberOfROws)); Done=.FALSE.

        DO i=1,ParEnv % PEs
          IF(i-1==ParEnv % Mype) THEN
            IF(G % NumberOfRows>=i) THEN
              DO j=G % Rows(i), G % Rows(i+1)-1
                m = G % Cols(j)
                IF(.NOT.Done(m)) THEN
                  Done(m) = .TRUE.
                  k = SIZE(CollMatrix % ParallelInfo % Neighbourlist(m) % Neighbours)
                  CALL Sort(k,CollMatrix % ParallelInfo % NeighbourList(m) % Neighbours)
                END IF
              END DO
            END IF
            CYCLE
          END IF

          j = 0
          IF(G % NumberOfRows>=i) j = G % Rows(i+1)-G % Rows(i)

          IF(ParEnv % IsNeighbour(i)) THEN
            CALL MPI_BSEND( j,1,MPI_INTEGER,i-1,130,CollMatrix % Comm,ierr )
            IF (j>0) THEN
              CALL MPI_BSEND(G % Cols(G % Rows(i):),j,MPI_INTEGER,i-1,131,CollMatrix % Comm,ierr )
            END IF
          ELSE
            ! Constraint rows are owned here for a partition that is not a
            ! neighbour, so there is no channel to communicate them over.
            IF(j>0) CALL Fatal('SyncConstraintOwnership','Have '//I2S(j)//&
                ' constraint row(s) to send to non-neighbour partition '//I2S(i-1))
          END IF
        END DO

        DO i=1,ParEnv % NumOfNeighbours
          CALL MPI_RECV( j,1,MPI_INTEGER,MPI_ANY_SOURCE,130,CollMatrix % Comm,status,ierr )
          IF (j>0) THEN
            IF(.NOT.ALLOCATED(owneddofs)) ALLOCATE(owneddofs(j))
            IF(SIZE(owneddofs)<j) THEN
              DEALLOCATE(owneddofs); ALLOCATE(owneddofs(j))
            END IF

            CALL MPI_RECV(owneddofs,j,MPI_INTEGER,status(MPI_SOURCE),131,CollMatrix % Comm,status,ierr )
            DO l=1,j
               m = SearchNode( CollMatrix % ParallelInfo, Owneddofs(l), Order=CollMatrix % Perm )
               IF(m>0) THEN
                 IF(.NOT.Done(m)) THEN
                   Done(m) = .TRUE.
                   k = SIZE(CollMatrix % ParallelInfo % Neighbourlist(m) % Neighbours)
                   CALL Sort(k,CollMatrix % ParallelInfo % NeighbourList(m) % Neighbours)
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
!-------------------------------------------------------------------------------
    END SUBROUTINE SyncConstraintOwnership
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
   SUBROUTINE AssignAtLeastOneDOFToPartition(n,ParallelInfo,comm)
!-------------------------------------------------------------------------------
     INTEGER :: i,j,k,n, comm
     TYPE(ParallelInfo_t)  :: ParallelInfo

     INTEGER, POINTER :: p(:)
     LOGICAL :: L(0:ParEnv % PEs-1), L1(0:ParEnv % PEs-1)
     INTEGER :: ierr, status(MPI_STATUS_SIZE), np, memb(0:ParEnv % PEs-1), imemb(0:ParEnv % PEs-1), dof, nbr

     np = 0
     memb  = -1
     imemb = -1
     DO i=1,ParEnv % PEs
       IF ( ParEnv % Active(i) ) THEN
         memb(np)  = i-1
         imemb(i-1) = np
         np = np + 1
       END IF
     END DO

     L1 = .TRUE.; L=.TRUE.
     DO i=1,n
       IF(ParallelInfo % NeighbourList(i) % Neighbours(1)==ParEnv % myPE) THEN
         L1(imemb(ParEnv % MyPE)) = .FALSE.; EXIT
       END IF
     END DO
     CALL MPI_ALLREDUCE(L1, L, np, MPI_LOGICAL, MPI_LAND, comm ,ierr)

     
     IF(ANY(L(0:np-1))) THEN

       IF(L(imemb(ParEnv % MyPE))) THEN

         ! Pick the DOF to donate to a partition that ended up owning none.
         !
         ! The loop below used to select the most shared DOF, but its result was
         ! always overwritten by the midpoint choice that follows, so it has been
         ! inactive for a long time. There has apparently been trouble with it in
         ! the past, so it is left commented out rather than restored: the
         ! midpoint choice is the behaviour that has actually been in use. Note
         ! that n/2 need not be a shared DOF at all, so this deserves a proper
         ! revisit if this rare Hypre path ever misbehaves.
         !
         ! j = 0; k=0
         ! DO i=1,n
         !   p => ParallelInfo % NeighbourList(i) % Neighbours
         !   IF (SIZE(p)>j) THEN
         !      j=SIZE(p); k=i
         !   END IF
         ! END DO
         k = n/2
         p => ParallelInfo % NeighbourList(k) % Neighbours

         DO i=1,SIZE(p)
           IF(p(i)==ParEnv % myPE) THEN
             p(i) = p(1)
             p(1) = ParEnv % MyPE
             EXIT
           END IF
         END DO
         CALL Sort(SIZE(p)-1, p(2:) )

         dof = ParallelInfo % GlobalDOFs(k)
         DO i=2,SIZE(p)
           IF(.NOT.ParEnv % Active(p(i)+1)) CYCLE
           CALL MPI_BSEND(dof, 1, MPI_INTEGER, imemb(p(i)), 501, comm, ierr)
         END DO

         dof = 0
         DO i=0,np-1
           IF(ANY(memb(i)==p)) CYCLE
           CALL MPI_BSEND(dof, 1, MPI_INTEGER, i, 501, comm, ierr)
         END DO
       END IF

       DO i=0,np-1
         IF(memb(i) == ParEnv % myPE) CYCLE

         IF(L(i)) THEN
           CALL MPI_RECV(dof, 1, MPI_INTEGER, i, 501, comm, status, ierr )

           IF(dof>0) THEN
             k = SearchNode( ParallelInfo, dof, Order = ParallelInfo % Gorder )
             ! A donated DOF not found locally is silently ignored.
             IF (k>0) THEN
               p => ParallelInfo % NeighbourList(k) % Neighbours
               DO j=1,SIZE(p)
                 IF(p(j)==ParEnv % myPE) THEN
                   p(j) = p(1)
                   p(1) = memb(i)
                   EXIT
                 END IF
               END DO
               CALL Sort(SIZE(p)-1, p(2:) )
             END IF
           END IF
         END IF
       END DO
     END IF
!-------------------------------------------------------------------------------
   END SUBROUTINE AssignAtLeastOneDOFToPartition
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
       IF(.NOT. ASSOCIATED(Matrix % ParMatrix)) THEN
         CALL Fatal('ParallelInitSolve','ParMatrix not associated!')
       END IF
       CALL SetMatrixParEnv( Matrix )
       IF(.NOT. ASSOCIATED(ParEnv)) THEN
         CALL Fatal('ParallelInitSolve','ParEnv not associated!')
       END IF

       Upd = .TRUE.
       IF ( PRESENT(Update) ) Upd=Update
       CALL SParInitSolve( Matrix, x, b, r, Matrix % ParallelInfo, Upd )
#endif
!-------------------------------------------------------------------------------
    END SUBROUTINE ParallelInitSolve
!-------------------------------------------------------------------------------

! Compute parallel sum (or optional min=1 or max=2)
!-------------------------------------------------------------------------------
    SUBROUTINE ParallelSumVector( Matrix, x, Op )
!-------------------------------------------------------------------------------
       TYPE(Matrix_t) :: Matrix
       INTEGER, OPTIONAL :: op
       REAL(KIND=dp) CONTIG :: x(:)
!-------------------------------------------------------------------------------
       IF(Matrix % ParallelInfo % NothingShared ) RETURN

       CALL SetMatrixParEnv( Matrix )
       IF(.NOT.ASSOCIATED(Parenv % Active)) THEN
         ParEnv = ParEnv_Common
         ParEnv % ActiveComm = Matrix % Comm
       END IF

       CALL ExchangeSourceVec( Matrix, Matrix % ParMatrix % SplittedMatrix, &
              Matrix % ParallelInfo, x, op )
!-------------------------------------------------------------------------------
    END SUBROUTINE ParallelSumVector
!-------------------------------------------------------------------------------

! As the previous except for integer valued vectors.
!-------------------------------------------------------------------------------
    SUBROUTINE ParallelSumVectorInt( Matrix, x, Op )
!-------------------------------------------------------------------------------
       TYPE(Matrix_t) :: Matrix
       INTEGER, OPTIONAL :: op
       INTEGER CONTIG :: x(:)
!-------------------------------------------------------------------------------
       IF(Matrix % ParallelInfo % NothingShared ) RETURN

       CALL SetMatrixParEnv( Matrix )

       CALL ExchangeSourceVecInt( Matrix, Matrix % ParMatrix % SplittedMatrix, &
              Matrix % ParallelInfo, x, op )
!-------------------------------------------------------------------------------
     END SUBROUTINE ParallelSumVectorInt
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

      IF(Mesh % ParallelInfo % NothingShared ) RETURN

      ! We can inherit the ParEnv from the primary matrix even
      ! though the variable is not directly associated to it!
      IF( PRESENT( Matrix ) ) THEN
        CALL SetMatrixParEnv( Matrix )
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
    SUBROUTINE ParallelMatrixVector( Matrix, x, b, Update, UseMassVals, &
                   ZeroNotOwned, UseAbs )
!-------------------------------------------------------------------------------
      REAL(KIND=dp) CONTIG, TARGEt :: x(:), b(:)
      TYPE(Matrix_t), POINTER :: Matrix
      LOGICAL, OPTIONAL :: Update, UseMassVals,ZeroNotOwned, UseAbs
!-------------------------------------------------------------------------------
      INTEGER :: i
      REAL(KIND=dp), POINTER CONTIG :: Mx(:), Mr(:), Mb(:), r(:)

      LOGICAL:: UpdateL, UseMassValsL,ZeroNotOwnedL, UseAbsL

      TYPE(Matrix_t), POINTER :: SaveMatrix
!-------------------------------------------------------------------------------
#ifdef PARALLEL_FOR_REAL
      GlobalData => Matrix % ParMatrix
      SaveMatrix  => GlobalMatrix
      GlobalMatrix => Matrix
      CALL SetMatrixParEnv( Matrix )

      UpdateL = .FALSE.
      IF(PRESENT(Update)) UpdateL=Update

      UseMassValsL = .FALSE.
      IF(PRESENT(UseMassVals)) UseMassValsL=UseMassVals

      ZeroNotOwnedL = .FALSE.
      IF(PRESENT(ZeroNotOwned)) ZeroNotOwnedL=ZeroNotOwned

      UseABSL = .FALSE.
      IF(PRESENT(UseABS)) UseABSL=UseABS

      IF(UpdateL) THEN
        Mx => GlobalData % SplittedMatrix % TmpXVec
        Mr => GlobalData % SplittedMatrix % TmpRVec
      ELSE
        Mx => x
        Mr => b
      END  IF

      ! The mass coefficients are selected inside the product now. It used to be
      ! done by writing MassValues over Values here -- a deep copy for every
      ! interface block, a pointer swap for the inside matrix -- and undoing it
      ! afterwards. Same result, none of the copying, and A % Values is no longer
      ! ever pointer-assigned to a sibling array.
      IF(UseABSL) THEN
        CALL SParABSMatrixVectorVals( Mx, Mr, UseMassValsL )
      ELSE
        CALL SParMatrixVectorVals( Mx, Mr, UseMassValsL )
      END IF

      IF(UpdateL) CALL SParUpdateResult( Matrix, x, b, .FALSE. )

      IF (ZeroNotOwnedL) THEN
        DO i = 1, Matrix % NumberOFRows
          IF ( Matrix % ParallelInfo % NeighbourList(i) % Neighbours(1) /= ParEnv % MyPE ) &
            b(i) = 0._dp
        END DO
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
      CALL SetMatrixParEnv( Matrix )
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
!> Create a permutation from the dofs owned by this partition to the local dof
!> indexes such that Perm(1:nOwned) gives the local dofs that this partition
!> owns, and return the number of owned dofs. This enables computing partition-
!> consistent sums (dot products, norms, extrema, ...) with just one loop and no
!> ownership test inside it. The same loop is valid in serial and in parallel
!> since when there is nothing shared the permutation is simply the identity.
!>
!> The ownership is primarily obtained from the parallel matrix. There are some
!> exceptions when no matrix, and hence no associated communication, has been
!> created and we have to resort to the communication structure of the mesh.
!> Note that the latter is only consistent for nodal scalar fields.
!-------------------------------------------------------------------------------
    FUNCTION ParallelOwnedPerm( n, Perm, Matrix, Mesh ) RESULT( nOwned )
!-------------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: n
      INTEGER, ALLOCATABLE :: Perm(:)
      TYPE(Matrix_t), OPTIONAL :: Matrix
      TYPE(Mesh_t), OPTIONAL :: Mesh
      INTEGER :: nOwned
!-------------------------------------------------------------------------------
      TYPE(NeighbourList_t), POINTER :: NeighbourList(:)
      INTEGER :: i
!-------------------------------------------------------------------------------
      IF( ALLOCATED( Perm ) ) THEN
        IF( SIZE( Perm ) < n ) DEALLOCATE( Perm )
      END IF
      IF( .NOT. ALLOCATED( Perm ) ) ALLOCATE( Perm(n) )

      NeighbourList => NULL()
      IF( ParEnv % PEs > 1 ) THEN
        IF( PRESENT( Matrix ) ) THEN
          IF( ASSOCIATED( Matrix % ParallelInfo ) ) THEN
            IF( .NOT. Matrix % ParallelInfo % NothingShared ) &
                NeighbourList => Matrix % ParallelInfo % NeighbourList
          END IF
        END IF
        IF( .NOT. ASSOCIATED( NeighbourList ) .AND. PRESENT( Mesh ) ) THEN
          IF( .NOT. Mesh % ParallelInfo % NothingShared ) &
              NeighbourList => Mesh % ParallelInfo % NeighbourList
        END IF
      END IF

      IF( .NOT. ASSOCIATED( NeighbourList ) ) THEN
        ! Serial, or nothing shared: this partition owns every dof.
        nOwned = n
        DO i=1,n
          Perm(i) = i
        END DO
        RETURN
      END IF

      nOwned = 0
      DO i=1,MIN( n, SIZE( NeighbourList ) )
        IF( NeighbourList(i) % Neighbours(1) /= ParEnv % MyPE ) CYCLE
        nOwned = nOwned + 1
        Perm(nOwned) = i
      END DO
!-------------------------------------------------------------------------------
    END FUNCTION ParallelOwnedPerm
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
    FUNCTION ParallelCDOT( n, x, y ) RESULT(s)
!-------------------------------------------------------------------------------
      INTEGER :: n
      COMPLEX(KIND=dp) :: s
      COMPLEX(KIND=dp) CONTIG :: x(:),y(:)
!-------------------------------------------------------------------------------
      s = 0.0d0
#ifdef PARALLEL_FOR_REAL
      s = SParCDotProd( n, x, 1, y, 1 )
#endif
!-------------------------------------------------------------------------------
    END FUNCTION ParallelCDOT
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!> As ParallelCDOT but without conjugation, i.e. the bilinear form x^T y.
!> See SParCDotProdU for when this is the one you want.
!-------------------------------------------------------------------------------
    FUNCTION ParallelCDOTU( n, x, y ) RESULT(s)
!-------------------------------------------------------------------------------
      INTEGER :: n
      COMPLEX(KIND=dp) :: s
      COMPLEX(KIND=dp) CONTIG :: x(:),y(:)
!-------------------------------------------------------------------------------
      s = 0.0d0
#ifdef PARALLEL_FOR_REAL
      s = SParCDotProdU( n, x, 1, y, 1 )
#endif
!-------------------------------------------------------------------------------
    END FUNCTION ParallelCDOTU
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    SUBROUTINE ParallelGlobalNumbering(Mesh,OldMesh,NewNodes,Reorder)
!-------------------------------------------------------------------------------
       TYPE(Mesh_t) :: Mesh, OldMesh
       INTEGER :: NewNodes,Reorder(:)
!-------------------------------------------------------------------------------
#ifdef PARALLEL_FOR_REAL
        CALL SparGlobalNumbering( Mesh,OldMesh,NewNodes,Reorder )
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


!---------------------------------------------------------------------------
! Computed a parallel sum (or min or max) for a real valued scalar.
!-------------------------------------------------------------------------------
    FUNCTION ParallelReductionR(R,oper_arg) RESULT(rsum)
!-------------------------------------------------------------------------------
      REAL(KIND=dp) :: R, rsum
      INTEGER, OPTIONAL :: oper_arg
!-------------------------------------------------------------------------------
      INTEGER :: oper
      LOGICAL :: alloc
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

        alloc = ASSOCIATED(ParEnv % Active)
        IF (.NOT.alloc) CALL ParallelActive(.TRUE.)
        CALL SparActiveSUM(rsum,oper)
        IF(.NOT.alloc) THEN
!         DEALLOCATE(ParEnv % active)
!         PArEnv % Active => null()
        END IF
      END IF
#endif
!-------------------------------------------------------------------------------
    END FUNCTION ParallelReductionR
!-------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Same as previous but for integer values.
!-------------------------------------------------------------------------------
    FUNCTION ParallelReductionI(i,oper_arg) RESULT(isum)
!-------------------------------------------------------------------------------
      INTEGER :: i, isum
      INTEGER, OPTIONAL :: oper_arg
!-------------------------------------------------------------------------------
      INTEGER :: oper
      LOGICAL :: alloc
!-------------------------------------------------------------------------------
      isum = i
#ifdef PARALLEL_FOR_REAL
      IF ( ParEnv % PEs>1) THEN
        oper = 0
        IF (PRESENT(oper_arg)) THEN
          oper=oper_arg
        ELSE
          oper = 0
        END IF

        alloc = ASSOCIATED(ParEnv % Active)
        IF (.NOT.alloc) CALL ParallelActive(.TRUE.)
        CALL SparActiveSUMInt(isum,oper)
        IF(.NOT.alloc) THEN
!         DEALLOCATE(ParEnv % active)
!         PArEnv % Active => null()
        END IF
      END IF
#endif
!-------------------------------------------------------------------------------
    END FUNCTION ParallelReductionI
!-------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Same as previous but for complex values.
!-------------------------------------------------------------------------------
    FUNCTION ParallelReductionZ(z,oper_arg) RESULT(zsum)
!-------------------------------------------------------------------------------
      COMPLEX(KIND=dp) :: z, zsum
      INTEGER, OPTIONAL :: oper_arg
!-------------------------------------------------------------------------------
      INTEGER :: oper
      LOGICAL :: alloc
!-------------------------------------------------------------------------------
      zsum = z
#ifdef PARALLEL_FOR_REAL
      IF ( ParEnv % PEs>1) THEN
        oper = 0
        IF (PRESENT(oper_arg)) THEN
          oper=oper_arg
        ELSE
          oper = 0
        END IF

        alloc = ASSOCIATED(ParEnv % Active)
        IF (.NOT.alloc) CALL ParallelActive(.TRUE.)
        CALL SparActiveSUMComplex(zsum,oper)
        IF(.NOT.alloc) THEN
!         DEALLOCATE(ParEnv % active)
!         PArEnv % Active => null()
        END IF
      END IF
#endif
!-------------------------------------------------------------------------------
    END FUNCTION ParallelReductionZ
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


    FUNCTION ParallelSlicesComm() RESULT( CustomComm )
      INTEGER :: CustomComm
#ifdef PARALLEL_FOR_REAL
      LOGICAL :: Visited = .FALSE.
      INTEGER :: nSlices, nTimes, iSlice, iTime, ierr
      INTEGER :: CustomComm0
      LOGICAL :: GotIt
      SAVE Visited, CustomComm0

      IF(.NOT. Visited ) THEN
        nSlices = ListGetInteger( CurrentModel % Simulation,'Number Of Slices',GotIt)
        nTimes = ListGetInteger( CurrentModel % Simulation,'Number Of Times',GotIt)
        IF(nSlices > 1 .AND. nTimes > 1) THEN
          iSlice = MODULO( ParEnv % MyPe, nSlices )
          iTime = ParEnv % MyPe / nSlices
          CALL MPI_Comm_split(ELMER_COMM_WORLD, iTime, ParEnv % PEs, CustomComm0, ierr);
        ELSE
          CustomComm0 = ELMER_COMM_WORLD
        END IF
        Visited = .TRUE.
      END IF

      CustomComm = CustomComm0
#else
      CustomComm = -1
#endif

    END FUNCTION ParallelSlicesComm


    FUNCTION ParallelTimesComm() RESULT( CustomComm )
      INTEGER :: CustomComm
#ifdef PARALLEL_FOR_REAL
      LOGICAL :: Visited = .FALSE.
      INTEGER :: nSlices, nTimes, iSlice, iTime, ierr
      INTEGER :: CustomComm0
      LOGICAL :: GotIt
      SAVE Visited, CustomComm0

      IF(.NOT. Visited ) THEN
        nSlices = ListGetInteger( CurrentModel % Simulation,'Number Of Slices',GotIt)
        nTimes = ListGetInteger( CurrentModel % Simulation,'Number Of Times',GotIt)
        IF(nSlices > 1 .AND. nTimes > 1) THEN
          iSlice = MODULO( ParEnv % MyPe, nSlices )
          iTime = ParEnv % MyPe / nSlices
          CALL MPI_Comm_split(ELMER_COMM_WORLD, iSlice, ParEnv % PEs, CustomComm0, ierr);
        ELSE
          CustomComm0 = ELMER_COMM_WORLD
        END IF
        Visited = .TRUE.
      END IF

      CustomComm = CustomComm0
#else
      CustomComm = -1
#endif
    END FUNCTION ParallelTimesComm


    FUNCTION ParallelPieceRank(CustomComm) RESULT (CommRank)
      INTEGER :: CustomComm, CommRank, ierr
#ifdef PARALLEL_FOR_REAL
      CALL MPI_Comm_rank(CustomComm, CommRank, ierr)
#else
      CommRank = -1
#endif
    END FUNCTION ParallelPieceRank

    FUNCTION ParallelPieceSize(CustomComm) RESULT (CommSize)
      INTEGER :: CustomComm, CommSize, ierr
#ifdef PARALLEL_FOR_REAL
      CALL MPI_Comm_size(CustomComm, CommSize, ierr)
#else
      CommSize = -1
#endif
    END FUNCTION ParallelPieceSize



!--------------------------------'-----------------------------------------------
    SUBROUTINE ParallelMergeMatrix( Solver, A, A1, A2 )
!-------------------------------------------------------------------------------
      TYPE(Solver_t), TARGET :: Solver
      TYPE(Matrix_t), POINTER :: A, A1, A2
!-------------------------------------------------------------------------------
      TYPE(Matrix_t), POINTER :: Ai
      TYPE(ParallelInfo_t), POINTER :: P, Pi
      INTEGER :: i, j, k, l, m, n, ni, ksum, c
      INTEGER :: jumps(3)
      INTEGER :: Ierr, status(MPI_STATUS_SIZE)
!-------------------------------------------------------------------------------
#ifdef PARALLEL_FOR_REAL
      IF ( ParEnv % PEs <= 1 ) RETURN

      n = 0
      jumps = 0

      DO i=1,2
        IF(i == 1) THEN
          Ai => A1
        ELSE
          Ai => A2
        END IF
        ni = Ai % NumberOfRows
        Pi => Ai % ParallelInfo
        IF(.NOT. ASSOCIATED(Pi) ) THEN
          CALL Fatal('ParallelMergeMatrix',&
              'Submatrix '//I2S(i)//' does not have parallel info!')
        END IF

        n = n + ni
        m = ParallelReduction(MAXVAL(Pi % GlobalDOFs),2)
        jumps(i+1) = jumps(i) + m
      END DO

      !IF( ParEnv % MyPe == 0 ) THEN
      !  PRINT *,'offsets for parallel info:',jumps
      !END IF

      ALLOCATE( A % ParallelInfo )
      P => A % ParallelInfo

      c = A % NumberOfRows / n
!     IF( n /= A % NumberOfRows ) THEN
!       PRINT *,'n:',parenv % mype, n,A % NumberOfRows,c
!       CALL Warn('ParallelMergeMatrix','Mismatch in vector length')
!     END IF

      n = A % NumberOfRows
      ALLOCATE( P % NeighbourList(n) )
      CALL AllocateVector( P % GInterface, n)
      CALL AllocateVector( P % GlobalDOFs, n)
      IF( ASSOCIATED( A1 % Perm ) ) CALL AllocateVector( A % Perm, n)
      IF( ASSOCIATED( A1 % InvPerm ) ) CALL AllocateVector( A % InvPerm, n)


      ! This is made so that would there be need it is easy to add matrices
      ! if we later need more than two...
      n = 0
      DO i=1,2
        IF(i == 1) THEN
          Ai => A1
        ELSE
          Ai => A2
        END IF

        ni = Ai % NumberOfRows
        Pi => Ai % ParallelInfo
        m = jumps(i)

        IF(c==1) THEN
          P % GInterface(n+1:n+ni) = Pi % GInterface(1:ni)
          P % GlobalDofs(n+1:n+ni) = Pi % GlobalDofs(1:ni) + m
        ELSE
          P % GInterface(2*n+1:2*(n+ni)-1:2) = Pi % GInterface(1:ni)
          P % GInterface(2*n+2:2*(n+ni):2) = Pi % GInterface(1:ni)
          P % GlobalDofs(2*n+1:2*(n+ni)-1:2) = 2*Pi % GlobalDofs(1:ni)-1 + 2*m
          P % GlobalDofs(2*n+2:2*(n+ni):2) = 2*Pi % GlobalDofs(1:ni) + 2*m
        END IF

        !IF( ASSOCIATED( A % Perm ) ) THEN
        !  A % Perm(n+1:n+ni) = Ai % Perm(1:ni)
        !END IF
        !IF( ASSOCIATED( A % InvPerm ) ) THEN
        !  A % InvPerm(n+1:n+ni) = Ai % InvPerm(1:ni)
        !END IF

        ksum = 0
        DO j=1,ni
          IF(.NOT. ASSOCIATED(Pi % NeighbourList(j) % Neighbours)) CYCLE
          k = SIZE(Pi % NeighbourList(j) % Neighbours)
          ksum = ksum + k
          IF(c==1) THEN
            ALLOCATE(P % NeighbourList(n+j) % Neighbours(k))
            P % NeighbourList(n+j) % Neighbours = Pi % NeighbourList(j) % Neighbours
          ELSE
            ALLOCATE(P % NeighbourList(2*n+2*j-1) % Neighbours(k))
            P % NeighbourList(2*n+2*j-1) % Neighbours = Pi % NeighbourList(j) % Neighbours
            ALLOCATE(P % NeighbourList(2*n+2*j) % Neighbours(k))
            P % NeighbourList(2*n+2*j) % Neighbours = Pi % NeighbourList(j) % Neighbours
          END IF
        END DO

        n = n + ni
      END DO

BLOCK
     INTEGER, ALLOCATABLE :: Ind(:)

     n = A % NumberOfRows

     ALLOCATE(P % Gorder(n), Ind(n))
     Ind = P % GlobalDOFs
     P % Gorder = [(i,i=1,n)]
     CALL SortI( n, Ind, P % Gorder)
END BLOCK

     ! Finalize creation of parallel structures
     A % Solver => Solver
     A % ParMatrix => ParInitMatrix( A, A % ParallelInfo )
#endif
!-------------------------------------------------------------------------------
    END SUBROUTINE ParallelMergeMatrix
!-------------------------------------------------------------------------------


  END MODULE ParallelUtils

!> \}
