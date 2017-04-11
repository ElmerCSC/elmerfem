!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
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
! ******************************************************************************
! *
! *  Authors: Olivier Gagliardini, Juha Ruokolainen
! *  Email:   Juha.Ruokolainen [at] csc.fi 
! *  Web:     http://elmerice.elmerfem.org
! *  Address: CSC - Scientific Computing Ltd.  
! *               Keilaranta 14                    
! *               02101 Espoo, Finland             
! *                                                 
! *       Original Date: 14 May 2007                
! * 
! *****************************************************************************
!>  Module containing solver for computation of surface normals
SUBROUTINE ComputeNormalSolver( Model, Solver, dt, TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh, materials, BCs, etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear & nonlinear equation solver options
!
!  REAL(KIND=dp) :: dt,
!     INPUT: Timestep size for time dependent simulations
!
!  LOGICAL :: TransientSimulation
!     INPUT: Steady state or transient simulation
!
!Modification have been done to deal with problems where we don't want to compute
!the normal on all the boundaries (problem with the corners)
!
!Keywords are: ComputeAll = True/False computing / or not the normal on every boundary
!              ComputeNormal = True to compute the normal on a given boundary
!******************************************************************************
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation

!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Element_t), POINTER :: Element
  TYPE(Variable_t), POINTER :: NormalSolution
  TYPE(ValueList_t), POINTER :: BC, SolverParams

  INTEGER :: i, j, k, l, m, n, cnt, ierr, t, DIM
  REAL(KIND=dp) :: u, v, w, s 
  REAL(KIND=dp), POINTER :: Nvector(:), PassNVector(:), RecvNVector(:)
  INTEGER, POINTER :: Permutation(:), Neighbours(:)
  INTEGER, ALLOCATABLE :: NeighbourPerm(:), PassCount(:), RecvCount(:),&
       PassIndices(:,:), RecvIndices(:,:), LocalPerm(:)
  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
  TYPE(Nodes_t), SAVE :: Nodes
  REAL(KIND=dp) :: Bu, Bv, Normal(3), NormalCond(4)

  LOGICAL :: CompAll = .TRUE., CompBC = .TRUE., Found, Parallel, &
       FirstTime = .TRUE.,UnFoundFatal=.TRUE.
  LOGICAL, ALLOCATABLE :: Hit(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName = 'ComputeNormalSolver'

  SAVE :: NeighbourPerm, PassCount, RecvCount, PassIndices, &
       RecvIndices, LocalPerm, Hit

  Parallel = (ParEnv % PEs > 1)
  Mesh => Solver % Mesh
  DIM = CoordinateSystemDimension()

  !---------------------------------------
  ! Setup pointer to the current solution:
  !---------------------------------------
  NormalSolution => VariableGet( Solver % Mesh % Variables, 'Normal Vector',UnFoundFatal=UnFoundFatal)
  Nvector => NormalSolution % Values
  Permutation => NormalSolution % Perm
  NVector = 0.0_dp

  CALL INFO(SolverName, 'Computing Normal Vector for Nodes', level=3)

  SolverParams => GetSolverParams()
  CompAll = GetLogical (SolverParams, 'ComputeAll', Found)

  IF (.NOT.Found) THEN
    WRITE(Message,'(A)') ('ComputeAll not found, Normal is computed for all the boundaries')
    CALL INFO(SolverName, Message, level = 20)
    CompAll = .TRUE.
    CompBC = .TRUE.
  ELSE
    IF (.NOT.CompAll) CompBC = .FALSE.
    IF (CompAll) CompBC = .TRUE.
  END IF

  IF((FirstTime .OR. Solver % Mesh % Changed) .AND. Parallel) THEN
     
     IF(.NOT. FirstTime) THEN
        DEALLOCATE(Hit, NeighbourPerm, PassCount, RecvCount, &
             PassIndices, RecvIndices, LocalPerm)
     END IF

     ALLOCATE(Hit(Mesh % NumberOfNodes), &
          NeighbourPerm(ParEnv % PEs),&
          PassCount(ParEnv % NumOfNeighbours),&
          RecvCount(ParEnv % NumOfNeighbours),&
          PassIndices(ParEnv % NumOfNeighbours, Mesh % NumberOfNodes),&
          RecvIndices(ParEnv % NumOfNeighbours, Mesh % NumberOfNodes),&
          LocalPerm( MAXVAL(Mesh % ParallelInfo % GlobalDOFs)))

     NeighbourPerm = 0
     LocalPerm = 0
     Hit = .FALSE.
     RecvCount = 0
     PassCount = 0
     PassIndices = 0

     !Create a list of partitions with which we share nodes
     cnt = 0
     DO i=1,ParEnv % PEs
        IF(ParEnv % MyPE == i-1) CYCLE
        IF(.NOT. ParEnv % IsNeighbour(i)) CYCLE
        cnt = cnt + 1
        NeighbourPerm(i) = cnt
     END DO

     !Perm to get back from global to local node numbering
     DO i=1, Mesh % NumberOfNodes
        LocalPerm(Mesh % ParallelInfo % GlobalDOFs(i)) = i
     END DO

  END IF !FirstTime

  DO t = 1, Solver % Mesh % NumberOfBoundaryElements
    Element => GetBoundaryElement(t)

    IF(Element % PartIndex /= ParEnv % MyPE) CYCLE

    n = GetElementNOFNodes(Element)
    BC => GetBC(Element)
    IF (n == 1) CYCLE

    CALL GetElementNodes( Nodes, Element )

    IF (.NOT.CompAll) THEN
      CompBC = GetLogical ( BC,'ComputeNormal',Found)
      IF (.NOT.Found) THEN
          NormalCond = 0.0
          NormalCond(1:n) = GetReal( BC, 'ComputeNormal Condition', Found )

	! If at least one value in NormalCond > 0, then CompBC=.true.
        IF (COUNT(NormalCond > 0.0) > 0) CompBC = .TRUE.
      END IF
    END IF

    IF (CompBC) THEN
       DO i = 1,n
          IF (NormalCond(i) .LT. 0) CYCLE
          j = Element % NodeIndexes( i )
          k = Permutation(j)

          Bu = Element % Type % NodeU(i)
          IF ( Element % Type % Dimension > 1 ) THEN
             Bv = Element % Type % NodeV(i)
          ELSE
             Bv = 0.0D0
          END IF
          Normal = NormalVector(Element, Nodes, Bu, Bv, .TRUE.)
          Nvector(DIM*(k-1)+1:DIM*k) = Nvector(DIM*(k-1)+1:DIM*k) +& 
               Normal(1:DIM)

          IF(Parallel) Hit(j) = .TRUE.
       END DO
    END IF
 END DO

  IF(Parallel) THEN

     IF(FirstTime .OR. Solver % Mesh % Changed) THEN

        !Find nodes on partition boundaries
        FirstTime = .FALSE.
        RecvCount = 0
        PassCount = 0
        PassIndices = 0
        DO i=1,Model % Mesh % NumberOfNodes
           IF(.NOT.Hit(i)) CYCLE !no normal computation
           IF(.NOT. Mesh % ParallelInfo % INTERFACE(i)) CYCLE !no partition interface
           Neighbours => Mesh % ParallelInfo % NeighbourList(i) % Neighbours

           DO j=1,SIZE(Neighbours)
              IF(Neighbours(j) == ParEnv % MyPE) CYCLE

              m = NeighbourPerm(Neighbours(j)+1)
              PassCount(m) = PassCount(m) + 1
              PassIndices(m, PassCount(m)) = i
           END DO
        END DO
     END IF

     !Send
     DO i=1,ParEnv % PEs
        IF(NeighbourPerm(i) == 0) CYCLE
        m = NeighbourPerm(i)

        !Send count
        CALL MPI_BSEND(PassCount(m), 1, MPI_INTEGER, i-1, 200, ELMER_COMM_WORLD, ierr)
        !Send (global) node numbers
        CALL MPI_BSEND(Mesh % ParallelInfo % GlobalDOFs(PassIndices(m,1:PassCount(m))),&
             PassCount(m), MPI_INTEGER, i-1, 201, ELMER_COMM_WORLD, ierr)

        !Construct normal vector array to pass
        ALLOCATE(PassNVector(PassCount(m) * DIM))
        DO j=1, PassCount(m)
           l = Permutation(PassIndices(m,j))
           PassNVector(DIM*(j-1)+1:DIM*j) = NVector(DIM*(l-1)+1:DIM*l)
        END DO

        !Send normal vectors
        CALL MPI_BSEND(PassNVector, PassCount(m)*DIM, MPI_DOUBLE_PRECISION, &
             i-1, 202, ELMER_COMM_WORLD, ierr)

        DEALLOCATE(PassNVector)
     END DO

     !Receive
     DO i=1,ParEnv % PEs

        IF(NeighbourPerm(i) == 0) CYCLE
        m = NeighbourPerm(i)

        !Receive count
        CALL MPI_RECV(RecvCount(m), 1, MPI_INTEGER, i-1, 200, ELMER_COMM_WORLD, status, ierr)

        !Receive (global) node numbers
        CALL MPI_RECV(RecvIndices(m,1:RecvCount(m)), RecvCount(m), MPI_INTEGER, &
             i-1, 201, ELMER_COMM_WORLD, status, ierr)

        !Receive normal vectors
        ALLOCATE(RecvNVector(RecvCount(m)*DIM))
        CALL MPI_RECV(RecvNVector, RecvCount(m)*DIM, MPI_DOUBLE_PRECISION, &
             i-1, 202, ELMER_COMM_WORLD, status, ierr)
        
        DO j=1, RecvCount(m)
           IF(.NOT. Hit(LocalPerm(RecvIndices(m,j)))) CYCLE !passed a halo node
           k = Permutation(LocalPerm(RecvIndices(m,j)))
           NVector(DIM*(k-1)+1:DIM*k) = NVector(DIM*(k-1)+1:DIM*k) + RecvNVector(DIM*(j-1)+1:DIM*j)
        END DO

        DEALLOCATE(RecvNVector)
     END DO

  END IF !Parallel
  
  DO i=1,Model % NumberOfNodes
    k = Permutation(i)
    IF ( k > 0 ) THEN
      s = SQRT( SUM( Nvector(DIM*(k-1)+1:DIM*k)**2 ) )

      IF ( s /= 0.0D0 ) THEN
	Nvector(DIM*(k-1)+1:DIM*k) = Nvector(DIM*(k-1)+1:DIM*k)/s 
      END IF
    END IF

  END DO

  CALL INFO(SolverName, 'End', level=3)

!------------------------------------------------------------------------------
END SUBROUTINE ComputeNormalSolver
!------------------------------------------------------------------------------
