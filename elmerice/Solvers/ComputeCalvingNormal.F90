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

!NB: This is just a local copy of ComputeNormalSolver, with modified keywords
! to allow proper normal computation on two adjacent boundaries.
SUBROUTINE ComputeCalvingNormalSolver( Model, Solver, dt, TransientSimulation )
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
  TYPE(ValueList_t), POINTER :: SolverParams

  INTEGER :: i, j, k, l, m, n, cnt, ierr, t, DIM
  REAL(KIND=dp) :: u, v, w, s 
  REAL(KIND=dp), POINTER :: Nvector(:), PassNVector(:), RecvNVector(:)
  INTEGER, POINTER :: Permutation(:), Neighbours(:)
  INTEGER, ALLOCATABLE :: NeighbourPerm(:), PassCount(:), RecvCount(:),&
       PassIndices(:,:), RecvIndices(:,:), LocalPerm(:)
  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
  TYPE(Nodes_t), SAVE :: Nodes
  REAL(KIND=dp) :: Bu, Bv, Normal(3)

  LOGICAL :: Found, Parallel, MeActive=.TRUE.
  LOGICAL, ALLOCATABLE :: Hit(:), PartActive(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName = 'ComputeCalvingNormalSolver'

  Parallel = (ParEnv % PEs > 1)
  Mesh => Solver % Mesh
  DIM = CoordinateSystemDimension()
 
  !---------------------------------------
  ! Setup pointer to the current solution:
  !---------------------------------------
  NormalSolution => Solver % Variable
  IF ( ASSOCIATED( NormalSolution ) ) THEN 
     Nvector => NormalSolution % Values
     Permutation => NormalSolution % Perm
     Nvector = 0.0_dp !wipe out previous
  ELSE
     PRINT *,'FATAL: Unable to set pointer to the current solution'
     STOP
  END IF


  CALL INFO(SolverName, 'Computing Normal Vector for Nodes', level=3)

  SolverParams => GetSolverParams()

  !Check if current partition has anything to do
  IF(Parallel) THEN
    MeActive = .TRUE.
    IF(Solver % NumberOfActiveElements <= 0) MeActive = .FALSE.
    ALLOCATE(PartActive(ParEnv % PEs))

    CALL MPI_ALLGATHER(MeActive,1,MPI_LOGICAL,&
         PartActive,1,MPI_LOGICAL, ELMER_COMM_WORLD, ierr)

    IF(.NOT. MeActive) RETURN

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
      IF(.NOT. PartActive(i)) CYCLE
      IF(.NOT. ParEnv % IsNeighbour(i)) CYCLE
      cnt = cnt + 1
      NeighbourPerm(i) = cnt
    END DO

    !Perm to get back from global to local node numbering
    DO i=1, Mesh % NumberOfNodes
      LocalPerm(Mesh % ParallelInfo % GlobalDOFs(i)) = i
    END DO

  END IF !Parallel

  DO t = 1, Solver % NumberOfActiveElements
    Element => GetActiveElement(t)

    IF(Element % PartIndex /= ParEnv % MyPE) CYCLE

    n = GetElementNOFNodes(Element)
    IF (n == 1) CYCLE

    CALL GetElementNodes( Nodes, Element )

    DO i = 1,n
      j = Element % NodeIndexes( i )
      k = Permutation(j)

      IF(k <= 0) CALL Fatal(SolverName, &
           "Encountered invalid perm, this is a programming error")

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
  END DO

  IF(Parallel) THEN

    !Find nodes on partition boundaries
    RecvCount = 0
    PassCount = 0
    PassIndices = 0
    DO i=1,Model % Mesh % NumberOfNodes
      IF(.NOT.Hit(i)) CYCLE !no normal computation
      IF(.NOT. Mesh % ParallelInfo % INTERFACE(i)) CYCLE !no partition interface
      Neighbours => Mesh % ParallelInfo % NeighbourList(i) % Neighbours

      DO j=1,SIZE(Neighbours)
        IF(Neighbours(j) == ParEnv % MyPE) CYCLE
        IF(.NOT. PartActive(Neighbours(j)+1)) CYCLE

        m = NeighbourPerm(Neighbours(j)+1)
        PassCount(m) = PassCount(m) + 1
        PassIndices(m, PassCount(m)) = i
      END DO
    END DO

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
END SUBROUTINE ComputeCalvingNormalSolver
!------------------------------------------------------------------------------
