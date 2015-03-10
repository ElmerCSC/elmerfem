SUBROUTINE DifferenceSolver(Model, Solver, dt, TransientSimulation)
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Computes the difference between two scalar fields.
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
!******************************************************************************
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
  CHARACTER(LEN=MAX_NAME_LEN) :: SenderName, F1name, F2name, MeshName
  TYPE(Variable_t), POINTER :: F1, F2, Difference
  LOGICAL :: Found, AllocationsDone
  INTEGER :: n, t
  TYPE(Element_t), POINTER :: Element
  REAL(KIND=dp), ALLOCATABLE :: F1loc(:), F2loc(:)
  REAL(KIND=dp) :: L2nrm, L2nrmloc, L2nrmf1, L2nrmf1loc, L2nrmf2, L2nrmf2loc
  REAL(KIND=dp) :: H1nrm, H1nrmloc, H1nrmf1, H1nrmf1loc, H1nrmf2, H1nrmf2loc
  SAVE F1loc, F2loc, AllocationsDone
!------------------------------------------------------------------------------
  SenderName = 'DifferenceSolve'

  ! Print mesh name:
  !------------------
  MeshName = Solver % Mesh % Name

  WRITE(Message, *) 'Functions are interpolated using mesh: ', TRIM(MeshName)

  CALL Info(SenderName, Message)

  ! Allocate local storage:
  !-------------------------
  IF(.NOT. AllocationsDone .OR. Solver % Mesh % Changed) THEN
     n = Solver % Mesh % MaxElementDOFs

     IF(AllocationsDone) THEN
        DEALLOCATE(F1loc, F2loc)
     END IF

     ALLOCATE(F1loc(n), F2loc(n))

     AllocationsDone = .TRUE.
  END IF

  ! Fetch first scalar function:
  !------------------------------
  SolverParams => GetSolverParams()
  
  F1name = GetString(SolverParams, "F1", Found)
  
  IF(.NOT.Found) THEN
     WRITE(Message, *) 'F1 not found (name unspecified)'
     CALL Fatal(SenderName, Message)
  END IF

  F1 => VariableGet(Solver % Mesh % Variables, F1name)

  IF(.NOT.ASSOCIATED(F1)) THEN
     WRITE(Message, *) 'F1 not found: ', TRIM(F1name)
     CALL Fatal(SenderName, Message)
  ELSE
     WRITE(Message, *) 'F1: "', TRIM(F1name), '" interpolated on ', MeshName
     CALL Info(SenderName, Message)     
  END IF

  ! Fetch second scalar function:
  !-------------------------------
  F2name = GetString(SolverParams, "F2", Found)

  IF(.NOT.Found) THEN
     WRITE(Message, *) 'F2 not found (name unspecified)'
     CALL Fatal(SenderName, Message)
  END IF

  F2 => VariableGet(Solver % Mesh % Variables, F2name)

  IF(.NOT.ASSOCIATED(F2)) THEN
     WRITE(Message, *) 'F2 not found: ', TRIM(F2name)
     CALL Fatal(SenderName, Message)
  ELSE
     WRITE(Message, *) 'F2: "', TRIM(F2name), '" interpolated on ', MeshName
     CALL Info(SenderName, Message)     
  END IF

  ! Compute the difference:
  !-------------------------
  Difference => Solver % Variable
  Difference % Values = F1 % Values - F2 % Values

  ! Compute norms:
  !----------------
  L2nrm = 0.0d0
  L2nrmf1 = 0.0d0
  L2nrmf2 = 0.0d0

  H1nrm = 0.0d0
  H1nrmf1 = 0.0d0
  H1nrmf2 = 0.0d0

  DO t = 1, GetNofActive()
     Element => GetActiveElement(t)
     n = GetElementNOFNodes(Element)

     CALL GetScalarLocalSolution(F1loc, F1name)
     CALL GetScalarLocalSolution(F2loc, F2name)
     
     CALL Compute(Element, n, F1loc, F2loc, &
          L2nrmloc, L2nrmf1loc, L2nrmf2loc, &
          H1nrmloc, H1nrmf1loc, H1nrmf2loc)

     L2nrm = L2nrm + L2nrmloc
     L2nrmf1 = L2nrmf1 + L2nrmf1loc
     L2nrmf2 = L2nrmf2 + L2nrmf2loc

     H1nrm = H1nrm + H1nrmloc
     H1nrmf1 = H1nrmf1 + H1nrmf1loc
     H1nrmf2 = H1nrmf2 + H1nrmf2loc
  END DO

  L2nrm = SQRT(L2nrm)
  L2nrmf1 = SQRT(L2nrmf1)
  L2nrmf2 = SQRT(L2nrmf2)

  H1nrm = SQRT(H1nrm)
  H1nrmf1 = SQRT(H1nrmf1)
  H1nrmf2 = SQRT(H1nrmf2)

  WRITE(Message, *) '|| F1 - F2 || =', L2nrm
  CALL Info(SenderName, Message)
  WRITE(Message, *) '|| F1 || =', L2nrmf1
  CALL Info(SenderName, Message)
  WRITE(Message, *) '|| F2 || =', L2nrmf2
  CALL Info(SenderName, Message)

  WRITE(Message, *) '|| grad(F1 - F2) || =', H1nrm
  CALL Info(SenderName, Message)
  WRITE(Message, *) '|| grad F1 || =', H1nrmf1
  CALL Info(SenderName, Message)
  WRITE(Message, *) '|| grad F2 || =', H1nrmf2
  CALL Info(SenderName, Message)

CONTAINS

  SUBROUTINE Compute(Element, n, F1, F2, &
       L2nrm, L2nrmf1, L2nrmf2, &
       H1nrm, H1nrmf1, H1nrmf2)

    TYPE(Element_t), POINTER :: Element
    INTEGER :: n
    REAL(KIND=dp) :: F1(:), F2(:)
    REAL(KIND=dp) :: L2nrm, L2nrmf1, L2nrmf2
    REAL(KIND=dp) :: H1nrm, H1nrmf1, H1nrmf2

    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: detJ, Basis(n), dBasisdx(n, 3)
    REAL(KIND=dp) :: F1atIP, F2atIP
    REAL(KIND=dp) :: gradF1atIP(n), gradF2atIP(n)
    INTEGER :: t, i
    LOGICAL :: stat
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes

    CALL GetElementNodes(Nodes)

    IP = GaussPoints(Element)

    L2nrm = 0.0d0
    L2nrmf1 = 0.0d0
    L2nrmf2 = 0.0d0

    H1nrm = 0.0d0
    H1nrmf1 = 0.0d0
    H1nrmf2 = 0.0d0

    DO t = 1, IP % n
       stat = ElementInfo(Element, Nodes, &
            IP % U(t), IP % V(t), IP % W(t), &
            detJ, Basis, dBasisdx)

       F1atIP = SUM( F1(1:n) * Basis(1:n) )
       F2atIP = SUM( F2(1:n) * Basis(1:n) )

       gradF1atIP(1:3) = MATMUL(TRANSPOSE(dBasisdx(1:n,1:3)), F1(1:n) )
       gradF2atIP(1:3) = MATMUL(TRANSPOSE(dBasisdx(1:n,1:3)), F2(1:n) )

       L2nrm = L2nrm + (F1atIP - F2atIP)**2 * IP % s(t) * detJ
       L2nrmf1 = L2nrmf1 + F1atIP**2 * IP % s(t) * detJ
       L2nrmf2 = L2nrmf2 + F2atIP**2 * IP % s(t) * detJ

       DO i = 1, 3
          H1nrm = H1nrm + (gradF1atIP(i) - gradF2atIP(i))**2 * IP % s(t) * detJ
          H1nrmf1 = H1nrmf1 + gradF1atIP(i)**2 * IP % s(t) * detJ
          H1nrmf2 = H1nrmf2 + gradF2atIP(i)**2 * IP % s(t) * detJ
       END DO

    END DO

  END SUBROUTINE Compute

END SUBROUTINE DifferenceSolver
