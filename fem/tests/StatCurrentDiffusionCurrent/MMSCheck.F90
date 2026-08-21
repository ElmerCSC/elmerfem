!------------------------------------------------------------------------------
!> Evidence gate for the three-ion diffusion-current manufactured solution.
!>
!> Every quantity is derived from the model at run time -- Faraday's constant
!> from Constants, the conductivity from the material, and the species contract
!> from the potential solver's own keywords -- so this checker cannot drift out
!> of step with case.sif.  The exact potential follows from the species fields
!> themselves rather than from a hard-coded amplitude:
!>
!>   sigma*phi(x) = -F * ( SUM_i z_i D_i c_i(x) - SUM_i z_i D_i c_i(0) )
!>
!> which is the solution of  div(sigma grad phi) = -div(F SUM_i z_i D_i grad c_i)
!> on a strip whose potential is pinned to zero at both ends.  That the ends
!> agree is a PRECONDITION of the case, so it is asserted rather than assumed.
!>
!> Two modes, DERIVED from whether the solver that owns Potential declared
!> `Diffusion Current Species` -- never from a keyword of this checker's own,
!> which could disagree with the case and then check the wrong thing quietly:
!>   opted in  - the manufactured field must be reproduced AND be present.
!>   opted out - the ohmic solver must return exactly zero AND therefore MISS
!>               the manufactured field.  That control lives in
!>               fem/tests/StatCurrentDiffusionCurrentOff.
!------------------------------------------------------------------------------
SUBROUTINE CheckDiffusionCurrentMMS(Model, Solver, dt, TransientSimulation)
  USE DefUtils
  IMPLICIT NONE

  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation

  TYPE(Variable_t), POINTER :: Potential, Current, Species
  TYPE(ValueList_t), POINTER :: Params, PotentialParams, Material
  REAL(KIND=dp), ALLOCATABLE :: WeightedSum(:)
  REAL(KIND=dp) :: Faraday, Conductivity, Charge, Diffusivity
  REAL(KIND=dp) :: Reference, ExactPotential, ExactPeak, EndMismatch
  REAL(KIND=dp) :: PotentialError, CurrentMaximum, PotentialPeak
  REAL(KIND=dp) :: PotentialTolerance, CurrentTolerance, ZeroTolerance
  REAL(KIND=dp) :: XMinimum, XMaximum
  INTEGER :: Node, Permutation, Component, Dim, SpeciesIndex, SpeciesCount
  INTEGER :: LowNode, HighNode
  LOGICAL :: Found, ExpectDiffusion
  CHARACTER(LEN=MAX_NAME_LEN) :: Keyword, SpeciesName

  Params => GetSolverParams()

  Potential => VariableGet(Model % Mesh % Variables, 'Potential')
  Current => VariableGet(Model % Mesh % Variables, 'Volume Current')
  IF (.NOT. ASSOCIATED(Potential)) &
      CALL Fatal('CheckDiffusionCurrentMMS', 'Potential variable not found')
  IF (.NOT. ASSOCIATED(Current)) &
      CALL Fatal('CheckDiffusionCurrentMMS', 'Volume Current variable not found')
  IF (.NOT. ASSOCIATED(Potential % Solver)) &
      CALL Fatal('CheckDiffusionCurrentMMS', 'Potential has no owning solver')

  ! The mode is DERIVED, never declared: whichever solver owns Potential is
  ! asked whether it opted in.  A separate keyword here could disagree with the
  ! case it is checking, and would then be checking the wrong thing quietly.
  PotentialParams => Potential % Solver % Values
  SpeciesCount = ListGetInteger(PotentialParams, 'Diffusion Current Species', &
      Found, minv=1)
  IF (.NOT. Found) SpeciesCount = 0
  ExpectDiffusion = (SpeciesCount > 0)

  Dim = CoordinateSystemDimension()
  Material => Model % Materials(1) % Values
  Faraday = ListGetConstReal(Model % Constants, 'Faraday Constant', Found)
  IF (.NOT. Found) CALL Fatal('CheckDiffusionCurrentMMS', &
      'Faraday Constant is missing from the Constants section')
  Conductivity = ListGetConstReal(Material, 'Electric Conductivity', Found)
  IF (.NOT. Found) CALL Fatal('CheckDiffusionCurrentMMS', &
      'Electric Conductivity is missing from Material 1')

!------------------------------------------------------------------------------
! Build SUM_i z_i D_i c_i node by node, from the fields the solver actually saw.
!------------------------------------------------------------------------------
  ALLOCATE(WeightedSum(Model % NumberOfNodes))
  WeightedSum = 0.0_dp

  DO SpeciesIndex = 1, MAX(SpeciesCount, 1)
    WRITE(Keyword, '(A,I0)') 'Diffusion Current Variable ', SpeciesIndex
    SpeciesName = ListGetString(PotentialParams, TRIM(Keyword), Found)
    IF (.NOT. Found) THEN
      ! The feature-off control declares no species on the potential solver; it
      ! still needs the exact field, so the same registered keyword is read from
      ! the checker's own section instead.
      SpeciesName = ListGetString(Params, TRIM(Keyword), Found)
    END IF
    IF (.NOT. Found) THEN
      IF (SpeciesIndex == 1) CALL Fatal('CheckDiffusionCurrentMMS', &
          'No species named by either the potential solver or the checker')
      EXIT
    END IF

    Species => VariableGet(Model % Mesh % Variables, TRIM(SpeciesName))
    IF (.NOT. ASSOCIATED(Species)) CALL Fatal('CheckDiffusionCurrentMMS', &
        'Species variable not found: '//TRIM(SpeciesName))

    Charge = ListGetConstReal(Material, TRIM(SpeciesName)//' Ion Charge', Found)
    IF (.NOT. Found) CALL Fatal('CheckDiffusionCurrentMMS', &
        'Missing ion charge for '//TRIM(SpeciesName))
    Diffusivity = ListGetConstReal(Material, TRIM(SpeciesName)//' Diffusivity', Found)
    IF (.NOT. Found) CALL Fatal('CheckDiffusionCurrentMMS', &
        'Missing diffusivity for '//TRIM(SpeciesName))

    DO Node = 1, Model % NumberOfNodes
      Permutation = Species % Perm(Node)
      IF (Permutation <= 0) CALL Fatal('CheckDiffusionCurrentMMS', &
          'Species undefined at a node: '//TRIM(SpeciesName))
      WeightedSum(Node) = WeightedSum(Node) + &
          Charge * Diffusivity * Species % Values(Permutation)
    END DO
  END DO

!------------------------------------------------------------------------------
! Precondition: the strip's two potential-pinned ends must carry the same
! SUM z D c, or the closed form below does not hold and nothing can be decided.
!------------------------------------------------------------------------------
  XMinimum = MINVAL(Model % Nodes % x(1:Model % NumberOfNodes))
  XMaximum = MAXVAL(Model % Nodes % x(1:Model % NumberOfNodes))
  LowNode = MINLOC(ABS(Model % Nodes % x(1:Model % NumberOfNodes) - XMinimum), 1)
  HighNode = MINLOC(ABS(Model % Nodes % x(1:Model % NumberOfNodes) - XMaximum), 1)
  Reference = WeightedSum(LowNode)
  EndMismatch = ABS(WeightedSum(HighNode) - Reference)

  ExactPeak = 0.0_dp
  DO Node = 1, Model % NumberOfNodes
    ExactPeak = MAX(ExactPeak, &
        ABS(Faraday * (WeightedSum(Node) - Reference) / Conductivity))
  END DO

  IF (ExactPeak <= 0.0_dp) CALL Fatal('CheckDiffusionCurrentMMS', &
      'The manufactured field is identically zero -- this case tests nothing')
  IF (EndMismatch * Faraday / Conductivity > 1.0e-3_dp * ExactPeak) CALL Fatal( &
      'CheckDiffusionCurrentMMS', &
      'Species differ between the pinned ends; the closed form does not apply')

  PotentialTolerance = 1.0e-6_dp * ExactPeak
  ZeroTolerance = 1.0e-9_dp * ExactPeak
  CurrentTolerance = 1.0e-6_dp * Faraday * ExactPeak * Conductivity

!------------------------------------------------------------------------------
! Measure.
!------------------------------------------------------------------------------
  PotentialError = 0.0_dp
  PotentialPeak = 0.0_dp
  CurrentMaximum = 0.0_dp

  DO Node = 1, Model % NumberOfNodes
    ExactPotential = -Faraday * (WeightedSum(Node) - Reference) / Conductivity

    Permutation = Potential % Perm(Node)
    IF (Permutation <= 0) &
        CALL Fatal('CheckDiffusionCurrentMMS', 'Potential undefined at a node')
    PotentialError = MAX(PotentialError, &
        ABS(Potential % Values(Permutation) - ExactPotential))
    PotentialPeak = MAX(PotentialPeak, ABS(Potential % Values(Permutation)))

    Permutation = Current % Perm(Node)
    IF (Permutation <= 0) &
        CALL Fatal('CheckDiffusionCurrentMMS', 'Volume Current undefined at a node')
    DO Component = 1, Dim
      CurrentMaximum = MAX(CurrentMaximum, &
          ABS(Current % Values(Current % DOFs * (Permutation - 1) + Component)))
    END DO
  END DO

  WRITE(*,'(A,ES16.8)') 'MMS analytic peak potential   = ', ExactPeak
  WRITE(*,'(A,ES16.8)') 'MMS max potential error       = ', PotentialError
  WRITE(*,'(A,ES16.8)') 'MMS max computed potential    = ', PotentialPeak
  WRITE(*,'(A,ES16.8)') 'MMS max total current         = ', CurrentMaximum

!------------------------------------------------------------------------------
! Gate.  Each branch must be failable by the OTHER branch's solution: the ohmic
! field trips "signal present", and the diffusion field trips "stays zero".
!------------------------------------------------------------------------------
  IF (ExpectDiffusion) THEN
    IF (PotentialError > PotentialTolerance) CALL Fatal( &
        'CheckDiffusionCurrentMMS', 'Manufactured potential is outside tolerance')
    IF (PotentialPeak < 0.99_dp * ExactPeak) CALL Fatal( &
        'CheckDiffusionCurrentMMS', &
        'Potential carries no manufactured signal; an ohmic solve would look like this')
    IF (CurrentMaximum > CurrentTolerance) CALL Fatal( &
        'CheckDiffusionCurrentMMS', 'Total current does not cancel')
    WRITE(*,'(A)') 'MMS verdict                   = PASS (diffusion current present)'
  ELSE
    IF (PotentialPeak > ZeroTolerance) CALL Fatal( &
        'CheckDiffusionCurrentMMS', &
        'Feature-off control is not ohmic; the diffusion current leaked into an opt-out case')
    IF (PotentialError < 0.5_dp * ExactPeak) CALL Fatal( &
        'CheckDiffusionCurrentMMS', &
        'Feature-off control did not miss the manufactured field; the case tests nothing')
    WRITE(*,'(A)') 'MMS verdict                   = PASS (ohmic control, feature off)'
  END IF

  DEALLOCATE(WeightedSum)
END SUBROUTINE CheckDiffusionCurrentMMS
