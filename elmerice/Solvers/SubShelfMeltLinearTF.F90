!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
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
! *****************************************************************************/
! ******************************************************************************
! *
! *  Linear thermal forcing basal melt parameterisation:
! *
! *    m = gammaT * (rhow * cp / (rhoi * Lf)) * (Tw - Tf)
! *
! *  where Tf = lambda1 * S + lambda2 + cc * z  is the salinity- and
! *  pressure-dependent freezing temperature (ISOMIP+ linearisation),
! *  Tw is the ocean temperature read from variable temp_oce,
! *  S  is the ocean salinity read from variable sal_oce, and
! *  gammaT is a heat exchange velocity (m/s).
! *  The result is converted to m/yr.
! *
! *  Required Constants (in sif):
! *    Ice Density            (kg/m3)
! *    Water Density          (kg/m3)
! *    Latent Heat SI         (J/kg)
! *    SW Cp                  (J/kg/K)
! *
! *  Required Solver keywords (in sif):
! *    lower surface variable name  = string "Zb"
! *    grounded mask name           = string "GroundedMask"
! *    grounding line melt          = logical False
! *    gamma T                      = real 1.0e-4
! *
! *  Optional Solver keywords:
! *    water column scaling         = logical True
! *    water column scaling factor  = real 75.0
! *    bedrock variable name        = string "bedrock"   ! required if water column scaling = True
! *
! *  Required Elmer variables:
! *    temp_oce      : nodal ocean temperature (degC)
! *    sal_oce       : nodal ocean salinity (PSU)
! *    ice_thickness : ice thickness (m)
! *
! *  Output variables (created automatically):
! *    <var>_flux : integrated melt flux per element (m3/yr), element variable
! *
! *  Example sif solver section:
! *
! *  Solver 1
! *    Exec Solver = "before timestep"
! *    Equation = "Shelf melt"
! *    Procedure = "ElmerIceSolvers" "SubShelfMeltLinearTF"
! *    Variable = bmb
! *    Variable DOFs = 1
! *    lower surface variable name = string "Zb"
! *    grounded mask name          = string "GroundedMask"
! *    grounding line melt         = logical False
! *    gamma T                     = real 1.0e-4
! *    water column scaling        = logical True
! *    water column scaling factor = real 75.0
! *    bedrock variable name       = string "bedrock"
! *  End
! *
! *****************************************************************************/

SUBROUTINE SubShelfMeltLinearTF (Model, Solver, dt, Transient)

  USE DefUtils

  IMPLICIT NONE

  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t)          :: Model
  REAL(KIND=dp)          :: dt
  LOGICAL                :: Transient

  ! Elmer variables
  TYPE(Variable_t), POINTER :: z_iceBase, z_bedrock, groundedMask, T_oce_var, sal_oce_var, ice_thickness_var
  INTEGER, POINTER          :: T_oce_Perm(:), sal_oce_Perm(:)
  REAL(KIND=dp), POINTER    :: T_oce_vals(:), sal_oce_vals(:)

  ! Solver parameters
  CHARACTER(LEN=MAX_NAME_LEN) :: lowerSurfName, groundedMaskName, bedrockName
  REAL(KIND=dp) :: gammaT
  LOGICAL       :: glMelt, wct_sc

  ! Physical constants
  REAL(KIND=dp) :: rhoi, rhoo, Lf, SWCp, seconds_per_year
  ! Freezing point parameters (from ISOMIP+)
  REAL(KIND=dp), PARAMETER :: lambda1       = -0.0573_dp   ! freezing point salinity coeff (degC/PSU)
  REAL(KIND=dp), PARAMETER :: lambda2       =  0.0832_dp   ! freezing point offset (degC)
  REAL(KIND=dp), PARAMETER :: cc            = 7.61e-4_dp

  ! Local variables
  TYPE(ValueList_t), POINTER :: SolverParams, Material
  REAL(KIND=dp) :: T_freeze, T_far, S_far, meltRate, meltScaling, wct, wct_factor, min_ice_thickness
  LOGICAL       :: found, T_oce_found, sal_oce_found
  INTEGER       :: ii, ii_mat, nZeroedNodes, ierr

  ! For element-based integrated melt flux
  TYPE(Variable_t), POINTER         :: meltFlux_var
  REAL(KIND=dp), POINTER            :: FluxValues(:)
  INTEGER, POINTER                  :: FluxPerm(:)
  TYPE(Element_t), POINTER          :: Element
  TYPE(Nodes_t)                     :: ElementNodes
  TYPE(GaussIntegrationPoints_t)    :: IntegStuff
  REAL(KIND=dp) :: Basis(MAX_ELEMENT_NODES)
  REAL(KIND=dp) :: detJ, U, V, W, Sw, meltRate_gp, elemFlux
  LOGICAL       :: stat
  INTEGER       :: element_id, n_el, jj, kk

  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName = 'SubShelfMeltLinearTF'

  !----------------------------------------------------------------------------
  ! Read solver parameters
  !----------------------------------------------------------------------------
  SolverParams => GetSolverParams()

  gammaT = GetConstReal( SolverParams, 'gamma T', Found )
  IF (.NOT. Found) CALL FATAL(SolverName, 'No keyword >gamma T< found')

  lowerSurfName = GetString( SolverParams, 'lower surface variable name', Found )
  IF (.NOT. Found) CALL FATAL(SolverName, 'No keyword >lower surface variable name< found')

  groundedMaskName = GetString( SolverParams, 'grounded mask name', Found )
  IF (.NOT. Found) CALL FATAL(SolverName, 'No keyword >grounded mask name< found')

  glMelt = GetLogical( SolverParams, 'grounding line melt', Found )
  IF (.NOT. Found) CALL FATAL(SolverName, 'No keyword >grounding line melt< found')

  wct_sc = GetLogical( SolverParams, 'water column scaling', Found )
  IF (wct_sc) THEN
     wct_factor = GetConstReal( SolverParams, 'water column scaling factor', Found )
     IF (.NOT. Found) CALL FATAL(SolverName, 'No keyword >water column scaling factor< found')
     bedrockName = GetString( SolverParams, 'bedrock variable name', Found )
     IF (.NOT. Found) CALL FATAL(SolverName, 'No keyword >bedrock variable name< found')
  END IF

  ! Read Min H from Material section (minimum ice thickness for melt computation)
  Found = .FALSE.
  DO ii_mat = 1, CurrentModel % NumberOfMaterials
     Material => CurrentModel % Materials(ii_mat) % Values
     min_ice_thickness = GetConstReal(Material, 'Min H', Found)
     IF (Found) EXIT
  END DO
  IF (.NOT. Found) CALL FATAL(SolverName, 'Minimal ice thickness >Min H< not found in any Material section')

  !----------------------------------------------------------------------------
  ! Read physical constants
  !----------------------------------------------------------------------------
  rhoi = GetConstReal( CurrentModel % Constants, 'Ice Density', Found )
  IF (.NOT. Found) CALL FATAL(SolverName, &
       'Ice Density not found in Constants; should be defined in the Constants section of the sif file')
  rhoo = GetConstReal( CurrentModel % Constants, 'Water Density', Found )
  IF (.NOT. Found) CALL FATAL(SolverName, &
       'Water Density not found in Constants; should be defined in the Constants section of the sif file')
  Lf   = GetConstReal( CurrentModel % Constants, 'Latent Heat SI', Found )
  IF (.NOT. Found) CALL FATAL(SolverName, &
       'Latent Heat SI not found in Constants; should be defined in the Constants section of the sif file')
  SWCp = GetConstReal( CurrentModel % Constants, 'SW Cp', Found )
  IF (.NOT. Found) CALL FATAL(SolverName, &
       'SW Cp not found in Constants; should be defined in the Constants section of the sif file')
  seconds_per_year = GetConstReal( CurrentModel % Constants, 'seconds_per_year', Found )
  IF (.NOT. Found) CALL FATAL(SolverName, &
       'seconds_per_year not found in Constants; should be defined in the Constants section of the sif file')
  !----------------------------------------------------------------------------
  ! Get Elmer variables
  !----------------------------------------------------------------------------
  z_iceBase => VariableGet( Solver % Mesh % Variables, TRIM(lowerSurfName) )
  IF (.NOT. ASSOCIATED(z_iceBase)) CALL FATAL(SolverName, 'Failed to find ice base variable')

  groundedMask => VariableGet( Solver % Mesh % Variables, TRIM(groundedMaskName) )
  IF (.NOT. ASSOCIATED(groundedMask)) CALL FATAL(SolverName, 'Failed to find grounded mask variable')

  T_oce_var   => VariableGet( Solver % Mesh % Variables, 'temp_oce' )
  T_oce_found = ASSOCIATED(T_oce_var)
  IF (.NOT. ASSOCIATED(T_oce_var)) CALL FATAL(SolverName, &
      'Variable temp_oce not found. Please make sure it is a defined variable in your .sif')
  T_oce_vals => T_oce_var % Values
  T_oce_Perm => T_oce_var % Perm

  sal_oce_var   => VariableGet( Solver % Mesh % Variables, 'sal_oce' )
  sal_oce_found = ASSOCIATED(sal_oce_var)
  IF (.NOT. ASSOCIATED(sal_oce_var)) CALL FATAL(SolverName, &
      'Variable sal_oce not found. Please make sure it is a defined variable in your .sif')
  sal_oce_vals => sal_oce_var % Values
  sal_oce_Perm => sal_oce_var % Perm

  IF (wct_sc) THEN
    z_bedrock => VariableGet( Solver % Mesh % Variables, TRIM(bedrockName) )
    IF (.NOT. ASSOCIATED(z_bedrock)) THEN
      CALL FATAL(SolverName, &
        'Failed to find bedrock variable' &
      )
    END IF
  END IF

  ice_thickness_var => VariableGet(Solver % Mesh % Variables, 'H')
  IF (.NOT. ASSOCIATED(ice_thickness_var)) THEN
    CALL FATAL(SolverName, &
      'Failed to find ice thickness variable >H<. Please make sure it is ' // &
      'a defined variable in your .sif' &
    )
  END IF


  !----------------------------------------------------------------------------
  ! Automatically create element-based flux variable if absent
  !----------------------------------------------------------------------------
  meltFlux_var => VariableGet( Solver % Mesh % Variables, &
       TRIM(Solver % Variable % Name) // '_flux' )
  IF (.NOT. ASSOCIATED(meltFlux_var)) THEN
     ALLOCATE( FluxValues(Solver % NumberOfActiveElements) )
     ALLOCATE( FluxPerm(Solver % Mesh % NumberOfBulkElements + &
                        Solver % Mesh % NumberOfBoundaryElements) )
     FluxValues = 0.0_dp
     FluxPerm   = 0
     DO element_id = 1, Solver % NumberOfActiveElements
        Element => GetActiveElement(element_id)
        FluxPerm(Element % ElementIndex) = element_id
     END DO
     CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, Solver, &
          TRIM(Solver % Variable % Name) // '_flux', 1, FluxValues, FluxPerm)
     meltFlux_var => VariableGet( Solver % Mesh % Variables, &
          TRIM(Solver % Variable % Name) // '_flux' )
  END IF

  !----------------------------------------------------------------------------
  ! Initialise melt and flux to zero
  !----------------------------------------------------------------------------
  Solver % Variable % Values = 0.0_dp
  meltFlux_var % Values      = 0.0_dp

  !----------------------------------------------------------------------------
  ! Loop over nodes to compute nodal melt rate
  !----------------------------------------------------------------------------
  nZeroedNodes = 0
  DO ii = 1, Model % NumberOfNodes

     IF (Solver % Variable % Perm(ii) .LE. 0) CYCLE

     IF (groundedMask % Values(groundedMask % Perm(ii)) .GT. 0) CYCLE
     IF (groundedMask % Values(groundedMask % Perm(ii)) .EQ. 0) THEN
        IF (.NOT. glMelt) CYCLE
     END IF

     IF (ice_thickness_var % Values(ice_thickness_var % Perm(ii)) .LE. min_ice_thickness) THEN
         Solver % Variable % Values(Solver % Variable % Perm(ii)) = 0.0_dp
         nZeroedNodes = nZeroedNodes + 1
         CYCLE
     END IF

     S_far = sal_oce_vals(sal_oce_Perm(ii))

     T_freeze = lambda1 * S_far + lambda2 + cc * z_iceBase % Values(z_iceBase % Perm(ii))

     T_far = T_oce_vals(T_oce_Perm(ii))

     !TODO: At the moment, we do not allow negative melt rates (i.e. freezing) at the ice-ocean interface.
     !      Should probably be revisited in the future.
     meltRate = MAX(gammaT * (rhoo * SWCp / (rhoi * Lf)) * (T_far - T_freeze) *  seconds_per_year, 0.0_dp)
     IF (wct_sc) THEN
        wct         = z_iceBase % Values(z_iceBase % Perm(ii)) &
                    - z_bedrock % Values(z_bedrock % Perm(ii))
        meltScaling = TANH(wct / (wct_factor / EXP(1.0_dp)))
     ELSE
        meltScaling = 1.0_dp
     END IF

     Solver % Variable % Values(Solver % Variable % Perm(ii)) = meltRate * meltScaling

  END DO
  IF (ParEnv % PEs > 1) CALL MPI_Allreduce(MPI_IN_PLACE, nZeroedNodes, 1, &
       MPI_INTEGER, MPI_SUM, ELMER_COMM_WORLD, ierr)
  IF (ParEnv % MyPE == 0) THEN
   CALL INFO(SolverName, &
     'SubShelfMeltLinearTF: nodes set to zero due to >Min H<:' // &
     I2S(nZeroedNodes) &
     )
  END IF

  !----------------------------------------------------------------------------
  ! Loop over active elements to integrate nodal melt rate -> element flux
  !----------------------------------------------------------------------------
  DO element_id = 1, Solver % NumberOfActiveElements
     Element => GetActiveElement(element_id)
     n_el = GetElementNOFNodes(Element)
     CALL GetElementNodes(ElementNodes)
     IntegStuff = GaussPoints(Element)

     ! This uses the nodal melt rate to compute the integrated melt flux for this element via FEM interpolation
     elemFlux = 0.0_dp
     DO jj = 1, IntegStuff % n
        U  = IntegStuff % u(jj)
        V  = IntegStuff % v(jj)
        W  = IntegStuff % w(jj)
        Sw = IntegStuff % s(jj)
        stat = ElementInfo(Element, ElementNodes, U, V, W, detJ, Basis)

        ! Interpolate nodal melt rate to this Gauss point
        meltRate_gp = 0.0_dp
        DO kk = 1, n_el
           ii = Element % NodeIndexes(kk)
           IF (Solver % Variable % Perm(ii) > 0) THEN
              meltRate_gp = meltRate_gp + Basis(kk) * &
                   Solver % Variable % Values(Solver % Variable % Perm(ii))
           END IF
        END DO

        ! Convert ice-equivalent melt rate to water-equivalent melt rate
        ! (rhoi/rhoo) before integrating, so that elemFlux is water-equivalent
        elemFlux = elemFlux + Sw * detJ * meltRate_gp * (rhoi / rhoo)
     END DO

     ! checks if the element is active in the flux variable and assigns the integrated flux
     IF (meltFlux_var % Perm(Element % ElementIndex) > 0) THEN
        meltFlux_var % Values(meltFlux_var % Perm(Element % ElementIndex)) = elemFlux
     END IF

  END DO

  ! Compute and store the solution norm for reference norm checking
  Solver % Variable % Norm = ComputeNorm(Solver, SIZE(Solver % Variable % Values), &
     Solver % Variable % Values)

END SUBROUTINE SubShelfMeltLinearTF
