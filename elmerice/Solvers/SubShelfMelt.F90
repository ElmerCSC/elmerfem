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
! *  Authors: Rupert Gladstone, Ben Galton-Fenzi
! *  Email:   rupertgladstone1972@gmail.com
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Modified July 2023 to incorporate ISMIP6 BMB parameterisation
! *  (these comments are in need of updating...)
! *
! *  Calculate sub ice shelf melting as a function of depth.  A number of 
! *  parameters may determine the melt function and are specified in the sif.
! *  
! *  sub shelf melt function (string) is the name of one of the hard coded melt 
! *  functions in this routine.  Currently available are:
! *  bgf - melt as a function of pressure and thermal forcing created by Ben  
! *  Galton-Fenzi.
! *  bgf_scaled - same as bgf but scales down melt near grounding line based on 
! *  water column thickness
! *  The above 2 are as in this paper: https://tc.copernicus.org/articles/11/319/2017/
! *  
! *  ismip6 - see ismip6 wiki
! *  https://www.climate-cryosphere.org/wiki/index.php?title=ISMIP6-Projections-Antarctica#Oceanic_forcing:_temperature.2C_salinity.2C_thermal_forcing_and_melt_rate_parameterization
! *  
! *  water column thickness (real) is the water column thickness (from ice base 
! *  to ocean floor).  
! *  
! *  water column thickness scaling (real) is an upper limit for the water column 
! *  thickness to scale the basal melt
! *  
! *  TODO: check for sea level defined in sif instead of assuming z=0 is sea level
! *  
! *  
!  example sif solver section (old one from approx. 2015):
!
!Solver 3
!  equation = SubShelfMeltCalculation
!  procedure = "ElmerIceSolvers" "SubShelfMelt"
!  variable = meltRate
!  melt function = String bgf_scaled
!  lower surface variable name = String Zb
!  bedrock variable name = String bedrock
!  grounded mask name = string GroundedMask
!  grounding line melt = logical True
!  omega = Real 0.6
!  far field ocean temperature = Real 2.0
!  water column scaling = logical true
!  water column scaling factor = Real 200.0
!  ice base scaling depth = Real 100.0
!  ice base scaling offset = Real 40.0
!  iceBaseScaling = logical True
!End
!
! example sif solver section (newer one for ISMIP6; July 2023):
!
!Solver 4
!   Exec Interval = 12 
!   Equation = "Get Meltrate"
!   Procedure = "SubShelfMelt" "SubShelfMelt"
!   Variable = bmb 
!   Variable DOFs = 1 
!   melt function = string "ismip6"
!   lower surface variable name = string "Zb"
!   bedrock variable name = string "bedrock"
!   grounded mask name = String "GroundedMask"
!   grounding line melt = logical False
!   water column scaling = logical True
!   water column scaling factor = real 75.0
!End  
!
! Here's an example for the "no curtain" setup (added November 2023):
!
!Solver #
!  Exec Solver = "before timestep"
!  Equation = "Shelf melt"
!  Procedure = "SubShelfMelt" "SubShelfMelt"
!  Variable = bmb 
!  Variable DOFs = 1 
!  melt function = string "linear"
!  minimum bmb = real 0.0
!  maximum bmb = real 25.0
!  depth threshold = real 700.0
!  lower surface variable name = string "Zb"
!  grounded mask name = String "GroundedMask"
!  grounding line melt = logical False
!End

 SUBROUTINE SubShelfMelt (Model, Solver, dt, Transient)

  USE DefUtils
  
  IMPLICIT NONE

  ! in args
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt ! timestep
  LOGICAL :: Transient

  ! variables to be found in the elmer data structures
  TYPE(Variable_t), POINTER :: z_iceBase, z_bedrock, groundedMask, abmb_var, TimeVar
  TYPE(Variable_t), POINTER :: TF_var, deltaTF_var 
  INTEGER, POINTER :: TF_Perm(:), deltaTF_Perm(:), abmb_perm(:)

  ! parameters to be read in from this solvers section in the sif 
  CHARACTER(LEN=MAX_NAME_LEN) :: meltFunc         ! which melt function to be used
  CHARACTER(LEN=MAX_NAME_LEN) :: lowerSurfName    ! name of the variable within Elmer
  CHARACTER(LEN=MAX_NAME_LEN) :: groundedMaskName ! name of the variable within Elmer
  CHARACTER(LEN=MAX_NAME_LEN) :: bedrockName      ! " "  (bedrock = bathymmetry under the ice shelf)
  REAL (KIND=dp)              :: omega            ! melt rate tuning parameter
  REAL (KIND=dp)              :: T_far            ! the far field ocean temperature
  LOGICAL                     :: glMelt, wct_sc
  
  REAL (KIND=dp)           :: anomFactor, preFactor
  REAL (KIND=dp)           :: TF, DeltaTF, gamma0
  REAL(KIND=dp), POINTER   :: abmb_vals(:), TF_vals(:), deltaTF_vals(:)
  LOGICAL                  :: applyAnomaly

  ! an optional additional scaling factor using depth of base of ice shelf
  REAL (KIND=dp)              :: iceBaseOffset, iceBaseRef, iceBaseFactor
  LOGICAL                     :: iceBaseScaling

  ! truly local variables!
  TYPE(ValueList_t), POINTER :: SolverParams
  REAL (KIND=dp) :: wct, wct_factor          ! water column thickness
  REAL (KIND=dp) :: meltScaling, meltRate, bmb_min, bmb_max, depth_thresh 
  REAL (KIND=dp) :: rhoi, Lf, rhoo, SWCp, gammaT, cc, cnst, secondstoyear, deltaT, T_freeze
  LOGICAL        :: found
  INTEGER        :: ii
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName = 'SubShelfMelt'

  SolverParams => GetSolverParams()

  meltFunc = GetString( SolverParams, 'melt function',  Found )
  IF (.NOT. Found) THEN
     CALL FATAL(SolverName,'No keyword >melt function< found in SubShelfMelt solver')
  ELSE
     WRITE(Message,'(A,A)') 'Melt function is ',meltFunc
     CALL INFO(SolverName, Message,Level=2)
  END IF

  iceBaseScaling = GetLogical( SolverParams, 'ice base scaling',  Found )
  IF (Found) THEN
     IF (iceBaseScaling) THEN
        CALL INFO(SolverName, 'Using icebase scaling', Level=4)
        iceBaseOffset = GetConstReal( SolverParams, 'ice base scaling offset',  Found )
        IF (.NOT.Found) THEN
           CALL INFO(SolverName, 'Setting ice base offset to zero', Level=4)
           iceBaseOffset = 0.0_dp
        END IF
        iceBaseRef = GetConstReal( SolverParams, 'ice base scaling depth',  Found )
        IF (.NOT.Found) THEN
           CALL FATAL(SolverName,'ice base scaling used but >ice base scaling depth< not set')
        END IF
     ELSE
        CALL INFO(SolverName, 'Not using icebase scaling', Level=4)
     END IF
  ELSE
     CALL INFO(SolverName, 'Not using icebase scaling', Level=4)
     iceBaseScaling = .FALSE.
  END IF

  wct_sc = GetLogical( SolverParams, 'water column scaling',  Found )
  IF (Found) THEN
     wct_factor = GetConstReal( SolverParams, 'water column scaling factor',  Found )
     IF (.NOT.Found) THEN
        CALL FATAL(SolverName,'Water column scaling factor is required')
     END IF
  ELSE
     CALL INFO(SolverName, 'Not using water column scaling', Level=4)
  END IF

  rhoi =  GetConstReal( CurrentModel % Constants,'Ice Density',Found)
  IF (.NOT.Found) CALL FATAL(Solvername,'Ice Density not found in Constants')
  rhoo =  GetConstReal( CurrentModel % Constants,'Ocean Water Density',Found)
  IF (.NOT.Found) CALL FATAL(Solvername,'Ocean Water Density not found in Constants')
  Lf =  GetConstReal( CurrentModel % Constants,'Latent Heat SI',Found)
  IF (.NOT.Found) CALL FATAL(Solvername,'Latent Heat SI not found in Constants')
  SWCp =  GetConstReal( CurrentModel % Constants,'SW Cp',Found) ! specific heat of sea water
  IF (.NOT.Found) CALL FATAL(Solvername,'SW Cp not found in Constants')
  
  SELECT CASE(meltFunc)
  CASE('ISMIP6','ismip6')
     gamma0 =  GetConstReal( CurrentModel % Constants,'gamma 0',Found)
     IF (.NOT.Found) CALL FATAL(Solvername,'gamma 0 not found in Constants')
     ! Note: this is specifically the InitMIP anomaly, see hard coded time dependency further down
     applyAnomaly = GetLogical( SolverParams,'Apply Anomaly',Found)
     
     IF (applyAnomaly) THEN
        abmb_var => VariableGet(Model %Mesh % Variables,"abmb")
        IF (ASSOCIATED(abmb_var)) THEN
           abmb_vals => abmb_var % Values
           abmb_Perm => abmb_var % Perm
           WRITE(Message,*) "Found >abmb<, applying InitMIP bmb anomaly"
           CALL INFO(SolverName, Message, Level=3)
        ELSE
           applyAnomaly = .FALSE.
           WRITE(Message,*) "No var >abmb<, not applying InitMIP bmb anomaly"
           CALL INFO(SolverName, Message, Level=3)
        END IF
     END IF
     
     TF_var => VariableGet(Model %Mesh % Variables,"t_forcing")
     IF (ASSOCIATED(TF_var)) THEN
        TF_vals => TF_var % Values
        TF_Perm => TF_var % Perm
     ELSE
        CALL FATAL(SolverName,'Cant find var t_forcing')
     END IF
    
     deltaTF_var => VariableGet(Model %Mesh % Variables,"deltat_basin")
     IF (ASSOCIATED(deltaTF_var)) THEN
        deltaTF_vals => deltaTF_var % Values
        deltaTF_Perm => deltaTF_var % Perm
        CALL INFO(SolverName,'read var deltat_basin successfully')
     ELSE
        CALL FATAL(SolverName,'Cant find var deltat_basin:')
     END IF
     
  CASE('bgf')  
     ! some hard coded things specific to Rupert's old settings (as in 2017 TC paper)
     gammaT = 1.0e-4
     cc = 7.61e-4          ! Freezing temp change with depth (T_freeze decrease with zeta)
     cnst = (rhoo*SWCp*gammaT)/(rhoi*Lf)
     secondstoyear = 60.0_dp * 60.0_dp * 24.0_dp * 365.25_dp
     ! and some user defined settings...
     omega = GetConstReal( SolverParams, 'omega',  Found )
     IF (.NOT. Found) THEN
        CALL FATAL(SolverName,'No keyword >omega< found in SubShelfMelt solver')
     END IF    
     T_far = GetConstReal( SolverParams, 'far field ocean temperature',  Found )
     IF (.NOT. Found) THEN
        CALL FATAL(SolverName,'No keyword >far field ocean temperature< found in SubShelfMelt solver')
     END IF

  CASE ('linear','Linear')
     
     bmb_min      = GetConstReal( SolverParams, 'minimum bmb',  Found )
     IF (.NOT. Found) THEN
        CALL FATAL(SolverName,'No keyword >minimum bmb< found in SubShelfMelt solver')
     END IF
     bmb_max      = GetConstReal( SolverParams, 'maximum bmb',  Found )
     IF (.NOT. Found) THEN
        CALL FATAL(SolverName,'No keyword >maximum bmb< found in SubShelfMelt solver')
     END IF
     depth_thresh = GetConstReal( SolverParams, 'depth threshold',  Found )
     IF (.NOT. Found) THEN
        CALL FATAL(SolverName,'No keyword >depth threshold< found in SubShelfMelt solver')
     END IF
     
  CASE DEFAULT
     CALL FATAL(SolverName,'melt function not recognised')
     
  END SELECT
  
  lowerSurfName = GetString( SolverParams, 'lower surface variable name',  Found )
  IF (.NOT. Found) THEN
     CALL FATAL(SolverName,'No keyword >lower surface variable name< found in SubShelfMelt solver')
  END IF

  IF (wct_sc) THEN
     bedrockName = GetString( SolverParams, 'bedrock variable name',  Found )
     IF (.NOT. Found) THEN
        CALL FATAL(SolverName,'No keyword >bedrock variable name< found in SubShelfMelt solver')
     END IF
  END IF
  
  groundedMaskName = GetString( SolverParams, 'grounded mask name',  Found )
  IF (.NOT. Found) THEN
     CALL FATAL(SolverName,'No keyword >grounded mask name< found in SubShelfMelt solver')
  END IF
  
  glMelt = GetLogical( SolverParams, 'grounding line melt',  Found )
  IF (.NOT. Found) THEN
     CALL FATAL(SolverName,'No keyword >grounding line melt< found in SubShelfMelt solver')
  END IF

  ! find Elmer variables essential to the melt rate calculation
  z_iceBase => VariableGet( Solver % Mesh % Variables, TRIM(lowerSurfName))
  IF (.NOT. ASSOCIATED(z_iceBase)) THEN
     CALL FATAL(SolverName,'Failed to find ice base variable') !TODO: report variable name
  END IF
  IF (wct_sc) THEN
     z_bedrock => VariableGet( Solver % Mesh % Variables, TRIM(bedrockName))
     IF (.NOT. ASSOCIATED(z_bedrock)) THEN
        CALL FATAL(SolverName,'Failed to find bedrock variable') !TODO: report variable name
     END IF
  END IF
  groundedMask => VariableGet( Solver % Mesh % Variables, TRIM(groundedMaskName))
  IF (.NOT. ASSOCIATED(groundedMask)) THEN
     CALL FATAL(SolverName,'Failed to find grounded mask variable') !TODO: report variable name
  END IF

! TODO: check perm is associated for these vars
!  IF ( ((.NOT. ASSOCIATED(Solver % Variable % Perm)) .OR. &
!       (.NOT. ASSOCIATED(z_iceBase % Perm(i)))) .OR.      &
!       (.NOT. ASSOCIATED(z_bedrock % Perm(i))) ) THEN
!     CALL FATAL(SolverName,'permutations not associated where expected at the lower surface...')
!  END IF

  ! set no melting over whole domain.  
  ! Melt will be calculated for the ungrounded portion of the lower surface.
  Solver%Variable%Values = 0.0_dp

  DO ii=1,Model % NumberOfNodes

     ! the permutation is only non-zero for the body on which this solver is defined 
     ! (which should be the lower ice surface)
     IF (Solver%Variable%Perm(ii).LE.0) THEN
        CYCLE
     END IF

     ! only apply melt under the floating shelf ...
     IF (groundedMask%values(groundedMask%Perm(ii)) .GT. 0) THEN
        CYCLE
        ! ...and at the grounding line if desired
     ELSE IF (groundedMask%values(groundedMask%Perm(ii)) .EQ. 0) THEN
        IF (.NOT. glMelt) THEN
           CYCLE
        END IF
     END IF

     ! calculate the sub shelf melt rate
     SELECT CASE (meltFunc)

     CASE('bgf')
        T_freeze =  -1.85_dp + cc * z_iceBase%values(z_iceBase%Perm(ii)) ! In situ freezing temperature
        deltaT = T_freeze - T_far ! Thermal driving
        meltRate = cnst * Omega * deltaT * secondstoyear

     CASE('ISMIP6','ismip6')
        prefactor = (rhoo * SWCp/(rhoi * Lf))**2.0_dp
        TF        = TF_vals(TF_Perm(ii))
        deltaTF   = deltaTF_vals(deltaTF_Perm(ii))
        meltrate  = gamma0 * prefactor * (MAX((TF + deltaTF),0.0_dp))**2.0_dp
        
        IF (applyAnomaly) THEN
           meltrate = meltrate - anomFactor*abmb_vals(abmb_perm(ii))
        END IF
        
     CASE ('linear','Linear')
        prefactor = z_iceBase%values(z_iceBase%Perm(ii))/depth_thresh
        IF (prefactor.LT.-1.0) prefactor = -1.0
        meltrate = bmb_min*(1.0-prefactor) + bmb_max*prefactor
        
     CASE DEFAULT
        CALL FATAL(SolverName,'Melt function not recognised in SubShelfMelt calculation')
     END SELECT

     ! apply scaling if required (motivated by avoiding high melt near the grounding line)
     IF (wct_sc) THEN
        wct = z_iceBase%values(z_iceBase%Perm(ii)) - z_bedrock%values(z_bedrock%Perm(ii)) 
        meltScaling = TANH((wct)/(wct_factor/exp(1.0_dp)))
     ELSE
        meltScaling = 1.0_dp
     END IF
     
     ! additional scaling may also be applied to reduce melting of thin shelves near sea level
     IF (iceBaseScaling) THEN
        iceBaseFactor = MAX(0.0_dp, &
             -TANH((z_iceBase%values(z_iceBase%Perm(ii))+iceBaseOffset)*exp(1.0_dp)/iceBaseRef))
     ELSE
        iceBaseFactor = 1.0_dp
     END IF
     meltScaling = meltScaling * iceBaseFactor
     
     Solver%Variable%Values(Solver%Variable%Perm(ii)) = meltRate * meltScaling
     
  END DO
  
  END SUBROUTINE SubShelfMelt
