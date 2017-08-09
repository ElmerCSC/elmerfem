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
! *  Calculate sub ice shelf melting as a function of depth.  A number of 
! *  parameters may determine the melt function and are specificed in the sif.
! *  
! *  sub shelf melt function (string) is the name of one of the hard coded melt 
! *  functions in this routine.  Currently available are:
! *  bgf - melt as a function of pressure and thermal forcing created by Ben  
! *  Galton-Fenzi.
! *  bgf_scaled - same as bgf but scales down melt near grounding line based on 
! *  water column thickness
! *  
! *  bgf_scaled requires the following properties (in solver sif section):
! *  
! *  water column thickness (real) is the water column thickness (from ice base 
! *  to ocean floor).  
! *  
! *  water column thickness scaling (real) is an upper limit for the water column 
! *  thickness to scale the basal melt
! *  
! *  TODO: implement other melt functions - Walker? MISMIP+?
! *  TODO: check for sea level defined in sif instead of assuming z=0 is sea level
! *  TODO: use physical constants set in the sif rather than repeating them here
! *  TODO: make omega a spatially varying field... or some way of defining basins...
! *  
! *  
! *  
! *  examples sif solver section:
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
!  water column scaling = Real 200.0
!  ice base scaling depth = Real 100.0
!  ice base scaling offset = Real 40.0
!  iceBaseScaling = logical True

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
  TYPE(Variable_t), POINTER :: z_iceBase, z_bedrock, groundedMask

  ! parameters to be read in from this solvers section in the sif 
  CHARACTER(LEN=MAX_NAME_LEN) :: meltFunc         ! which melt function to be used
  CHARACTER(LEN=MAX_NAME_LEN) :: lowerSurfName    ! name of the variable within Elmer
  CHARACTER(LEN=MAX_NAME_LEN) :: groundedMaskName ! name of the variable within Elmer
  CHARACTER(LEN=MAX_NAME_LEN) :: bedrockName      ! " "  (bedrock = bathymmetry under the ice shelf)
  REAL (KIND=dp)              :: omega            ! melt rate tuning parameter
  REAL (KIND=dp)              :: T_far            ! the far field ocean temperature
  REAL (KIND=dp)              :: wct_sc           ! a scaling factor using wct (water column thickness) 
  LOGICAL                     :: glMelt           ! determines whether melt occurs actually at the grounding line itself

  ! an optional additional scaling factor using depth of base of ice shelf
  REAL (KIND=dp)              :: iceBaseOffset, iceBaseRef, iceBaseFactor
  LOGICAL                     :: iceBaseScaling

  ! truly local variables!
  TYPE(ValueList_t), POINTER :: SolverParams
  REAL (KIND=dp) :: wct             ! water column thickness
  REAL (KIND=dp) :: meltScaling, meltRate 
  REAL (KIND=dp) :: rhoi, L, rhow, cw, gammaT, c, cnst, secondstoyear, deltaT, T_freeze
  LOGICAL        :: found
  INTEGER        :: i
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName = 'SubShelfMelt'

  SolverParams => GetSolverParams()


  ! Physical constants
  rhoi = 916.0_dp
  L = 3.34e5
  rhow = 1025.0_dp
  cw = 3148.0_dp
  gammaT = 1.0e-4
  c = 7.61e-4          ! Freezing temp change with depth (T_freeze decrease with zeta)
  cnst = (rhow*cw*gammaT)/(rhoi*L)
  secondstoyear = 60.0_dp * 60.0_dp * 24.0_dp * 365.25_dp
  
  ! read parameters from sif keywords

  meltFunc = GetString( SolverParams, 'melt function',  Found )
  IF (.NOT. Found) THEN
     CALL FATAL(SolverName,'No keyword >melt function< found in SubShelfMelt solver')
  ELSE
     WRITE(Message,'(A,A)') 'Melt function is ',meltFunc
     CALL INFO(SolverName, Message,Level=2)
  END IF
  
  lowerSurfName = GetString( SolverParams, 'lower surface variable name',  Found )
  IF (.NOT. Found) THEN
     CALL FATAL(SolverName,'No keyword >lower surface variable name< found in SubShelfMelt solver')
  END IF
  
  bedrockName = GetString( SolverParams, 'bedrock variable name',  Found )
  IF (.NOT. Found) THEN
     CALL FATAL(SolverName,'No keyword >bedrock variable name< found in SubShelfMelt solver')
  END IF
  
  groundedMaskName = GetString( SolverParams, 'grounded mask name',  Found )
  IF (.NOT. Found) THEN
     CALL FATAL(SolverName,'No keyword >grounded mask name< found in SubShelfMelt solver')
  END IF
  
  omega = GetConstReal( SolverParams, 'omega',  Found )
  IF (.NOT. Found) THEN
     CALL FATAL(SolverName,'No keyword >omega< found in SubShelfMelt solver')
  END IF
  
  T_far = GetConstReal( SolverParams, 'far field ocean temperature',  Found )
  IF (.NOT. Found) THEN
     CALL FATAL(SolverName,'No keyword >far field ocean temperature< found in SubShelfMelt solver')
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

  wct_sc = GetConstReal( SolverParams, 'water column scaling',  Found )
  IF (.NOT. Found) THEN
     SELECT CASE (meltFunc)
     CASE('bgf_scaled')
        CALL FATAL(SolverName,'No keyword >water column scaling< found in SubShelfMelt solver')
     CASE('bgf')
     CASE DEFAULT
        CALL FATAL(SolverName,'Melt function not recognised whilst reading water column scaling')
     END SELECT
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
  z_bedrock => VariableGet( Solver % Mesh % Variables, TRIM(bedrockName))
  IF (.NOT. ASSOCIATED(z_bedrock)) THEN
     CALL FATAL(SolverName,'Failed to find bedrock variable') !TODO: report variable name
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

  DO i=1,Model % NumberOfNodes

     ! the permutation is only non-zero for the body on which this solver is defined 
     ! (which should be the lower ice surface)
     IF (Solver%Variable%Perm(i).LE.0) THEN
        CYCLE
     END IF

     ! only apply melt under the floating shelf ...
     IF (groundedMask%values(groundedMask%Perm(i)) .GT. 0) THEN
        CYCLE
        ! ...and at the grounding line if desired
     ELSE IF (groundedMask%values(groundedMask%Perm(i)) .EQ. 0) THEN
        IF (.NOT. glMelt) THEN
           CYCLE
        END IF
     END IF

     ! calculate the sub shelf melt rate
     SELECT CASE (meltFunc)
     CASE('bgf','bgf_scaled')
        T_freeze =  -1.85_dp + c * z_iceBase%values(z_iceBase%Perm(i)) ! In situ freezing temperature
        deltaT = T_freeze - T_far ! Thermal driving
        meltRate = cnst * Omega * deltaT * secondstoyear
     CASE DEFAULT
        CALL FATAL(SolverName,'Melt function not recognised in SubShelfMelt calculation')
     END SELECT
     
     ! apply scaling if required (motivated by avoiding high melt near the grounding line)
     SELECT CASE (meltFunc)
     CASE('bgf')
        meltScaling = 1.0_dp
     CASE('bgf_scaled')
        wct = z_iceBase%values(z_iceBase%Perm(i)) - z_bedrock%values(z_bedrock%Perm(i)) 
        meltScaling = TANH((wct)/(wct_sc/exp(1.0_dp)))
     CASE DEFAULT
        CALL FATAL(SolverName,'Melt function not recognised in SubShelfMelt scaling')
     END SELECT

     ! an additional scaling term may also be applied to reduce melting of thin shelves near sea level
     IF (iceBaseScaling) THEN
        iceBaseFactor = MAX(0.0_dp, &
             -TANH((z_iceBase%values(z_iceBase%Perm(i))+iceBaseOffset)*exp(1.0_dp)/iceBaseRef))

     ELSE
        iceBaseFactor = 1.0_dp
     END IF
     meltScaling = meltScaling * iceBaseFactor

     Solver%Variable%Values(Solver%Variable%Perm(i)) = - meltRate * meltScaling

  END DO

END SUBROUTINE SubShelfMelt
