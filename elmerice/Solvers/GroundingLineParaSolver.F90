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
! *  Authors: Cheng Gong
! *  Email:   cheng.gong@it.uu.se
! *  Web:      
! *
! *  Original Date: 2016-09-05
! * 
! *****************************************************************************
!> Solver for Gronging Line parameterization according to the differences of loads 
!> at the bottom of ice.

SUBROUTINE GroundingLineParaSolver( Model,Solver,dt,TransientSimulation )

!------------------------------------------------------------------------------
!******************************************************************************
! 
!
!******************************************************************************
  USE ElementDescription
  USE DefUtils

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t), INTENT(INOUT)     :: Solver
  TYPE(Model_t), INTENT(INOUT)      :: Model
  REAL(KIND=dp), INTENT(IN)         :: dt
  LOGICAL, INTENT(IN)               :: TransientSimulation

  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------

  TYPE(Element_t), POINTER :: Element
  TYPE(ValueList_t), POINTER :: Material, SolverParams
  TYPE(Variable_t), POINTER :: PointerToVariable
  TYPE(variable_t), POINTER :: NormalVar, VarSurfResidual, GroundedMaskVar, HydroVar
  TYPE(Nodes_t), SAVE :: Nodes

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName = 'GroundingLinePara'
  REAL(KIND=dp), POINTER :: VariableValues(:), GroundedMask(:)
  REAL(KIND=dp), POINTER :: NormalValues(:), ResidValues(:), HydroValues(:)
  REAL(KIND=dp), ALLOCATABLE :: Normal(:), Fwater(:), Fbwater(:), Fbase(:)
  REAL(KIND=dp) :: comp

  INTEGER, POINTER :: Permutation(:), GroundedMaskPerm(:)
  INTEGER, POINTER :: NormalPerm(:), ResidPerm(:), HydroPerm(:)

  INTEGER :: DIM, HydroDIM, tt, ii, jj, n

  LOGICAL:: FirstTime = .TRUE., bedPComputed = .FALSE., UnFoundFatal


  SAVE HydroDIM, bedPComputed, DIM, FirstTime
  SAVE Normal, Fwater, Fbwater, Fbase
  !------------------------------------------------------------------------------
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! First time step for the First time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF (FirstTime) THEN
    ! Read in Solver Parameter 
    SolverParams => GetSolverParams ()
    IF ( .NOT. ASSOCIATED(SolverParams) ) &
          CALL FATAL(SolverName, 'No Solver Section found' )

    DIM = CoordinateSystemDimension()
    FirstTime = .FALSE.
    ALLOCATE( Normal(DIM), Fwater(DIM), Fbwater(DIM), Fbase(DIM) )

    bedPComputed = GetLogical( SolverParams, 'Compute Bed Pressure' )
  END IF

  ! Initialize temp arrays
  Normal = 0.0_dp
  Fwater = 0.0_dp
  Fbwater = 0.0_dp
  Fbase = 0.0_dp

  ! Pointer for the current solver
  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values

  ! Ground Mask
  GroundedMaskVar => VariableGet( Model % Mesh % Variables, 'GroundedMask', UnFoundFatal=UnFoundFatal)
  GroundedMask => GroundedMaskVar % Values
  GroundedMaskPerm => GroundedMaskVar % Perm

  ! Load from flow solver
  VarSurfResidual => VariableGet( Model % Mesh % Variables, 'Flow Solution Loads',UnFoundFatal=UnFoundFatal)
  ResidPerm => VarSurfResidual  % Perm
  ResidValues => VarSurfResidual % Values

  ! Normal Vector
  NormalVar => VariableGet(Model % Variables,'Normal Vector',UnFoundFatal=UnFoundFatal)
  NormalPerm => NormalVar % Perm
  NormalValues => NormalVar % Values

  ! Water Pressure
  HydroVar => VariableGet( Model % Mesh % Variables, 'Fw',UnFoundFatal=UnFoundFatal)
  HydroPerm => HydroVar % Perm
  HydroValues => HydroVar % Values

  ! Check for bed pressure
  IF (((2*DIM) == HydroVar % dofs) .OR. (bedPComputed) ) THEN
    HydroDIM = 2 * DIM
    bedPComputed = .TRUE.
  ELSE 
    HydroDIM = DIM
    bedPComputed = .FALSE.
  END IF

  ! Loop over each boundary element
  DO tt = 1, Model % NumberOfBoundaryElements   
    ! Get boundary element
    Element => GetBoundaryElement(tt)
    IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
    n = GetElementNOFNodes()

    ! Get Nodes of the current element 
    CALL GetElementNodes( Nodes )
        
    ! Look up only within bottom element with groundmasks  
    IF (ANY(GroundedMaskPerm(Element % NodeIndexes(1:n))==0)) CYCLE

    ! For each node
    DO ii = 1, n
      jj = Element % NodeIndexes(ii)
      ! Normal vectors on the bottom
      Normal = NormalValues(DIM*(NormalPerm(jj)-1)+1 : DIM*NormalPerm(jj))
      ! Water Pressure on the bottom of ice
      Fwater = HydroValues(HydroDIM*(HydroPerm(jj)-1)+1 : HydroDIM*(HydroPerm(jj)-1)+DIM)
      ! Ice load on the bottom 
      Fbase = ResidValues((DIM+1)*(ResidPerm(jj)-1)+1 : (DIM+1)*ResidPerm(jj)-1)

      ! comparison between water pressure and bed action
      comp = ABS( SUM( Fwater * Normal ) ) - ABS( SUM( Fbase * Normal ) )

      ! Compute water pressure at the bedrock, it could be different from Fwater as
      ! the node is floating
      IF (bedPComputed) THEN
        Fbwater = HydroValues(HydroDIM*(HydroPerm(jj)-1)+DIM+1 : HydroDIM*HydroPerm(jj))
        comp = ABS( SUM( Fbwater * Normal ) ) - ABS( SUM( Fbase * Normal ) )
      END IF
      ! Save the difference of loads
      VariableValues( Permutation(Element % NodeIndexes(ii)) ) = comp
    END DO
     

  END DO

  IF ( ParEnv % PEs>1 ) CALL ParallelSumVector( Solver % Matrix, VariableValues, 1 )
 
  CALL INFO( SolverName , 'Done')
 

END SUBROUTINE GroundingLineParaSolver 

