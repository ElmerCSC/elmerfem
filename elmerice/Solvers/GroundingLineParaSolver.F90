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
  REAL(KIND=dp) :: comp, GLNodeX, FFNodeX, GLstressSum, FFstressSum, cond
  REAL(KIND=dp) :: GLParaPosition, ratio

  INTEGER, POINTER :: Permutation(:), GroundedMaskPerm(:)
  INTEGER, POINTER :: NormalPerm(:), ResidPerm(:), HydroPerm(:)

  INTEGER :: DIM, HydroDIM, tt, ii, jj, n, GLnodenumber

  LOGICAL:: FirstTime = .TRUE., bedPComputed = .FALSE., UnFoundFatal


!=========================================================================
    ! TYPE(Nodes_t) :: tempElementNodes
    TYPE(Variable_t), POINTER :: CurrentTimeVar
    ! TYPE(variable_t), POINTER :: GroundedMaskVar, GroundingLineVar

    ! INTEGER, POINTER :: GroundedMaskPerm(:), GroundingLineParaPerm(:)
    INTEGER :: GLNodeIndex, FFNodeIndex
    REAL(KIND=dp) :: Time 

    LOGICAL :: GLparaSaveData = .FALSE., GotIt
    ! CHARACTER(LEN=MAX_NAME_LEN) :: Format
    CHARACTER(LEN=MAX_NAME_LEN) :: GLParaFileName='GLPressureData.dat'
    ! REAL(KIND=dp), POINTER :: LGParaData(:), FFParaData(:), GroundedMask(:), GroundingLinePara(:)
!=========================================================================

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

  ! Save pressure differences at GL element
  GLparaSaveData = ListGetLogical( Solver % Values, &
                      'Save Data for GL Parameterization', GotIt )
  IF ( .NOT. GotIt ) GLparaSaveData = .FALSE.

  IF (GLparaSaveData) THEN
    ! Get File name
    GLParaFileName= ListgetString( Solver % Values, 'Save Data File Name', GotIt )
    IF ( .NOT. GotIt ) GLParaFileName='GLPressureData.dat'
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

 
!======================= Compute parameterized GL position ==========================
! currently only work for 2D problem
!====================================================================================
      GLNodeX = 0.0
      FFNodeX = 0.0
      GLNodeIndex = 0
      FFNodeIndex = 0
      GLParaPosition = 0.0

      DO tt = 1, Model % NumberOfBoundaryElements
        Element => GetBoundaryElement(tt)
        IF (ParEnv % myPe .NE. Element % partIndex) CYCLE

        CALL GetElementNodes(Nodes, Element)

      ! For GL element which contains GL and FF nodes
        IF ( ALL(Permutation(Element % NodeIndexes) > 0) ) THEN
          IF (((ANY(VariableValues(Permutation(Element % NodeIndexes)) >= 0))) .AND. &
              ((ANY(VariableValues(Permutation(Element % NodeIndexes)) < 0)))) THEN
            n = GetElementNOFNodes(Element)
            ! Total stress on the elment
            GLstressSum = 0.0_dp
            FFstressSum = 0.0_dp
            DO ii = 1, n           
              GLnodenumber = Permutation(Element % NodeIndexes(ii))

              IF (GLnodenumber == 0) CYCLE


              cond = VariableValues(GLnodenumber)
              ! Check for GL nodes
              IF (cond < 0) THEN 
                GLstressSum = GLstressSum + VariableValues(GLnodenumber)
                GLNodeIndex = Element % NodeIndexes(ii)
                GLNodeX = Nodes % x(ii)
              ! Floating Nodes
              ELSE IF (cond >= 0) THEN
                FFstressSum = FFstressSum + VariableValues(GLnodenumber)
                FFNodeIndex = Element % NodeIndexes(ii)
                FFNodeX = Nodes % x(ii)
              END IF
            END DO

            IF ( (GLstressSum*FFstressSum) < 0.0 ) THEN
              ratio =  ABS(GLstressSum) / ( ABS(GLstressSum) + ABS(FFstressSum) )
              GLParaPosition = GLNodeX + ratio * ABS(FFNodeX - GLNodeX)
              WRITE (Message, '(A, g15.10)') '============== GL parameterization position at x =', GLParaPosition
              CALL Info(SolverName, Message, Level=3)
              CALL ListAddConstReal( Model % Constants, 'GroundingLine Position', GLParaPosition )
          
              IF (GLparaSaveData) THEN
                ! Get the current Time
                CurrentTimeVar => VariableGet( Solver % Mesh % Variables, 'Time')
                Time = CurrentTimeVar % Values(1) 
                    
                ! Save Data
                OPEN(unit=134, file=GLParaFileName, POSITION='APPEND')

                WRITE(134, *) Time, FFNodeIndex, GLNodeIndex, GLParaPosition, &
                      VariableValues(Permutation(FFNodeIndex)), VariableValues(Permutation(GLNodeIndex))

                CLOSE(134)
              END IF
            ELSE
              CALL Fatal(SolverName, 'GL parameterization error!')
            END IF
          END IF
        END IF
      END DO


      
!=====================================================================================





  CALL INFO( SolverName , 'Done')
 

END SUBROUTINE GroundingLineParaSolver 

