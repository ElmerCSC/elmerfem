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
! *  Authors: Thomas Zwinger, Denis Cohen, Juha Hartikainen
! *  Email:  thomas Zwinger [at] csc.fi 
! *  Web:     http://elmerice.elmerfem.org
! *  Address: CSC - Scientific Computing Ltd.  
! *               Keilaranta 14                    
! *               02101 Espoo, Finland             
! *                                                 
! *       Original Date:  January 2017  -               
! * 
! *****************************************************************************
!>  Solvers und utilities needed for solid ground deformation in permafrost model
!---------------------------------------------------------------------------------------------
!==============================================================================
!>  initialization of Porosity to given reference value in material
!> \ingroup Solvers
!==============================================================================
SUBROUTINE PorosityInit(Model, Solver, Timestep, TransientSimulation )
  !==============================================================================

  USE DefUtils
  USE PermaFrostMaterials
  IMPLICIT NONE

  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: Timestep
  LOGICAL :: TransientSimulation

  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: CurrentElement
  TYPE(Variable_t), POINTER :: PorosityVariable
  TYPE(ValueList_t), POINTER :: SolverParams,Material
  INTEGER, POINTER :: PorosityPerm(:), NodeIndexes(:)
  REAL(KIND=dp), POINTER :: PorosityValues(:)
  INTEGER(KIND=dp), ALLOCATABLE :: NodalHits(:)
  INTEGER :: DIM, i, j, k, N, NumberOfRockRecords,RockMaterialID,CurrentNode,Active,totalunset,totalset
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName="PorosityInit"
  CHARACTER(LEN=MAX_NAME_LEN) :: PorosityName,ElementRockMaterialName
  LOGICAL :: Visited = .FALSE., Found, GotIt,ElementWiseRockMaterial,IsNodalVariable=.FALSE.

  SAVE Visited,ElementWiseRockMaterial,IsNodalVariable,NumberOfRockRecords
  !,DIM,NumberOfRockRecords

  !------------------------------------------------------------------------------

  ! Execute solver only once at beginning
  IF (Visited) RETURN

  CALL INFO(SolverName, '-----------------------------------', Level=4)
  CALL INFO(SolverName, 'Initializing porosity to reference ', Level=4)
  CALL INFO(SolverName, 'levels in material file            ', Level=4)
  CALL INFO(SolverName, '-----------------------------------', Level=4)

  ! Get variables
  DIM = CoordinateSystemDimension()

  ! Get info
  SolverParams => GetSolverParams()

  PorosityName = ListGetString(SolverParams, &
       'Porosity Variable', GotIt )
  IF (.NOT.GotIt) THEN
    PorosityName = "Porosity"
    CALL WARN(SolverName, ' "Porosity Variable" not found - trying default "Porosity"')
  END IF
  PorosityVariable => VariableGet( Solver % Mesh % Variables, PorosityName,GotIt )

  IF ( ASSOCIATED( PorosityVariable ) ) THEN
    PorosityPerm    => PorosityVariable % Perm
    PorosityValues  => PorosityVariable % Values
    CALL INFO(SolverName,'Found porosity variable:'//TRIM(PorosityName),Level=7)
  ELSE
    CALL FATAL(SolverName, 'Could not find "Porosity Variable"')
  END IF

  IsNodalVariable = GetLogical(SolverParams, &
       'Nodal Porosity', GotIt)
  PRINT *,"IsNodalVariable", IsNodalVariable, GotIt
  IF (.NOT.GotIt) THEN
    CALL WARN(SolverName,'Keyword "Nodal Porosity" not found. Assuming element-wise porosity variable')
  ELSE IF (IsNodalVariable) THEN
    CALL INFO(SolverName,'Assigning porosity to nodal variable',Level=3)
    ALLOCATE(NodalHits(Solver % Mesh % NumberOfNodes))
    NodalHits = 0
  ELSE
    CALL INFO(SolverName,'Assigning porosity to element-wise variable',Level=3)
  END IF
  !==============================================================================
  ! Loop over elements
  Active = Solver % NumberOFActiveElements
  PorosityValues = 0.0000001_dp ! some not complete insane default

  DO i = 1, Active
    CurrentElement => GetActiveElement(i)

    !PRINT *,"PorosityInit", i,"/",Active
    
    IF (ParEnv % myPe .NE. CurrentElement % partIndex) CYCLE
    N = GetElementNOFNodes(CurrentElement)
    NodeIndexes => CurrentElement % NodeIndexes
    Material => GetMaterial(CurrentElement)
    
    IF (.NOT.ASSOCIATED(Material)) CALL FATAL(SolverName,'No Material pointer found')
    IF (.NOT.Visited) THEN
      ! check, whether we have globally or element-wise defined values of rock-material parameters
      ElementRockMaterialName = &
           ListGetString(Material,"Element Rock Material File",ElementWiseRockMaterial)
      !PRINT *,"PorosityInit:",TRIM(ElementRockMaterialName),ElementWiseRockMaterial
      IF (ElementWiseRockMaterial) THEN
        WRITE (Message,*) 'Found "Element Rock Material File"'
        CALL INFO(SolverName,Message,Level=5)
        CALL INFO(SolverName,'Using element-wise rock material definition',Level=5)
        ! read element-wise material parameter (GlobalRockMaterial will have one entry each element)
        NumberOfRockRecords = &
             ReadPermafrostElementRockMaterial(ElementRockMaterialName,Solver,DIM)
        PRINT *, "Partition:", ParEnv % myPe, "NumberOfRockRecords:", NumberOfRockRecords
      ELSE
        NumberOfRockRecords =  ReadPermafrostRockMaterial( Material )
        WRITE (Message,*) "NumberOfRockRecords in File: ", TRIM(ElementRockMaterialName), ":", NumberOfRockRecords
        CALL INFO(SolverName, Message, Level=5)
      END IF
      IF (NumberOfRockRecords < 1) THEN
        CALL FATAL(SolverName,'No Rock Material specified')
      ELSE
        CALL INFO(SolverName,'Permafrost Rock Material read',Level=6)
      END IF
      Visited=.True.
    END IF

    IF (ElementWiseRockMaterial) THEN
      RockMaterialID = i
    ELSE      
      RockMaterialID = ListGetInteger(Material,'Rock Material ID', GotIt,UnfoundFatal=.TRUE.)
      IF (.NOT.GotIt) CALL FATAL(SolverName,"Rock Material ID not found")
    END IF

    IF (IsNodalVariable) THEN
      DO j=1,N
        CurrentNode = PorosityPerm(NodeIndexes(j))
        IF (CurrentNode <= 0) CYCLE
        PorosityValues(CurrentNode) = PorosityValues(CurrentNode) + GlobalRockMaterial % eta0(RockMaterialID)
        NodalHits(CurrentNode) = NodalHits(CurrentNode) + 1
      END DO
    ELSE
      PorosityValues(PorosityPerm(i)) = GlobalRockMaterial % eta0(RockMaterialID)
    END IF
  END DO

  IF (IsNodalVariable) THEN
    DO j = 1,Model % Mesh % NumberOfNodes
      CurrentNode = PorosityPerm(j)
      IF (CurrentNode <= 0) CYCLE
      PorosityValues(CurrentNode) = PorosityValues(CurrentNode)/NodalHits(CurrentNode)    
    END DO
    !DEALLOCATE(NodalHits)
  END IF
  
  CALL INFO(SolverName, '-----------------------------------', Level=4)
  CALL INFO(SolverName, 'Initializing porosity done         ', Level=4)
  CALL INFO(SolverName, '-----------------------------------', Level=4)
END SUBROUTINE PorosityInit
!==============================================================================
!>  Evolution of Porosity
!> \ingroup Solvers
!==============================================================================
SUBROUTINE PermafrostPorosityEvolution( Model, Solver, Timestep, TransientSimulation )
 USE DefUtils
  USE PermaFrostMaterials
  IMPLICIT NONE

  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: Timestep
  LOGICAL :: TransientSimulation

  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: CurrentElement
  TYPE(Variable_t), POINTER :: PorosityVariable, TemperatureVar, PressureVar, StrainVar
  TYPE(ValueList_t), POINTER :: SolverParams,Material
  INTEGER, POINTER :: PorosityPerm(:), StrainPerm(:), TemperaturePerm(:), PressurePerm(:), NodeIndexes(:)
  REAL(KIND=dp), POINTER :: PorosityValues(:), Strain(:), Temperature(:), Pressure(:),&
       NodalStrain(:), NodalTemperature(:), NodalPressure(:),&
       PrevNodalTemperature(:), PrevNodalPressure(:),&
       PrevTemperature(:), PrevPressure(:)
  REAL(KIND=dp), ALLOCATABLE :: PrevStrainInvariant(:)
  REAL(KIND=dp) :: aux, Nodalrhos, PrevNodalrhos, StrainInvariant,&
            GasConstant, N0, DeltaT, T0, p0, eps, Gravity(3)
  INTEGER :: DIM, i, j, k, N, NumberOfRockRecords,RockMaterialID,CurrentNode,Active,&
       StrainDOFs,TemperatureDOFS,PressureDOFs,totalunset,totalset,istat
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName="PermafrostPorosityEvolution"
  CHARACTER(LEN=MAX_NAME_LEN) :: PorosityName,PressureName,TemperatureName,StrainVarName,ElementRockMaterialName
  LOGICAL :: FirstTime=.TRUE.,FirstVisit=.TRUE.,Found,GotIt,ElementWiseRockMaterial,&
       StrainVarExists, TemperatureVarExists, PressureVarExists,ConstVal,ConstantTemp
  !------------------------------
  SAVE FirstTime,FirstVisit,ElementWiseRockMaterial,&
       NodalStrain, NodalTemperature, NodalPressure,&
       PrevNodalTemperature, PrevNodalPressure,&
       StrainVar, TemperatureVar, PressureVar,&
       Strain, PrevStrainInvariant, Temperature, Pressure,&
       StrainDOFs,TemperatureDOFs,PressureDOFs,&
       PrevTemperature, PrevPressure, &
       StrainPerm, TemperaturePerm, PressurePerm,&
       StrainVarExists, TemperatureVarExists, PressureVarExists, &
       DIM, N0, GasConstant, DeltaT, T0, p0, eps, Gravity,&
       NumberOfRockRecords, ConstantTemp

  !------------------------------------------------------------------------------
  CALL INFO(SolverName, '-----------------------------------', Level=4)
  CALL INFO(SolverName, ' computing evolution of porosity   ', Level=4)
  CALL INFO(SolverName, '-----------------------------------', Level=4)

  ! Get info and solver variable
  SolverParams => GetSolverParams()
  PorosityVariable => Solver % Variable
  IF (.NOT.ASSOCIATED(PorosityVariable)) CALL FATAL(SolverName,'No variable for "Porosity" associated')
  PorosityName = TRIM( Solver % Variable % Name )
  IF (TRIM(PorosityName) .NE. "porosity") THEN
    WRITE (Message,*) TRIM(PorosityName),' is not the expected "Porosity" - hopefully on purpose'
    CALL WARN(SolverName, Message)
  END IF

  PorosityPerm    => PorosityVariable % Perm
  PorosityValues  => PorosityVariable % Values
  IF (.NOT.ASSOCIATED( PorosityVariable % PrevValues )) THEN
    ALLOCATE(PorosityVariable % PrevValues(SIZE(PorosityValues), 1))
  END IF
  PorosityVariable % PrevValues(1:SIZE(PorosityValues),1) =&
       PorosityVariable % Values(1:SIZE(PorosityValues))
  
  IF(FirstTime) THEN
    IF(.NOT.(&
         ReadPermafrostConstants(Model, SolverName, DIM, GasConstant, N0, DeltaT, T0, p0, eps, Gravity)))&
         CALL FATAL(SolverName,"Errors in reading constants")
  END IF
  
  ! assign needed variables
  !(NB: we rather skip AssignVar routine, as ONLY Temperature and Pressure are needed)
  IF (FirstTime .OR. (Solver % Mesh % Changed ) ) THEN
    CALL INFO(SolverName,"Initialisation",Level=5)
    StrainVarName = GetString(SolverParams,'Strain Variable',Found)
    IF (.NOT.Found) THEN
      WRITE(StrainVarName,*) 'Strain'
      WRITE(Message,*) '"Strain Variable" not found - assuming default value: ',TRIM(StrainVarName)
      CALL WARN(SolverName,Message)
    ELSE
      WRITE(Message,*) '"Strain Variable" found and set to: ',TRIM(StrainVarName)
      CALL INFO(SolverName,Message,Level=5)
    END IF
    CALL AssignSingleVar(Solver,Model,NodalStrain,StrainVar,StrainPerm, Strain, &
         StrainVarName,StrainDOFs,StrainVarExists)
    IF (ASSOCIATED(StrainVar)) THEN
      IF (.NOT.FirstTime) DEALLOCATE(PrevStrainInvariant)
      ALLOCATE(PrevStrainInvariant(SIZE(StrainPerm)),stat=istat)
    ELSE
      CALL FATAL(SolverName,'No "Strain Varaible" associated')
    END IF
    ConstantTemp = GetLogical(SolverParams,'Constant Temperature',Found)
    TemperatureName = GetString(SolverParams,'Temperature Variable',Found)
    IF (.NOT.Found) THEN
      CALL WARN(SolverName,' "Temperature Variable" not found - assuming default value "Temperature" ')
      WRITE(TemperatureName,*) 'Temperature'
          ELSE
      WRITE(Message,*) '"Temperature Variable" found and set to: ',TRIM(TemperatureName)
      CALL INFO(SolverName,Message,Level=5)
    END IF
    IF (ConstantTemp) THEN
      CALL AssignSingleVar(Solver,Model,NodalTemperature,TemperatureVar,&
           TemperaturePerm, Temperature, &
           TemperatureName,TemperatureDOFS,TemperatureVarExists)
      ALLOCATE(PrevNodalTemperature(Solver % Mesh % MaxElementNodes))
    ELSE
      CALL AssignSingleVar(Solver,Model,NodalTemperature,TemperatureVar,&
           TemperaturePerm, Temperature, &
           TemperatureName,TemperatureDOFS,TemperatureVarExists,&
           PrevNodalVariable=PrevNodalTemperature, PrevVariable=PrevTemperature)
    END IF
    
    PressureName = GetString(SolverParams,'Pressure Variable',Found)
    IF (.NOT.Found) THEN
      CALL WARN(SolverName,' "Pressure Variable" not found - assuming default value "Pressure" ')
      WRITE(PressureName,*) 'Pressure'
    ELSE
      WRITE(Message,*) ' "Pressure Variable" found and set to: ',PressureName
      CALL INFO(SolverName,Message,Level=5)
    END IF
    CALL AssignSingleVar(Solver,Model,NodalPressure,PressureVar,&
         PressurePerm, Pressure, &
         PressureName,PressureDOFs,PressureVarExists,&
         PrevNodalVariable=PrevNodalPressure, PrevVariable=PrevPressure)

    FirstVisit = GetLogical(SolverParams,'Initialize Time Derivatives',Found)
    IF (.NOT.Found) FirstVisit = .FALSE.
  ELSE
    IF (ConstantTemp) &
         ALLOCATE(PrevNodalTemperature(Solver % Mesh % MaxElementNodes))
  END IF
  ! some sanity check
  IF (.NOT.ASSOCIATED(Temperature) .OR. .NOT.ASSOCIATED(TemperaturePerm))&
       CALL FATAL(SolverName,'Values of temperature variable not found')
  IF (.NOT.ConstantTemp .AND. .NOT.ASSOCIATED(PrevTemperature)) &
    CALL FATAL(SolverName,'Previous values of temperature variable not found') 
  IF (.NOT.ASSOCIATED(Pressure) .OR. .NOT.ASSOCIATED(PressurePerm))&
       CALL FATAL(SolverName,'Values of pressure variable not found')
  IF (.NOT.ASSOCIATED(PrevPressure))&
       CALL FATAL(SolverName,'Previous values of pressure variable not found')
  ! Loop over elements
  Active = Solver % NumberOFActiveElements

  DO i = 1, Active
    CurrentElement => GetActiveElement(i)
    NodeIndexes => CurrentElement % NodeIndexes
    Material => GetMaterial(CurrentElement)
    ConstVal = GetLogical(Material,'Constant Permafrost Properties',Found)
    IF (ConstVal) &
        CALL INFO(SolverName,'"Constant Permafrost Properties" set to true',Level=5)
    IF (.NOT.ASSOCIATED(Material)) CALL FATAL(SolverName,'No Material pointer found')
    IF (FirstTime) THEN
      ! check, whether we have globally or element-wise defined values of rock-material parameters
      ElementRockMaterialName = ListGetString(Material,"Element Rock Material File",ElementWiseRockMaterial)
      IF (ElementWiseRockMaterial) THEN
        CALL INFO(SolverName,'Found "Element Rock Material File"',Level=5)
        CALL INFO(SolverName,'Using element-wise rock material definition',Level=5)
      END IF
      IF (ElementWiseRockMaterial) THEN
        ! read element-wise material parameter (GlobalRockMaterial will have one entry each element)
        NumberOfRockRecords = &
             ReadPermafrostElementRockMaterial(ElementRockMaterialName,Solver,DIM)
        PRINT *, "NumberOfRockRecords", NumberOfRockRecords
      ELSE
        NumberOfRockRecords =  ReadPermafrostRockMaterial( Material )
      END IF
      IF (NumberOfRockRecords < 1) THEN
        CALL FATAL(SolverName,'No Rock Material specified')
      ELSE
        CALL INFO(SolverName,'Permafrost Rock Material read',Level=6)
      END IF
      dim = CoordinateSystemDimension()
      FirstTime = .FALSE.
    END IF

    IF (ElementWiseRockMaterial) THEN
      RockMaterialID = i
    ELSE      
      RockMaterialID = ListGetInteger(Material,'Rock Material ID', GotIt,UnfoundFatal=.TRUE.)
      IF (.NOT.GotIt) CALL FATAL(SolverName,"Rock Material ID not found")
    END IF
    N = GetElementNOFNodes(CurrentElement)
    CALL ReadSingleVar(N,CurrentElement,TemperaturePerm,NodalTemperature,Temperature,TemperatureDOFs)
    CALL ReadSingleVar(N,CurrentElement,PressurePerm,NodalPressure,Pressure,PressureDOFs)
    IF (FirstVisit) THEN ! write current values if starting
      IF (.NOT.ConstantTemp) THEN
        CALL ReadSingleVar(N,CurrentElement,TemperaturePerm,PrevNodalTemperature,Temperature,TemperatureDOFs)
      ELSE
        PrevNodalTemperature(1:N) = NodalTemperature(1:N)
      END IF
      CALL ReadSingleVar(N,CurrentElement,PressurePerm,PrevNodalPressure,Pressure,PressureDOFs)
    ELSE
      IF (.NOT.ConstantTemp) THEN
        CALL ReadSingleVar(N,CurrentElement,TemperaturePerm,PrevNodalTemperature,PrevTemperature,TemperatureDOFs)
      ELSE
        PrevNodalTemperature(1:N) = NodalTemperature(1:N)
      END IF
      CALL ReadSingleVar(N,CurrentElement,PressurePerm,PrevNodalPressure,PrevPressure,PressureDOFs)
    END IF

    ! Loop over nodes of element
    DO k = 1, N
      CurrentNode = CurrentElement % NodeIndexes(k)
      PrevNodalrhos = rhos(RockMaterialID,T0,p0,&
           PrevNodalTemperature(k),&
           PrevNodalPressure(k),ConstVal)
      IF (PrevNodalrhos .NE. PrevNodalrhos) THEN
        PRINT *,"PermafrostPorosityEvolution: ","Found weird number for PrevNodalrhos"
        PRINT *,"PermafrostPorosityEvolution: ",&
           PrevNodalTemperature(k),&
           PrevNodalPressure(k),&
           RockMaterialID,T0,p0,ConstVal,&
           CurrentElement % NodeIndexes(k)
        CALL FATAL(SolverName,'Exiting')
      END IF
      Nodalrhos = rhos(RockMaterialID,T0,p0,&
           NodalTemperature(k),&
           NodalPressure(k),ConstVal)
      IF (Nodalrhos .NE. Nodalrhos) THEN
        PRINT *,"PermafrostPorosityEvolution: ","Found weird number for Nodalrhos"
        PRINT *,"PermafrostPorosityEvolution: ",&
           NodalTemperature(k),&
           NodalPressure(k),&
           RockMaterialID,T0,p0,ConstVal,&
           CurrentElement % NodeIndexes(k)
        CALL FATAL(SolverName,'Exiting')
      END IF
      StrainInvariant = 0.0
      ! first 1..DIM elements of StrainRate variable are th ediagonal entries
      DO J=1,DIM
        PrevStrainInvariant(StrainPerm(CurrentNode)) = StrainInvariant
        StrainInvariant = StrainInvariant &
             + Strain((StrainPerm(CurrentNode)-1)*StrainDOFs + J)
      END DO
      aux = 1.0_dp - (PrevNodalrhos/Nodalrhos)*&
           (1.0_dp + PrevStrainInvariant(StrainPerm(CurrentNode)))/(1.0_dp + StrainInvariant)
      PorosityValues(PorosityPerm(CurrentNode)) = &
           PorosityVariable % PrevValues(PorosityPerm(CurrentNode),1)*(1.0_dp - aux) + aux
      IF (PorosityValues(PorosityPerm(CurrentNode)) <= 0.0) THEN
        WRITE(Message,*) "Reset Invalid value: ",&
             PorosityValues(PorosityPerm(CurrentNode)),&
             "=", PorosityVariable % PrevValues(PorosityPerm(CurrentNode),1),&
             "*", (1.0_dp - aux), "+", aux
        !CALL FATAL(SolverName, Message)
        CALL WARN(SolverName, Message)
        PorosityValues(PorosityPerm(CurrentNode)) = 1.0d-08 !!!!! REPLACE 
      END IF
      IF (PorosityValues(PorosityPerm(CurrentNode)) .NE. PorosityValues(PorosityPerm(CurrentNode))) THEN
        PRINT *,"PermafrostPorosityEvolution: ","Found weird number for PorosityValues"
        PRINT *,"PermafrostPorosityEvolution: PorosityValues(",PorosityPerm(CurrentNode),")=",&
             PorosityValues(PorosityPerm(CurrentNode))
        PRINT *,"PermafrostPorosityEvolution:  Prev=",PorosityVariable % PrevValues(PorosityPerm(CurrentNode),1)
        PRINT *,"PermafrostPorosityEvolution: ",PrevNodalrhos,Nodalrhos,StrainInvariant
        CALL FATAL(SolverName,'Exiting')
      END IF
    END DO
  END DO
  IF (ConstantTemp) &
       DEALLOCATE(PrevNodalTemperature)
  FirstVisit = .FALSE.
END SUBROUTINE PermafrostPorosityEvolution
!-----------------------------------------------------------------------------
!  Compute invariant of stress tensor and its time-derivative
!> \ingroup Solvers
!-----------------------------------------------------------------------------
SUBROUTINE PermafrostStressInvariant( Model,Solver,dt,TransientSimulation )
  !-----------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostStressInvariant'
  CHARACTER(LEN=MAX_NAME_LEN) :: VariableName, StressVariableName, PressureName
  LOGICAL :: Found, FirstTime=.TRUE., NoPressure, UpdatePrev=.FALSE.,SteadyState=.FALSE.
  TYPE(Variable_t), POINTER :: InvariantVar, InvariantVeloVar, StressVariableVar, PressureVar
  INTEGER, POINTER :: InvariantPerm(:), InvariantVeloPerm(:), StressVariablePerm(:), PressurePerm(:)
  INTEGER :: I, DIM, StressVariableDOFs, CurrentTime, activenodes
  REAL (KIND=dp), POINTER :: Invariant(:), InvariantVelo(:), StressVariable(:),&
       InvariantPrev(:,:), Pressure(:), PrevPressure(:)
  REAL (KIND=dp) :: AverageInvariant
  !------------------------------------------------------------------------------
  SAVE FirstTime, CurrentTime, DIM, NoPressure,&
       PressureVar, Pressure, PressurePerm, PressureName,&
       StressVariableVar, StressVariable, StressVariablePerm, StressVariableDOFs, &
       StressVariableName, InvariantVeloVar, InvariantVeloPerm, InvariantVelo,SteadyState
 
  
  CALL INFO( SolverName, '-----------------------------------------',Level=4 )
  CALL INFO( SolverName, ' Computing Permafrost Stress Invariant',Level=4 )
  CALL INFO( SolverName, '-----------------------------------------',Level=4 )

  SolverParams => GetSolverParams()


  IF (FirstTime) THEN
    CurrentTime=GetTimeStep()
    UpdatePrev = .FALSE.
    DIM=CoordinateSystemDimension()
  ELSE
    IF (CurrentTime .NE. GetTimeStep()) THEN
      UpdatePrev = .TRUE.
      CurrentTime=GetTimeStep()
    ELSE
      UpdatePrev = .FALSE.
    END IF
  END IF

  IF (FirstTime .OR. Model % Mesh % Changed) THEN
    CALL INFO( SolverName, 'Initialization step:',Level=9 )
    SteadyState = GetLogical(SolverParams,'Steady State',Found)
    IF (SteadyState) THEN
      CALL Info (SolverName,'Computing steady state only',Level=4)
      UpdatePrev = .FALSE.
    END IF
    PressureName = ListGetString(SolverParams, &
         'Pressure Variable', Found )
    IF (.NOT.Found) THEN
      CALL WARN(SolverName," 'Pressure Variable' not found. Using default 'Pressure' ")
      WRITE(PressureName,'(A)') 'Pressure'
    ELSE
      WRITE(Message,'(A,A)') "'Pressure Variable' found and set to: ", PressureName
      CALL INFO(SolverName,Message,Level=9)
    END IF
    PressureVar => VariableGet(Solver % Mesh % Variables,PressureName)
    IF (.NOT.ASSOCIATED(PressureVar)) THEN
      NULLIFY(Pressure)
      WRITE(Message,'(A,A,A)') 'Pressure Variable "', TRIM(PressureName), '" not associated'
      CALL FATAL(SolverName,Message)
    ELSE
      Pressure => PressureVar % Values
      IF (ASSOCIATED(Pressure)) THEN
        PressurePerm => PressureVar % Perm
        PrevPressure => PressureVar % PrevValues(:,1)
        NoPressure = .FALSE.
        WRITE(Message,'(A,A,A)') 'Pressure Variable "', TRIM(PressureName), '" associated'
        CALL INFO(SolverName,Message,Level=9)
      ELSE
        CALL FATAL(SolverName,'Pressure values not associated')
      END IF
    END IF
    CALL INFO( SolverName, 'Initialization completed',Level=9 )
  END IF
  
  StressVariableName = ListGetString(SolverParams,'Stress Variable Name',Found)
  IF (.NOT.Found) CALL FATAL(SolverName,' "Stress Variable Name" not found')
  StressVariableVar => VariableGet( Solver % Mesh % Variables, StressVariableName )
  IF ( ASSOCIATED(StressVariableVar)) THEN
    StressVariablePerm => StressVariableVar % Perm
    StressVariable => StressVariableVar % Values
    StressVariableDOFs = StressVariableVar % DOFs
  ELSE
    CALL FATAL(SolverName, TRIM(StressVariableName)//' not found')
  END IF
    
  InvariantVar =>  Solver % Variable
  IF ( ASSOCIATED(InvariantVar)) THEN    
    VariableName = InvariantVar % Name
    Invariant => Solver % Variable % Values
    InvariantPerm => Solver % Variable % Perm
    InvariantPrev => Solver % Variable % PrevValues
    WRITE (Message,*) 'Solver variable ', TRIM(VariableName),', found'
    CALL INFO(SolverName,Message,Level=9)
  ELSE
    CALL FATAL(SolverName,'Solver variable not associated')
  END IF

  InvariantVeloVar => VariableGet( Solver % Mesh % Variables, TRIM(VariableName) // " Velocity" )
  IF ( ASSOCIATED(InvariantVeloVar)) THEN
    InvariantVeloPerm => InvariantVeloVar % Perm
    InvariantVelo =>InvariantVeloVar  % Values
     CALL INFO(SolverName, TRIM(VariableName) // " Velocity"//' found',Level=9)
  ELSE
    CALL FATAL(SolverName, TRIM(VariableName) // " Velocity"//' not found')
  END IF
  !!!!!!!!!!!!!
  !!!! VariableAdd for TRIM(VariableName) // " Velocity"
  
  activenodes = 0
  InvariantVelo = 0.0_dp  
  DO I = 1,Solver % Mesh % Nodes % NumberOfNodes
    IF (InvariantPerm(I) == 0) CYCLE
    IF (PressurePerm(I) == 0) THEN
      WRITE (Message,*) 'No entry for pressure variable',PressureName,'at point',I
      CALL FATAL(SolverName,Message)
    END IF
    IF ( UpdatePrev .AND. .NOT.(SteadyState) ) THEN
      InvariantPrev(InvariantPerm(I),1) =  Invariant(InvariantPerm(I))  
    END IF
    !PRINT *,"PermafrostStressInvariant:",Pressure(PressurePerm(I))
    !PRINT *,InvariantPerm(I)
    !PRINT *,Invariant(InvariantPerm(I))
    !PRINT *, (StressVariablePerm(I)-1)*StressVariableDOFs
    Invariant(InvariantPerm(I)) = &
         GetFirstInvariant(StressVariable,Pressure(PressurePerm(I)),(StressVariablePerm(I)-1)*StressVariableDOFs,DIM)
    IF (SteadyState) THEN
      InvariantVelo(InvariantVeloPerm(I)) = 0.0
    ELSE
      InvariantVelo(InvariantVeloPerm(I)) = &
           (Invariant(InvariantPerm(I))-InvariantPrev(InvariantPerm(I),1))/dt
      !PRINT *,"PermafrostStressInvariant:",InvariantVelo(InvariantVeloPerm(I))
    END IF
    !IF (InvariantVelo(InvariantVeloPerm(I)) /= 0.0_dp) THEN
    !  PRINT *, "InvariantVelo:", InvariantVelo(InvariantVeloPerm(I)), Invariant(InvariantPerm(I)), InvariantPrev(InvariantPerm(I),1)
    !END IF
    activenodes = activenodes + 1
    AverageInvariant = AverageInvariant + Invariant(InvariantPerm(I))
    IF (FirstTime .AND. .NOT.(SteadyState)) THEN
      InvariantPrev(InvariantPerm(I),1) = Invariant(InvariantPerm(I))
    END IF
  END DO
  AverageInvariant = AverageInvariant/DBLE(activenodes)
  WRITE(Message,'(A,I0,A,I0,A,ES12.3)') 'Average invariant of ', activenodes,&
      ' out of ', Solver % Mesh % Nodes % NumberOfNodes,&
      ' active nodes:', AverageInvariant
  CALL INFO(SolverName,Message,Level=5)
  FirstTime = .FALSE.
  CONTAINS
    FUNCTION GetFirstInvariant(Stress,PressureAtPoint,Position,DIM) RESULT(FirstInvariant)
      REAL (KIND=dp) ::  FirstInvariant
      REAL (KIND=dp), POINTER :: Stress(:)
      REAL (KIND=dp) :: PressureAtPoint
      INTEGER :: DIM, Position
      !--------
      INTEGER :: I
      
      FirstInvariant = 0.0_dp
      ! tr(sigma - 1p)
      DO I=1,DIM
        FirstInvariant = FirstInvariant + Stress(Position+I) - PressureAtPoint
      END DO
      
    END FUNCTION GetFirstInvariant  
END SUBROUTINE PermafrostStressInvariant
!---------------------------------------------------------------------------------------------
! Functions needed for permafrost model (might be shifted to USF-directory)
!---------------------------------------------------------------------------------------------
FUNCTION GetKGuu(Model,IPNo,PorosityAtIP) RESULT(KGuuAtIP)
  USE DefUtils
  USE PermaFrostMaterials
  IMPLICIT NONE
  
  TYPE(Model_t) :: Model
  INTEGER, INTENT(IN) :: IPNo
  REAL(KIND=dp) :: PorosityAtIP, KGuuAtIP(6,6)
  !-----------
  TYPE(Solver_t) :: DummySolver
  TYPE(ValueList_t), POINTER :: Material
  TYPE(Element_t),POINTER :: Element
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: RockMaterialID, NumberOfRockRecords, DIM, t, i, IPPerm
  REAL(KIND=dp) :: EGAtIP, nuGAtIP
  TYPE(Variable_t), POINTER :: XiAtIPVar
  INTEGER, POINTER :: XiAtIPPerm(:)
  REAL(KIND=dp), POINTER :: XiAtIP(:)
  LOGICAL :: FirstTime = .TRUE., ElementWiseRockMaterial, Found
  CHARACTER(LEN=MAX_NAME_LEN) :: ElementRockMaterialName
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName = 'PermafrostMaterials (GetKGuu)'
  !-----------
  SAVE FirstTime,NumberOfRockRecords,CurrentSolventMaterial,DIM,ElementWiseRockMaterial

  IF (IPNo >= 0) THEN
    CALL FATAL(FunctionName,'Invalid IP number in argument')
  END IF
  Element => Model % CurrentElement
  IF (.NOT.ASSOCIATED(Element)) CALL FATAL(FunctionName,'Element not associated')
  t = Element % ElementIndex
  Material => GetMaterial(Element)

  XiAtIPVar => VariableGet( Model % Mesh % Variables, 'Xi')
  IF (.NOT.ASSOCIATED(XiAtIPVar)) THEN
    CALL FATAL(FunctionName,'Variable Xi is not associated')
  END IF
  XiAtIPPerm => XiAtIPVar % Perm
  XiAtIp => XiAtIPVar % Values
  
  IPPerm = XiAtIPPerm(t) + ABS(IPNo)

  IF (FirstTime .OR. (Model % Mesh % Changed)) THEN
    DIM =  CoordinateSystemDimension()

    ! check, whether we have globally or element-wise defined values of rock-material parameters
    ElementRockMaterialName = GetString(Material,'Element Rock Material File',ElementWiseRockMaterial)
    IF (ElementWiseRockMaterial) THEN
      CALL INFO(FunctionName,'Found "Element Rock Material File"',Level=5)
      CALL INFO(FunctionName,'Using element-wise rock material definition',Level=5)
    END IF
    IF (ElementWiseRockMaterial) THEN
      ! read element-wise material parameter (GlobalRockMaterial will have one entry each element)
      NumberOfRockRecords = &
           ReadPermafrostElementRockMaterial(ElementRockMaterialName,DummySolver,DIM,SkipInit=.TRUE.)
    ELSE
      NumberOfRockRecords =  ReadPermafrostRockMaterial( Material )
    END IF

    IF (NumberOfRockRecords < 1) THEN
      CALL FATAL(FunctionName,'No Rock Material specified')
    ELSE
      CALL INFO(FunctionName,'Permafrost Rock Material read',Level=6)
      FirstTime = .FALSE.
    END IF
    CALL SetPermafrostSolventMaterial( CurrentSolventMaterial )
  END IF

  IF (ElementWiseRockMaterial) THEN
    RockMaterialID = t  ! each element has it's own set of parameters
  ELSE
    RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found, UnfoundFatal=.TRUE.)
  END IF

  EGAtIP = EG(CurrentSolventMaterial,RockMaterialID,XiAtIP(IPPerm),PorosityAtIP)
  nuGAtIP = nuG(CurrentSolventMaterial,RockMaterialID,XiAtIP(IPPerm),PorosityAtIP)
  KGuuAtIP = KGuu(EGAtIP,nuGAtIP,DIM)
END FUNCTION GetKGuu
!---------------------------------------------------------------------------------------------
FUNCTION GetBetaG(Model,DummyIPNo,ArgumentsAtIP) RESULT(betaGAtIP)
  USE DefUtils
  USE PermaFrostMaterials
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER, INTENT(IN) :: DummyIPNo
  REAL(KIND=dp) :: ArgumentsAtIP(2), betaGAtIP
  !--------------
  TYPE(Solver_t) :: DummySolver
  TYPE(ValueList_t), POINTER :: Material
  TYPE(Element_t),POINTER :: Element
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: RockMaterialID, NumberOfRockRecords, DIM, t, i, IPPerm
  REAL(KIND=dp) :: EGAtIP, nuGAtIP, PorosityAtIP, XiAtIP
  LOGICAL :: FirstTime = .TRUE., ElementWiseRockMaterial, Found
  CHARACTER(LEN=MAX_NAME_LEN) :: ElementRockMaterialName
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName = 'PermafrostMaterials (GetBetaG)'
  !-----------
  SAVE FirstTime,NumberOfRockRecords,DIM,ElementWiseRockMaterial

  Element => Model % CurrentElement
  IF (.NOT.ASSOCIATED(Element)) CALL FATAL(FunctionName,'Element not associated')
  t = Element % ElementIndex
  Material => GetMaterial(Element)

  IF (FirstTime .OR. (Model % Mesh % Changed)) THEN
    DIM =  CoordinateSystemDimension()

    ! check, whether we have globally or element-wise defined values of rock-material parameters
    ElementRockMaterialName = GetString(Material,'Element Rock Material File',ElementWiseRockMaterial)
    IF (ElementWiseRockMaterial) THEN
      CALL INFO(FunctionName,'Found "Element Rock Material File"',Level=5)
      CALL INFO(FunctionName,'Using element-wise rock material definition',Level=5)
    END IF
    IF (ElementWiseRockMaterial) THEN
      ! read element-wise material parameter (GlobalRockMaterial will have one entry each element)
      NumberOfRockRecords = &
           ReadPermafrostElementRockMaterial(ElementRockMaterialName,DummySolver,DIM,SkipInit=.TRUE.)
    ELSE
      NumberOfRockRecords =  ReadPermafrostRockMaterial( Material )
    END IF

    IF (NumberOfRockRecords < 1) THEN
      CALL FATAL(FunctionName,'No Rock Material specified')
    ELSE
      CALL INFO(FunctionName,'Permafrost Rock Material read',Level=6)
      FirstTime = .FALSE.
    END IF
    CALL SetPermafrostSolventMaterial( CurrentSolventMaterial )
  END IF

  IF (ElementWiseRockMaterial) THEN
    RockMaterialID = t  ! each element has it's own set of parameters
  ELSE
    RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found, UnfoundFatal=.TRUE.)
  END IF
  PorosityAtIP = ArgumentsAtIP(1)
  XiAtIP = ArgumentsAtIP(2)
  betaGAtIP = betaG(CurrentSolventMaterial,RockMaterialID,XiAtIP,PorosityAtIP)
END FUNCTION GetBetaG
  !---------------------------------------------------------------------------------------------
FUNCTION GetNuG(Model,DummyIPNo,ArgumentsAtIP) RESULT(nuGAtIP)
  USE DefUtils
  USE PermaFrostMaterials
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER, INTENT(IN) :: DummyIPNo
  REAL(KIND=dp) :: ArgumentsAtIP(2), nuGAtIP
  
  !-----
  TYPE(Solver_t) :: DummySolver
  TYPE(ValueList_t), POINTER :: Material
  TYPE(Element_t),POINTER :: Element
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: RockMaterialID, NumberOfRockRecords, DIM, t, i, IPPerm
  REAL(KIND=dp) :: PorosityAtIP, XiAtIP
  LOGICAL :: Found, FirstTime = .TRUE., ElementWiseRockMaterial
  CHARACTER(LEN=MAX_NAME_LEN) :: ElementRockMaterialName
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName = 'PermafrostMaterials (GetNuG)'
  !-----------
  SAVE FirstTime,NumberOfRockRecords,CurrentSolventMaterial,DIM,ElementWiseRockMaterial
  
  IF (FirstTime) CALL INFO("Permafrost(GetNuG)","Initializing",Level=4)
  PorosityAtIP=ArgumentsAtIP(1)
  XiAtIP=ArgumentsAtIP(2)
  !XiAtIP=1.0_dp 
  Element => Model % CurrentElement
  IF (.NOT.ASSOCIATED(Element)) CALL FATAL(FunctionName,'Element not associated')
  t = Element % ElementIndex
  Material => GetMaterial(Element)

  IF (FirstTime .OR. (Model % Mesh % Changed)) THEN
    DIM =  CoordinateSystemDimension()

    ! check, whether we have globally or element-wise defined values of rock-material parameters
    ElementRockMaterialName = GetString(Material,'Element Rock Material File',ElementWiseRockMaterial)
    IF (ElementWiseRockMaterial) THEN
      CALL INFO(FunctionName,'Found "Element Rock Material File"',Level=5)
      CALL INFO(FunctionName,'Using element-wise rock material definition',Level=5)
    END IF
    IF (ElementWiseRockMaterial) THEN
      ! read element-wise material parameter (GlobalRockMaterial will have one entry each element)
      NumberOfRockRecords = &
           ReadPermafrostElementRockMaterial(ElementRockMaterialName,DummySolver,DIM,SkipInit=.TRUE.)
    ELSE
      NumberOfRockRecords =  ReadPermafrostRockMaterial( Material )
    END IF

    IF (NumberOfRockRecords < 1) THEN
      CALL FATAL(FunctionName,'No Rock Material specified')
    ELSE
      CALL INFO(FunctionName,'Permafrost Rock Material read',Level=6)
      FirstTime = .FALSE.
    END IF
    CALL SetPermafrostSolventMaterial( CurrentSolventMaterial )
  END IF

  IF (ElementWiseRockMaterial) THEN
    RockMaterialID = t  ! each element has it's own set of parameters
  ELSE
    RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
  END IF

  nuGAtIP = nuG(CurrentSolventMaterial,RockMaterialID,XiAtIP,PorosityAtIP)
  !PRINT *,"getNuG:", nuGAtIP, XiAtIp(IPPerm),PorosityAtIP
END FUNCTION GetNuG
!---------------------------------------------------------------------------------------------
FUNCTION GetEG(Model,DummyIPNo,ArgumentsAtIP) RESULT(EGAtIP)
  USE DefUtils
  USE PermaFrostMaterials
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER, INTENT(IN) :: DummyIPNo
  REAL(KIND=dp) :: ArgumentsAtIP(2), EGAtIP
  !-----
  TYPE(Solver_t) :: DummySolver
  TYPE(ValueList_t), POINTER :: Material
  TYPE(Element_t),POINTER :: Element
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: RockMaterialID, NumberOfRockRecords, DIM, t, i, IPPerm
  REAL(KIND=dp) :: PorosityAtIP, XiAtIP
  LOGICAL :: Found,FirstTime = .TRUE., ElementWiseRockMaterial
  CHARACTER(LEN=MAX_NAME_LEN) :: ElementRockMaterialName
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName = 'PermafrostMaterials (GetNuG)'
  !-----------
  SAVE FirstTime,NumberOfRockRecords,CurrentSolventMaterial,DIM,ElementWiseRockMaterial

  PorosityAtIP=ArgumentsAtIP(1)
  XiAtIP=ArgumentsAtIP(2)
  !PRINT *, "GetEG:", PorosityAtIP, XiAtIP
  
  Element => Model % CurrentElement
  IF (.NOT.ASSOCIATED(Element)) CALL FATAL(FunctionName,'Element not associated')
  !t = Element % ElementIndex
  Material => GetMaterial(Element)
  
  IF (FirstTime .OR. (Model % Mesh % Changed)) THEN
    DIM =  CoordinateSystemDimension()

    ! check, whether we have globally or element-wise defined values of rock-material parameters
    ElementRockMaterialName = GetString(Material,'Element Rock Material File',ElementWiseRockMaterial)
    IF (ElementWiseRockMaterial) THEN
      WRITE (Message,*) 'Found "Element Rock Material File"'
      CALL INFO(FunctionName,Message,Level=5)
      CALL INFO(FunctionName,'Using element-wise rock material definition',Level=5)
    END IF
    IF (ElementWiseRockMaterial) THEN
      ! read element-wise material parameter (GlobalRockMaterial will have one entry each element)
      NumberOfRockRecords = &
           ReadPermafrostElementRockMaterial(ElementRockMaterialName,DummySolver,DIM,SkipInit=.TRUE.)
    ELSE
      NumberOfRockRecords =  ReadPermafrostRockMaterial( Material )
    END IF

    IF (NumberOfRockRecords < 1) THEN
      CALL FATAL(FunctionName,'No Rock Material specified')
    ELSE
      CALL INFO(FunctionName,'Permafrost Rock Material read',Level=6)
      FirstTime = .FALSE.
    END IF
    CALL SetPermafrostSolventMaterial( CurrentSolventMaterial )
  END IF

  IF (ElementWiseRockMaterial) THEN
    RockMaterialID = t  ! each element has it's own set of parameters
  ELSE
    RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
  END IF
  EGAtIP = EG(CurrentSolventMaterial,RockMaterialID,XiAtIP,PorosityAtIP)
  !IF ((EGAtIP < 1.0d05) .OR. (EGAtIP > 1.0d07)) PRINT *,"GetEG",EGAtIP,XiAtIP,PorosityAtIP
END FUNCTION GetEG
!---------------------------------------------------------------------------------------------
FUNCTION GetElasticityForce(Model,IPNo,ArgumentsAtIP) RESULT(EforceAtIP) ! needs arguments Temperature, Pressure, Porosity, Salinity
  USE DefUtils
  USE PermaFrostMaterials
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER, INTENT(IN) :: IPNo
  REAL(KIND=dp) :: ArgumentsAtIP(5), EforceAtIP
  !--------------
  REAL(KIND=dp) :: TemperatureAtIP, PressureAtIP, PorosityAtIP, SalinityAtIP,XiAtIP,&
       rhogwAtIP, rhosAtIP, rhowAtIP,rhocAtIP, rhoiAtIP,rhoGAtIP,&
       GasConstant, N0, DeltaT, T0, p0, eps, Gravity(3),InitialOffsetRhoAtIP
  TYPE(Variable_t), POINTER :: RhoOffsetAtIPVar
  INTEGER, POINTER :: RhoOffsetAtIPPerm(:)
  REAL(KIND=dp), POINTER :: RhoOffsetAtIP(:)
  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: Material
  INTEGER ::  DIM, t,NumberOfRockRecords, RockMaterialID, IPPerm
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
  TYPE(Solver_t) :: DummySolver
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName = 'Permafrost (GetElasticityForce)'
  CHARACTER(LEN=MAX_NAME_LEN) :: ElementRockMaterialName
  LOGICAL :: Found,FirstTime=.TRUE.,ConstVal=.FALSE., ConstantsRead = .FALSE.,&
       ElementWiseRockMaterial, OffsetDensity=.FALSE.

  SAVE ConstantsRead,ElementWiseRockMaterial,GasConstant, DIM, N0, DeltaT, T0, p0, eps, Gravity,&
       NumberOfRockRecords,FirstTime,CurrentSoluteMaterial,CurrentSolventMaterial,&
       OffsetDensity

  IF (IPNo >= 0) THEN
    WRITE(Message,*) 'IP number invalid:', IPNo, ' - is this an IP valid call to this routine?'
    CALL FATAL(FunctionName,Message)
  END IF
  IF (.NOT.ConstantsRead) THEN
    ConstantsRead = &
         ReadPermafrostConstants(Model, FunctionName, DIM, GasConstant, N0, DeltaT, T0, p0, eps, Gravity)
    OffsetDensity = GetLogical(Model %Constants,'Permafrost Offset Density', Found)
    IF (.NOT.Found) THEN
      OffsetDensity = .FALSE.
      CALL WARN(FunctionName,'No offset for groundwater pressure included - might lead to artifial high compression')
    ELSE
      CALL INFO(FunctionName,'Offset groundwater pressure is activated',Level=5)
    END IF
  END IF
  Element => Model % CurrentElement
  IF (.NOT.ASSOCIATED(Element)) CALL FATAL(FunctionName,'Element not asssociated')
  t = Element % ElementIndex
  IF (t <= 0) THEN
    WRITE(Message,*) 'Element number invalid:', t
    CALL FATAL(FunctionName,Message)
  END IF
  IF (OffsetDensity) THEN
    RhoOffsetAtIPVar => VariableGet( Model % Mesh % Variables, 'Reference Offset Density')
    IF (.NOT.ASSOCIATED(RhoOffsetAtIPVar)) THEN
      WRITE(Message,*) '"Permafrost Offset Density" is set, but variable "Reference Offset Density" is not associated'
      CALL FATAL(FunctionName,Message)
    END IF
    RhoOffsetAtIPPerm => RhoOffsetAtIPVar % Perm
    IPPerm = RhoOffsetAtIPPerm(t) + ABS(IPNo)
    RhoOffsetAtIP => RhoOffsetAtIPVar % Values
    InitialOffsetRhoAtIP = RhoOffsetAtIP(IPPerm)
  ELSE
    InitialOffsetRhoAtIP = 0.0_dp
  END IF

  Material => GetMaterial(Element)
  ConstVal = GetLogical(Material,'Constant Permafrost Properties',Found)
  IF (FirstTime) THEN
    ! check, whether we have globally or element-wise defined values of rock-material parameters
    ElementRockMaterialName = GetString(Material,'Element Rock Material File',ElementWiseRockMaterial)
    IF (ElementWiseRockMaterial) THEN
      CALL INFO(FunctionName,'Found "Element Rock Material File"',Level=5)
      CALL INFO(FunctionName,'Using element-wise rock material definition',Level=5)
    END IF
    IF (ElementWiseRockMaterial) THEN
      ! read element-wise material parameter (GlobalRockMaterial will have one entry each element)
      NumberOfRockRecords = &
           ReadPermafrostElementRockMaterial(ElementRockMaterialName,DummySolver,DIM,SkipInit=.TRUE.)
    ELSE
      NumberOfRockRecords =  ReadPermafrostRockMaterial( Material )
    END IF
    IF (NumberOfRockRecords < 1) THEN
      PRINT *, "NumberOfRockRecords=", NumberOfRockRecords
      CALL FATAL(FunctionName,'No Rock Material specified')
    ELSE
      CALL INFO(FunctionName,'Permafrost Rock Material read',Level=5)
      FirstTime = .FALSE.
    END IF
    CALL ReadPermafrostSoluteMaterial( Material,Model % Constants,CurrentSoluteMaterial )
    CALL SetPermafrostSolventMaterial( CurrentSolventMaterial )
  END IF

  IF (ElementWiseRockMaterial) THEN
    RockMaterialID = t
  ELSE      
    RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found)
    IF (.NOT.Found) CALL FATAL(FunctionName,"Rock Material ID not found")
  END IF

  TemperatureAtIP = ArgumentsAtIP(1)
  PressureAtIP    = ArgumentsAtIP(2)
  PorosityAtIP    = ArgumentsAtIP(3)
  SalinityAtIP    = ArgumentsAtIP(4)
  XiAtIP          = ArgumentsAtIP(5)
  
  rhocAtIP =  rhoc(CurrentSoluteMaterial,T0,p0,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal)
  !IF (rhocAtIP .NE. rhocAtIP) CALL FATAL(FunctionName,'rhocAtIP is NaN')  
  rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)
  !IF (rhowAtIP .NE. rhocAtIP) CALL FATAL(FunctionName,'rhowAtIP is NaN') 
  rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)
  rhosAtIP = rhos(RockMaterialID,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)
  rhogwAtIP = rhogw(rhowAtIP,rhocAtIP,XiAtIP,SalinityAtIP)
  !!!! IF (Time == 0) store that one at IP

  
  rhoGAtIP = rhoG(rhosAtIP,rhogwAtIP,rhoiAtIP,PorosityAtIP,SalinityAtIP,XiAtIP)
  IF (rhoGAtIP .NE. rhoGAtIP) THEN
    PRINT *,rhosAtIP,rhogwAtIP,rhoiAtIP,PorosityAtIP,SalinityAtIP,XiAtIP
    CALL FATAL(FunctionName,'rhoGAtIP is NaN')
  END IF
  
  !EforceAtIP = -rhoGAtIP * SQRT(SUM(Gravity(1:3)*Gravity(1:3)))
  EforceAtIP = (rhoGAtIP - InitialOffsetRhoAtIP)* Gravity(DIM)
  !IF (t==100 )PRINT *, "GetElasticityForce: EforceAtIP=", EforceAtIP, "=(",rhoGAtIP, "-", InitialOffsetRhoAtIP,")*", Gravity(DIM)
  !     PorosityAtIP,SalinityAtIP,XiAtIP, "   *********"  
END FUNCTION GetElasticityForce
