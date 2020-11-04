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
! *  Authors: Samuel Cook 
! *  Email:   sc690@cam.ac.uk
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - Scientific Computing Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 18 January 2018
! *
! *****************************************************************************/
!> Interpolates ice variables necessary for hydrology from ice mesh to hydro
!> mesh as part of combined calving-hydrology simulations.
   SUBROUTINE IceToHydroInterp( Model,Solver,Timestep,TransientSimulation )
!******************************************************************************
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh,materials,BCs,etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear equation solver options
!
!  REAL(KIND=dp) :: Timestep
!     INPUT: Timestep size for time dependent simulations
!
!******************************************************************************
     USE Differentials
     USE MaterialModels
     USE DefUtils
     USE InterpVarToVar
!------------------------------------------------------------------------------
     IMPLICIT NONE
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!    External variables
!------------------------------------------------------------------------------
     TYPE(Model_t)  :: Model
     TYPE(Solver_t), TARGET :: Solver
     LOGICAL :: TransientSimulation
     REAL(KIND=dp) :: Timestep
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Mesh_t), POINTER :: IceMesh, HydroMesh
     TYPE(Solver_t), POINTER ::  HydroSolver, TempSolver
     TYPE(Variable_t), POINTER :: HydroVar=>NULL(),WorkVar=>NULL(),&
                                  InterpVar1=>NULL(),InterpVar2=>NULL(),&
                                  InterpVar3=>NULL(),InterpVar4=>NULL(),&
                                  InterpVar5=>NULL(),WorkVar2=>NULL()
     TYPE(Variable_t), TARGET :: InterpVar1Copy, InterpVar2Copy,&
                                 InterpVar3Copy, InterpVar4Copy, InterpVar5Copy
     TYPE(ValueList_t), POINTER :: Params, BC
     TYPE(Element_t), POINTER :: Element=>NULL()
     CHARACTER(MAX_NAME_LEN) :: Variable1, Variable2, Variable3, Variable4,&
                                Variable5, VarName, iString
     LOGICAL :: Found=.FALSE., FirstTime=.TRUE., LoadReaders, Hit=.FALSE.
     LOGICAL, POINTER :: BasalLogical(:)
     INTEGER, POINTER :: InterpDim(:), IceMeshBasePerm(:)=>NULL(),&
                         NPP(:)=>NULL(), RefNode(:)=>NULL(), NewPerm1(:),&
                         ReaderPerm1(:), NewPerm2(:), NewPerm3(:),&
                         NewPerm4(:), NewPerm5(:), ReaderPerm2(:),&
                         ReaderPerm3(:), ReaderPerm4(:), ReaderPerm5(:),&
                         ReaderPerm6(:), ReaderPerm7(:), ReaderPerm8(:),&
                         ReaderPerm9(:), ReaderPerm10(:)
     INTEGER, ALLOCATABLE, TARGET :: DefaultPerm(:), ZeroNodes(:)
     INTEGER :: i, j, k, HPSolver, Reader1, Reader2, Reader3,&
                Reader4, Reader5, Reader, NumVar, DummyInt, ierr,&
                ElementBC, BasalBC, n, NearestNode, SideBC1, SideBC2, HitCount,&
                ZeroCounter, TSolver
     REAL(KIND=dp), POINTER :: NVP(:)=>NULL(), NewValues1(:), NewValues2(:),&
                               NewValues3(:), NewValues4(:), NewValues5(:),&
                               ReaderValues1(:), ReaderValues2(:),&
                               ReaderValues3(:), ReaderValues4(:),&
                               ReaderValues5(:), ReaderValues6(:),&
                               ReaderValues7(:), ReaderValues8(:),&
                               ReaderValues9(:), ReaderValues10(:)  
     REAL(KIND=dp) :: IceTempResSum, HydroTempResSum, ScaleFactor,&
                      ElementTempResSum, Dist, Threshold, MinDist, x, y,&
                      NewNS, MeanNS
     REAL(KIND=dp), ALLOCATABLE :: NSValues(:)
     REAL(KIND=dp), ALLOCATABLE, TARGET :: ParITRS(:), ParHTRS(:), NewValues6(:)
     SAVE HydroMesh, FirstTime, HPSolver, TSolver
          !NewPerm1,&
          !NewValues1,NewValues2, NewValues3, NewValues4, NewValues5, DefaultPerm
!------------------------------------------------------------------------------
    Params => GetSolverParams()

    !For loading variables that remain constant and are read in from files
    !E.g. Zb or surface melting
    !Various keywords need to be specified in solver section of SIF
    IF(FirstTime) THEN
      DO i=1,Model % NumberOfSolvers
        IF(Model % Solvers(i) % Variable % Name == 'hydraulic potential') THEN
          HPSolver = i
          HydroSolver => Model % Solvers(HPSolver)
          EXIT
        END IF
      END DO
      DO i=1,Model % NumberOfSolvers
        IF(Model % Solvers(i) % Variable % Name == 'temp') THEN
          TSolver = i
          TempSolver => Model % Solvers(TSolver)
          EXIT
        END IF
      END DO
      
      LoadReaders = GetLogical(Params,'Load Reader Variables',Found)
      IF(.NOT. Found) LoadReaders = .FALSE.
      IF(LoadReaders) THEN
        NumVar = GetInteger(Params, 'Number Of Variables To Read',Found)
        IF(.NOT. Found) CALL Fatal('Ice2Hydro', 'Number of variables to read &
        not specified')
        IF(NumVar > 10) CALL Info('Ice2Hydro', 'Not set up for more than 10&
        reader variables - increase maximum limit in source code', Level=1)

        !A default perm in case the variables to be loaded don't have a perm
        !associated. Assumes node 1 = value 1
        ALLOCATE(DefaultPerm(HydroSolver % Mesh % NumberOfNodes))
        DO i=1, HydroSolver % Mesh % NumberOfNodes
          DefaultPerm(i) = i                
        END DO

        DO i=1, NumVar
          iString = STR(i)
          Reader = GetInteger(Params, 'Reader Solver '//iString, Found)
          IF(Found) THEN
            VarName = GetString(Params, 'Reader V'//iString, Found)
            IF(.NOT. Found) PRINT *, 'No reader '//iString//' variable specified'
            !For the case if you're restarting from a mesh where this variable
            !has already been added
            WorkVar => VariableGet(HydroSolver % Mesh % Variables,&
                       VarName, ThisOnly=.TRUE.)
            IF(ASSOCIATED(WorkVar)) CYCLE

            WorkVar => VariableGet(Model % Solvers(Reader) % Mesh % Variables,&
                       VarName, ThisOnly=.TRUE.)

            !Allocate the variable perm to default value if not associated
            IF(.NOT. ASSOCIATED(WorkVar % Perm)) THEN
              ALLOCATE(WorkVar % Perm(HydroSolver % Mesh % NumberOfNodes))
              WorkVar % Perm(1:SIZE(WorkVar % Perm)) = DefaultPerm(1:SIZE(WorkVar % Perm))
            END IF
           
            !There has to be a better way to do this, but you can't use indices
            !in variable names, so I can't see a way of allocating the arrays in
            !the loop, such that each variable gets its own perm and values
            IF(i==1) ALLOCATE(ReaderPerm1(SIZE(WorkVar % Perm)), ReaderValues1(SIZE(WorkVar % Values)))
            IF(i==2) ALLOCATE(ReaderPerm2(SIZE(WorkVar % Perm)), ReaderValues2(SIZE(WorkVar % Values)))
            IF(i==3) ALLOCATE(ReaderPerm3(SIZE(WorkVar % Perm)), ReaderValues3(SIZE(WorkVar % Values)))
            IF(i==4) ALLOCATE(ReaderPerm4(SIZE(WorkVar % Perm)), ReaderValues4(SIZE(WorkVar % Values)))
            IF(i==5) ALLOCATE(ReaderPerm5(SIZE(WorkVar % Perm)), ReaderValues5(SIZE(WorkVar % Values)))
            IF(i==6) ALLOCATE(ReaderPerm6(SIZE(WorkVar % Perm)), ReaderValues6(SIZE(WorkVar % Values)))
            IF(i==7) ALLOCATE(ReaderPerm7(SIZE(WorkVar % Perm)), ReaderValues7(SIZE(WorkVar % Values)))
            IF(i==8) ALLOCATE(ReaderPerm8(SIZE(WorkVar % Perm)), ReaderValues8(SIZE(WorkVar % Values)))
            IF(i==9) ALLOCATE(ReaderPerm9(SIZE(WorkVar % Perm)), ReaderValues9(SIZE(WorkVar % Values)))
            IF(i==10) ALLOCATE(ReaderPerm10(SIZE(WorkVar % Perm)), ReaderValues10(SIZE(WorkVar % Values)))
            IF(i==1) THEN
              ReaderPerm1 = WorkVar % Perm
              ReaderValues1 = WorkVar % Values
              CALL VariableAdd(HydroSolver % Mesh % Variables, HydroSolver % & 
                Mesh, Model % Solvers(HPSolver), VarName, 1, ReaderValues1, &
                ReaderPerm1)
            ELSEIF(i==2) THEN
              ReaderPerm2 = WorkVar % Perm
              ReaderValues2 = WorkVar % Values
              CALL VariableAdd(HydroSolver % Mesh % Variables, HydroSolver % & 
                Mesh, Model % Solvers(HPSolver), VarName, 1, ReaderValues2, &
                ReaderPerm2)
            ELSEIF(i==3) THEN
              ReaderPerm3 = WorkVar % Perm
              ReaderValues3 = WorkVar % Values
              CALL VariableAdd(HydroSolver % Mesh % Variables, HydroSolver % & 
                Mesh, Model % Solvers(HPSolver), VarName, 1, ReaderValues3, &
                ReaderPerm3)
            ELSEIF(i==4) THEN
              ReaderPerm4 = WorkVar % Perm
              ReaderValues4 = WorkVar % Values
              CALL VariableAdd(HydroSolver % Mesh % Variables, HydroSolver % & 
                Mesh, Model % Solvers(HPSolver), VarName, 1, ReaderValues4, &
                ReaderPerm4)
            ELSEIF(i==5) THEN
              ReaderPerm5 = WorkVar % Perm
              ReaderValues5 = WorkVar % Values
              CALL VariableAdd(HydroSolver % Mesh % Variables, HydroSolver % & 
                Mesh, Model % Solvers(HPSolver), VarName, 1, ReaderValues5, &
                ReaderPerm5)
            ELSEIF(i==6) THEN
              ReaderPerm6 = WorkVar % Perm
              ReaderValues6 = WorkVar % Values
              CALL VariableAdd(HydroSolver % Mesh % Variables, HydroSolver % & 
                Mesh, Model % Solvers(HPSolver), VarName, 1, ReaderValues6, &
                ReaderPerm6)
            ELSEIF(i==7) THEN
              ReaderPerm7 = WorkVar % Perm
              ReaderValues7 = WorkVar % Values
              CALL VariableAdd(HydroSolver % Mesh % Variables, HydroSolver % & 
                Mesh, Model % Solvers(HPSolver), VarName, 1, ReaderValues7, &
                ReaderPerm7)
            ELSEIF(i==8) THEN
              ReaderPerm8 = WorkVar % Perm
              ReaderValues8 = WorkVar % Values
              CALL VariableAdd(HydroSolver % Mesh % Variables, HydroSolver % & 
                Mesh, Model % Solvers(HPSolver), VarName, 1, ReaderValues8, &
                ReaderPerm8)
            ELSEIF(i==9) THEN
              ReaderPerm9 = WorkVar % Perm
              ReaderValues9 = WorkVar % Values
              CALL VariableAdd(HydroSolver % Mesh % Variables, HydroSolver % & 
                Mesh, Model % Solvers(HPSolver), VarName, 1, ReaderValues9, &
                ReaderPerm9)
            ELSEIF(i==10) THEN
              ReaderPerm10 = WorkVar % Perm
              ReaderValues10 = WorkVar % Values
              CALL VariableAdd(HydroSolver % Mesh % Variables, HydroSolver % & 
                Mesh, Model % Solvers(HPSolver), VarName, 1, ReaderValues10, &
                ReaderPerm10)
            END IF

          END IF
        END DO
        
        NULLIFY(WorkVar)

      END IF!LoadReaders
    END IF!FirstTime

    !Time to interpolate variables that are solved on ice mesh
    HydroSolver => Model % Solvers(HPSolver)
    TempSolver => Model % Solvers(TSolver)
    ALLOCATE(InterpDim(1)); InterpDim = (/3/)
    ALLOCATE(BasalLogical(Model % Mesh % NumberOfNodes))
    ALLOCATE(IceMeshBasePerm(Model % Mesh % NumberOfNodes))
 
    CALL MakePermUsingMask(Model, Solver, Model % Mesh, "Bottom Surface Mask",&
         .FALSE., IceMeshBasePerm, DummyInt)

    !True nodes ignored by InterpolateVarToVarReduced 
    DO i=1, Model % Mesh % NumberOfNodes
      BasalLogical(i) = (IceMeshBasePerm(i) <=0)
    END DO

    !Set up list of variables needed by GlaDS that have to be interpolated
    InterpVar1 => VariableGet(Model % Mesh % Variables, "normalstress", ThisOnly=.TRUE., UnfoundFatal=.TRUE.)
    InterpVar2 => VariableGet(Model % Mesh % Variables, "velocity 1", ThisOnly=.TRUE., UnfoundFatal=.TRUE.)
    InterpVar3 => VariableGet(Model % Mesh % Variables, "velocity 2", ThisOnly=.TRUE., UnfoundFatal=.TRUE.)
    InterpVar4 => VariableGet(Model % Mesh % Variables, "groundedmask", ThisOnly=.TRUE., UnfoundFatal=.FALSE.)
    InterpVar5 => VariableGet(Model % Mesh % Variables, "gmcheck", ThisOnly=.TRUE., UnfoundFatal=.FALSE.)

    !Make copies of the relevant variables to save messing around with the mesh
    !variable list - only need perms and values
    ALLOCATE(InterpVar1Copy % Values(SIZE(InterpVar1 % Values)), InterpVar1Copy % Perm(SIZE(InterpVar1 % Perm)))
    InterpVar1Copy % Values = InterpVar1 % Values
    InterpVar1Copy % Perm = InterpVar1 % Perm
    InterpVar1Copy % Next => InterpVar2Copy
    InterpVar1Copy % Name = InterpVar1 % Name

    ALLOCATE(InterpVar2Copy % Values(SIZE(InterpVar2 % Values)), InterpVar2Copy % Perm(SIZE(InterpVar2 % Perm)))
    InterpVar2Copy % Values = InterpVar2 % Values
    InterpVar2Copy % Perm = InterpVar2 % Perm
    InterpVar2Copy % Next => InterpVar3Copy
    InterpVar2Copy % Name = InterpVar2 % Name

    ALLOCATE(InterpVar3Copy % Values(SIZE(InterpVar3 % Values)), InterpVar3Copy % Perm(SIZE(InterpVar3 % Perm)))
    InterpVar3Copy % Values = InterpVar3 % Values
    InterpVar3Copy % Perm = InterpVar3 % Perm
    IF(ASSOCIATED(InterpVar4)) THEN
      InterpVar3Copy % Next => InterpVar4Copy
    ELSE
      InterpVar3Copy % Next => NULL()
    END IF
    InterpVar3Copy % Name = InterpVar3 % Name

    IF(ASSOCIATED(InterpVar4)) THEN
      ALLOCATE(InterpVar4Copy % Values(SIZE(InterpVar4 % Values)), InterpVar4Copy % Perm(SIZE(InterpVar4 % Perm)))
      InterpVar4Copy % Values = InterpVar4 % Values
      InterpVar4Copy % Perm = InterpVar4 % Perm
      IF(ASSOCIATED(InterpVar5)) THEN
        InterpVar4Copy % Next => InterpVar5Copy
      ELSE
        InterpVar4Copy % Next => NULL()
      END IF
      InterpVar4Copy % Name = InterpVar4 % Name
    END IF

    IF(ASSOCIATED(InterpVar5)) THEN
      ALLOCATE(InterpVar5Copy % Values(SIZE(InterpVar5 % Values)), InterpVar5Copy % Perm(SIZE(InterpVar5 % Perm)))
      InterpVar5Copy % Values = InterpVar5 % Values
      InterpVar5Copy % Perm = InterpVar5 % Perm
      InterpVar5Copy % Next => NULL()
      InterpVar5Copy % Name = InterpVar5 % Name
    END IF

    InterpVar1 => InterpVar1Copy
    InterpVar2 => InterpVar2Copy
    InterpVar3 => InterpVar3Copy
    IF(ASSOCIATED(InterpVar4)) InterpVar4 => InterpVar4Copy
    IF(ASSOCIATED(InterpVar5)) InterpVar5 => InterpVar5Copy
    IF(FirstTime) THEN
      FirstTime=.FALSE.
      WorkVar => VariableGet(HydroSolver % Mesh % Variables,&
                 'hydraulic potential', ThisOnly=.TRUE.)
      ALLOCATE(NewPerm1(SIZE(WorkVar % Perm)),NewValues1(SIZE(WorkVar % Values)))
      ALLOCATE(NewPerm2(SIZE(WorkVar % Perm)),NewValues2(SIZE(WorkVar % Values)))
      ALLOCATE(NewPerm3(SIZE(WorkVar % Perm)),NewValues3(SIZE(WorkVar % Values)))
      ALLOCATE(NewPerm4(SIZE(WorkVar % Perm)),NewValues4(SIZE(WorkVar % Values)))
      ALLOCATE(NewPerm5(SIZE(WorkVar % Perm)),NewValues5(SIZE(WorkVar % Values)))
      NewPerm1 = WorkVar % Perm
      NewPerm2 = WorkVar % Perm
      NewPerm3 = WorkVar % Perm
      NewPerm4 = WorkVar % Perm
      NewPerm5 = WorkVar % Perm
      !NPP => NewPerm1
      NewValues1 = 0.0_dp
      NewValues2 = 0.0_dp
      NewValues3 = 0.0_dp
      NewValues4 = 0.0_dp
      NewValues5 = 0.0_dp
      !NVP => NewValues1
      CALL VariableAdd(HydroSolver % Mesh % Variables, HydroSolver % & 
           Mesh, HydroSolver, InterpVar1 % Name, 1,&
           NewValues1, NewPerm1)
      !NVP => NewValues2
      CALL VariableAdd(HydroSolver % Mesh % Variables, HydroSolver % & 
           Mesh, HydroSolver, InterpVar2 % Name, 1,&
           NewValues2, NewPerm2)
      !NVP => NewValues3
      CALL VariableAdd(HydroSolver % Mesh % Variables, HydroSolver % & 
           Mesh, HydroSolver, InterpVar3 % Name, 1,&
           NewValues3, NewPerm3)
      IF(ASSOCIATED(InterpVar4)) THEN
        !NVP => NewValues4
        CALL VariableAdd(HydroSolver % Mesh % Variables, HydroSolver % & 
             Mesh, HydroSolver, InterpVar4 % Name, 1,&
             NewValues4, NewPerm4)
      END IF
      IF(ASSOCIATED(InterpVar5)) THEN
        !NVP => NewValues5
        CALL VariableAdd(HydroSolver % Mesh % Variables, HydroSolver % & 
             Mesh, HydroSolver, InterpVar5 % Name, 1,&
             NewValues5, NewPerm5)
      END IF
    END IF

    !Have to divide temp residual by ice boundary weights before interpolating.
    !I initially did this by creating a new variable and interpolating that, but
    !it did some very odd things, such as randomly crashing the plume solver, so
    !this just alters the values of temp residual in place, interpolates them, 
    !and then restores the old values.
    WorkVar => VariableGet(Model % Mesh % Variables, 'temp residual', ThisOnly=.TRUE., UnfoundFatal=.TRUE.)
    ALLOCATE(NewValues6(SIZE(WorkVar % Values)))
    NewValues6 = 0.0_dp
    NewValues6 = WorkVar % Values
    !CALL CalculateNodalWeights(TempSolver, .TRUE., WorkVar % Perm, 'IceWeights')
    WorkVar2 => VariableGet(Model % Mesh % Variables, 'IceWeights', ThisOnly=.TRUE., UnfoundFatal=.TRUE.)
    !IF(ParEnv % PEs > 1) CALL ParallelSumVector(TempSolver % Matrix, WorkVar2 % Values)
    !DO i=1, Model % Mesh % NumberOfBoundaryElements!SIZE(WorkVar % Perm)
      !Element => Model % Mesh % Elements(Model % Mesh % NumberOfBulkElements+i)
      !n = GetElementNOFNodes(Element)
      !DO j=1, n
        !IF(WorkVar2 % Values(WorkVar2 % Perm(Element % NodeIndexes(j)))==0) CYCLE
        !WorkVar % Values(WorkVar % Perm(Element % NodeIndexes(j))) =&
        !WorkVar % Values(WorkVar % Perm(Element % NodeIndexes(j)))/&
        !WorkVar2 % Values(WorkVar2 % Perm(Element % NodeIndexes(j)))
      !END DO
    !END DO
    DO i=1, SIZE(WorkVar % Perm)
      IF(WorkVar2 % Values(WorkVar2 % Perm(i)) .NE. 0.0) THEN
        WorkVar % Values(WorkVar % Perm(i)) = WorkVar % Values(WorkVar % Perm(i))/WorkVar2 % Values(WorkVar2 % Perm(i))
      ELSE
        WorkVar % Values(WorkVar % Perm(i)) = 0.0
      END IF
    END DO

    CALL ParallelActive(.TRUE.)
    CALL InterpolateVarToVarReduced(Model % Mesh, HydroSolver % Mesh, 'temp residual',&
         InterpDim, OldNodeMask=BasalLogical, Variables=InterpVar1)
    CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

    !To restore temp residual values on ice mesh
    WorkVar % Values = NewValues6
    DEALLOCATE(NewValues6)

    !Interpolating GroundedMask will inevitably create values that aren't -1, 0
    !or 1, so need to round them to the nearest integer
    !Also, enforce grounded on upstream areas to deal with boundary
    !interpolation artefacts
    IF(ASSOCIATED(InterpVar5)) THEN
      WorkVar => VariableGet(HydroSolver % Mesh % Variables, "gmcheck", ThisOnly=.TRUE., UnfoundFatal=.TRUE.)
      DO i=1, SIZE(WorkVar % Values)
        WorkVar % Values(i) = ANINT(WorkVar % Values(i))
      END DO
    END IF
    IF(ASSOCIATED(InterpVar4)) THEN
      WorkVar => VariableGet(HydroSolver % Mesh % Variables, "groundedmask", ThisOnly=.TRUE., UnfoundFatal=.TRUE.)
      DO i=1, SIZE(WorkVar % Values)
        WorkVar % Values(i) = ANINT(WorkVar % Values(i))
      END DO

      RefNode => ListGetIntegerArray(Params, 'Reference Node', Found)
      IF(Found) THEN
        Threshold = GetConstReal(Params, 'Threshold Distance', Found)
        IF(.NOT. Found) Threshold = 10000.0

        DO i=1, SIZE(WorkVar % Perm)
          Dist = (HydroSolver % Mesh % Nodes % x(WorkVar % Perm(i)) -&
                  RefNode(1))**2
          Dist = Dist + (HydroSolver % Mesh % Nodes % y(WorkVar % Perm(i)) -&
                  RefNode(2))**2
          Dist = SQRT(Dist)
          IF(Dist > Threshold) WorkVar % Values(WorkVar % Perm(i)) = 1.0
        END DO
      END IF

      Hit = .FALSE.
      SideBC1 = 0
      SideBC2 = 0
      HitCount = 0
      DO i=1, Model % NumberOfBCs
        BC => Model % BCs(i) % Values
        Hit = ListGetLogical(BC, "Side", Found)
        IF(Hit) THEN
          IF(HitCount == 0) THEN
            SideBC1 = i
          ELSEIF(HitCount == 1) THEN
            SideBC2 = i
          ELSE
            CALL Info('CalvingHydroInterp', 'You appear to have more than two&
                       lateral boundaries. Are you sure about this?')
          END IF
          Hit = .FALSE.
          HitCount = HitCount+1
          IF(HitCount == 2) EXIT
        END IF
      END DO

      DO i=1, HydroSolver % Mesh % NumberOfBoundaryElements
        Element => GetBoundaryElement(i, HydroSolver)
        ElementBC = GetBCId(Element)
        IF(ElementBC == SideBC1 .OR. ElementBC == SideBC2) THEN
          n = GetElementNOFNodes(Element)
          DO j=1,n
            IF(WorkVar % Values(WorkVar % Perm(Element % NodeIndexes(j))) == -1.0) THEN
              WorkVar % Values(WorkVar % Perm(Element % NodeIndexes(j))) = 1.0
            END IF
          END DO
        END IF
      END DO
    END IF

    !This is to remove boundary artefacts in normalstress, which, otherwise,
    !really mess the hydrology up. Artefacts do exist in velocity and temp
    !residual, but have much less impact
    !Routine will search for nearest node that has a non-zero value of
    !normalstress and substitute that value for the erroneous zero value.
    !Ideally, should interpolate from nearest non-zero nodes, but not worth the
    !extra faff for a minimal gain

    !This first section is due to SolveLinearSystem invalidating the
    !Force2Stress variable (i.e. normalstress) on the hydro mesh. Not really
    !sure why it does this, but this fix seems to work without any knock-on
    !effects.
    WorkVar2 => HydroSolver % Mesh % Variables
    DO WHILE (ASSOCIATED(WorkVar2))
      IF (TRIM(WorkVar2 % Name) == 'normalstress') THEN
        IF (.NOT. WorkVar2 % Valid) THEN
          WorkVar2 % Valid = .TRUE.
          WorkVar2 % PrimaryMesh => HydroSolver % Mesh
        END IF
        EXIT
      END IF
      WorkVar2 => WorkVar2 % Next
    END DO

    WorkVar2 => VariableGet(HydroSolver % Mesh % Variables, "normalstress", ThisOnly=.TRUE., UnfoundFatal=.TRUE.)
    MeanNS = 0.0_dp
    !ZeroCounter = 0

    !DO i=1, HydroSolver % Mesh % NumberOfBoundaryElements
    !  Element => GetBoundaryElement(i, HydroSolver)
    !  n = GetElementNOFNodes(Element)
    !  IF(ANY(WorkVar % Values(WorkVar % Perm(Element % NodeIndexes(1:n))) == -1.0)) CYCLE
    !  DO j=1,n
    !    MeanNS = MeanNS+WorkVar2 % Values(WorkVar2 % Perm(Element %&
    !    NodeIndexes(j)))
     !   IF(WorkVar2 % Values(WorkVar2 % Perm(Element % NodeIndexes(j)))==0.0)&
     !     THEN
     !     ZeroCounter = ZeroCounter + 1
     !   END IF
     ! END DO
    !END DO
    !MeanNS = MeanNS/(HydroSolver % Mesh % NumberOfBoundaryElements-ZeroCounter)

    DO i=1, HydroSolver % Mesh % NumberOfBoundaryElements
      Element => GetBoundaryElement(i, HydroSolver)
      n = GetElementNOFNodes(Element)
      IF(ANY(WorkVar % Values(WorkVar % Perm(Element % NodeIndexes(1:n))) == -1.0)) CYCLE
      ALLOCATE(NSValues(n))
      DO j=1,n
        NSValues(j) = WorkVar2 % Values(WorkVar2 % Perm(Element % NodeIndexes(j)))
      END DO
      ZeroCounter = 0
      DO j=1,n
        IF(NSValues(j) .NE. 0.0) CYCLE
        ZeroCounter = ZeroCounter + 1
      END DO
      !This will currently only work for elements with 2 or 3 nodes
      SELECT CASE (n-ZeroCounter)
        CASE (0)
          CALL Info('CalvingHydroInterp', 'No non-0 values of NormalStress. &
                   Making a guess')
          NewNS = (SUM(WorkVar2 % Values)/SIZE(WorkVar2 % Values))+3.0_dp !MeanNS
        CASE (1)
          NewNS = SUM(NSValues) 
        CASE (2)
          NewNS = SUM(NSValues)/2
        CASE DEFAULT
          CALL Info('CalvingHydroInterp', 'No NormalStress correction possible')
          PRINT *, 'Diagnostics: ',n,ZeroCounter,Element % NodeIndexes
      END SELECT
      DO j=1,n
        IF(NSValues(j) .NE. 0.0) CYCLE
        WorkVar2 % Values(WorkVar2 % Perm(Element % NodeIndexes(j))) = NewNS
      END DO
      DEALLOCATE(NSValues)
    END DO

    !Finally, the area of the hydromesh beyond the calving front has to be
    !forced to be ungrounded, so interpolation artefacts there have to be
    !cleared away. The main job is done by InterpVarToVar, but there are always
    !a few places where it messes up, which are corrected here.
    DO i=1, HydroSolver % Mesh % NumberOfBulkElements
      Element => HydroSolver % Mesh % Elements(i)
      n = GetElementNOFNodes(Element)
      DO j=1, n
        IF(WorkVar2 % Values(WorkVar2 % Perm(Element % NodeIndexes(j))) == 0.0) THEN
          WorkVar % Values(WorkVar % Perm(Element % NodeIndexes(j))) = -1.0
        END IF
        
      END DO
    END DO    

    !Temp residual needs to be conserved. Here, just integrate across all
    !elements and compare totals, then scale values on hydromesh uniformly to
    !bring in line with ice mesh
    WorkVar => VariableGet(Model % Mesh % Variables, "temp residual", ThisOnly=.TRUE., UnfoundFatal=.TRUE.)

    IceTempResSum = 0.0_dp
    IceTempResSum = SUM(WorkVar % Values)

    !WorkVar => VariableGet(HydroSolver % Mesh % Variables, "hydraulic potential", ThisOnly=.TRUE., UnfoundFatal=.TRUE.)
    WorkVar => VariableGet(HydroSolver % Mesh % Variables, "temp residual", ThisOnly=.TRUE., UnfoundFatal=.TRUE.)
    !WorkVar => VariableGet(HydroSolver % Mesh % Variables, "WeightedTR", ThisOnly=.TRUE., UnfoundFatal=.TRUE.)
    !CALL CalculateNodalWeights(HydroSolver, .FALSE., WorkVar % Perm, 'HydroWeights')
    !WorkVar => VariableGet(HydroSolver % Mesh % Variables, "temp residual", ThisOnly=.TRUE., UnfoundFatal=.TRUE.)
    WorkVar2 => VariableGet(HydroSolver % Mesh % Variables, "HydroWeights", ThisOnly=.TRUE., UnfoundFatal=.TRUE.)
    !IF(ParEnv % PEs > 1) CALL ParallelSumVector(HydroSolver % Matrix, WorkVar2 % Values)
    DO i=1,SIZE(WorkVar % Perm)
      !Element => HydroSolver % Mesh % Elements(i)
      !n = GetElementNOFNodes(Element)
      !DO j=1, n
        !WorkVar % Values(WorkVar % Perm(Element % NodeIndexes(j))) =&
        !WorkVar % Values(WorkVar % Perm(Element % NodeIndexes(j)))*&
        !WorkVar2 % Values(WorkVar2 % Perm(Element % NodeIndexes(j)))
      !END DO
      WorkVar % Values(WorkVar % Perm(i)) =&
      WorkVar % Values(WorkVar % Perm(i))*WorkVar2 % Values(WorkVar2 % Perm(i))
    END DO
    HydroTempResSum = 0.0_dp
    HydroTempResSum = SUM(WorkVar % Values)

    ALLOCATE(ParITRS(ParEnv % PEs), ParHTRS(ParEnv % PEs))
    ParITRS = 0.0_dp
    ParHTRS = 0.0_dp

    IF(ParEnv % PEs > 1) THEN
      CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)
      CALL MPI_Gather(IceTempResSum, 1, MPI_DOUBLE_PRECISION, ParITRS, 1, MPI_DOUBLE_PRECISION, 0, ELMER_COMM_WORLD, ierr)
      CALL MPI_Gather(HydroTempResSum, 1, MPI_DOUBLE_PRECISION, ParHTRS, 1, MPI_DOUBLE_PRECISION, 0, ELMER_COMM_WORLD, ierr)
      IF(ParEnv % myPE == 0) THEN
        IF(ANINT(SUM(ParITRS)) .NE. ANINT(SUM(ParHTRS))) THEN
          ScaleFactor = SUM(ParITRS)/SUM(ParHTRS)
        END IF
      END IF
      CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)
      CALL MPI_Bcast(ScaleFactor, 1, MPI_DOUBLE_PRECISION, 0, ELMER_COMM_WORLD, ierr)
      DO i=1, SIZE(WorkVar % Values)
        WorkVar % Values(i) = WorkVar % Values(i)*ScaleFactor
      END DO
    ELSE
      IF(ANINT(IceTempResSum) .NE. ANINT(HydroTempResSum)) THEN
        ScaleFactor = IceTempResSum/HydroTempResSum
        DO i=1, SIZE(WorkVar % Values)
          WorkVar % Values(i) = WorkVar % Values(i)*ScaleFactor
        END DO
      END IF
    END IF

    DEALLOCATE(InterpDim, BasalLogical, IceMeshBasePerm, ParITRS, ParHTRS)
    DEALLOCATE(InterpVar1Copy % Values, InterpVar2Copy % Values, InterpVar3Copy % Values)
    DEALLOCATE(InterpVar1Copy % Perm, InterpVar2Copy % Perm, InterpVar3Copy % Perm)
    IF(ASSOCIATED(InterpVar4)) DEALLOCATE(InterpVar4Copy % Values, InterpVar4Copy % Perm)
    IF(ASSOCIATED(InterpVar5)) DEALLOCATE(InterpVar5Copy % Values, InterpVar5Copy % Perm)
    NULLIFY(InterpVar1, InterpVar2, InterpVar3, InterpVar4, InterpVar5,&
           WorkVar, HydroSolver, NVP, NPP, BC, Element, RefNode, WorkVar2,&
           TempSolver)
!-------------------------------------------------------------------------------
  CONTAINS
!-------------------------------------------------------------------------------
    CHARACTER(len=20) FUNCTION STR(k)
      !Converts and integer to a string
      INTEGER, INTENT(IN) :: k
      WRITE (STR, *) k
      STR = adjustl(STR)
    END FUNCTION STR
!-------------------------------------------------------------------------------
  END SUBROUTINE

! *****************************************************************************/
!> Interpolates hydrology variables necessary for calving from hydro mesh to ice
!> mesh as part of combined calving-hydrology simulations.
   SUBROUTINE HydroToIceInterp( Model,Solver,Timestep,TransientSimulation )
!******************************************************************************
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh,materials,BCs,etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear equation solver options
!
!  REAL(KIND=dp) :: Timestep
!     INPUT: Timestep size for time dependent simulations
!
!******************************************************************************
     USE Differentials
     USE MaterialModels
     USE DefUtils
     USE InterpVarToVar
!------------------------------------------------------------------------------
     IMPLICIT NONE
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!    External variables
!------------------------------------------------------------------------------
     TYPE(Model_t)  :: Model
     TYPE(Solver_t), TARGET :: Solver
     LOGICAL :: TransientSimulation
     REAL(KIND=dp) :: Timestep
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Mesh_t), POINTER :: IceMesh, HydroMesh
     TYPE(Solver_t), POINTER ::  HydroSolver
     TYPE(Variable_t), POINTER :: WorkVar=>NULL(), InterpVar1=>NULL(),&
                                  InterpVar2=>NULL(), InterpVar3=>NULL()
     TYPE(Variable_t), TARGET :: InterpVar1Copy, InterpVar2Copy, InterpVar3Copy
     TYPE(ValueList_t), POINTER :: Params
     LOGICAL :: FirstTime=.TRUE.
     LOGICAL, POINTER :: BasalLogical(:)
     INTEGER, POINTER :: InterpDim(:), IceMeshBasePerm(:)=>NULL(),&
                         NPP(:)=>NULL(), NewPerm1(:)
     INTEGER, ALLOCATABLE, TARGET :: NewPerm2(:), NewPerm3(:)
     INTEGER :: i, HPSolver, DummyInt, ierr
     REAL(KIND=dp), POINTER :: NVP(:)=>NULL(), NewValues1(:)
     REAL(KIND=dp), ALLOCATABLE, TARGET :: NewValues2(:),&
                                           NewValues3(:)
     SAVE HydroMesh, FirstTime, HPSolver! NewPerm1,&
          !NewValues1, NewPerm2, NewValues2, NewPerm3, NewValues3

!------------------------------------------------------------------------------

    Params => GetSolverParams()

    !For finding main hydrology solver so know where to interpolate from.
    IF(FirstTime) THEN
      DO i=1,Model % NumberOfSolvers
        IF(Model % Solvers(i) % Variable % Name == 'hydraulic potential') THEN
          HPSolver = i
          EXIT
        END IF
      END DO
    END IF!FirstTime

    !Time to interpolate variables that are solved on ice mesh
    HydroSolver => Model % Solvers(HPSolver)
    ALLOCATE(InterpDim(1)); InterpDim = (/3/)
    ALLOCATE(BasalLogical(Model % Mesh % NumberOfNodes))
    ALLOCATE(IceMeshBasePerm(Model % Mesh % NumberOfNodes))
 
    CALL MakePermUsingMask(Model, Solver, Model % Mesh, "Bottom Surface Mask",&
         .FALSE., IceMeshBasePerm, DummyInt)

    !True nodes ignored by InterpolateVarToVarReduced 
    DO i=1, Model % Mesh % NumberOfNodes
      BasalLogical(i) = (IceMeshBasePerm(i) <=0)
    END DO

    !Set up list of variables needed by GlaDS that have to be interpolated
    InterpVar1 => VariableGet(HydroSolver % Mesh % Variables, "water pressure", ThisOnly=.TRUE., UnfoundFatal=.TRUE.)
    !InterpVar2 => VariableGet(HydroSolver % Mesh % Variables, "sheet discharge 1", ThisOnly=.TRUE., UnfoundFatal=.TRUE.)
    !InterpVar3 => VariableGet(HydroSolver % Mesh % Variables, "sheet discharge 2", ThisOnly=.TRUE., UnfoundFatal=.TRUE.)

    !Make copies of the relevant variables to save messing around with the mesh
    !variable list - only need perms and values
    ALLOCATE(InterpVar1Copy % Values(SIZE(InterpVar1 % Values)), InterpVar1Copy % Perm(SIZE(InterpVar1 % Perm)))
    InterpVar1Copy % Values = InterpVar1 % Values
    InterpVar1Copy % Perm = InterpVar1 % Perm
    InterpVar1Copy % Next => NULL() !InterpVar2Copy
    InterpVar1Copy % Name = InterpVar1 % Name

    !ALLOCATE(InterpVar2Copy % Values(SIZE(InterpVar2 % Values)), InterpVar2Copy % Perm(SIZE(InterpVar2 % Perm)))
    !InterpVar2Copy % Values = InterpVar2 % Values
    !InterpVar2Copy % Perm = InterpVar2 % Perm
    !InterpVar2Copy % Next => InterpVar3Copy
    !InterpVar2Copy % Name = InterpVar2 % Name

    !ALLOCATE(InterpVar3Copy % Values(SIZE(InterpVar3 % Values)), InterpVar3Copy % Perm(SIZE(InterpVar3 % Perm)))
    !InterpVar3Copy % Values = InterpVar3 % Values
    !InterpVar3Copy % Perm = InterpVar3 % Perm
    !InterpVar3Copy % Next => NULL()
    !InterpVar3Copy % Name = InterpVar3 % Name

    InterpVar1 => InterpVar1Copy
    !InterpVar2 => InterpVar2Copy
    !InterpVar3 => InterpVar3Copy

    !Variables need to be added to mesh before interpolated as list
    IF(FirstTime) THEN
      FirstTime = .FALSE.
      WorkVar => VariableGet(Model % Mesh % Variables,&
                 'velocity 1', ThisOnly=.TRUE.)
      ALLOCATE(NewPerm1(SIZE(WorkVar % Perm)),NewValues1(SIZE(WorkVar % Values)))
      !ALLOCATE(NewPerm2(SIZE(WorkVar % Perm)),NewValues2(SIZE(WorkVar % Values)))
      !ALLOCATE(NewPerm3(SIZE(WorkVar % Perm)),NewValues3(SIZE(WorkVar % Values)))
      NewPerm1 = WorkVar % Perm
      !NewPerm2 = NewPerm1
      !NewPerm3 = NewPerm1
      !NPP => NewPerm1
      NewValues1 = 0.0_dp
      !NewValues2 = 0.0_dp
      !NewValues3 = 0.0_dp
      !NVP => NewValues1
      CALL VariableAdd(Model % Mesh % Variables, Model % & 
           Mesh, CurrentModel % Solver, InterpVar1 % Name, 1,&
           NewValues1, NewPerm1)
      !NPP => NewPerm2
      !NVP => NewValues2
      !CALL VariableAdd(Model % Mesh % Variables, Model % & 
      !     Mesh, CurrentModel % Solver, InterpVar2 % Name, 1,&
      !     NVP, NPP)
      !NPP => NewPerm3
      !NVP => NewValues3
      !CALL VariableAdd(Model % Mesh % Variables, Model % & 
      !     Mesh, CurrentModel % Solver, InterpVar3 % Name, 1,&
      !     NVP, NPP)
    END IF

    CALL ParallelActive(.TRUE.)
    CALL InterpolateVarToVarReduced(HydroSolver % Mesh, Model % Mesh,&
         'effective pressure', InterpDim, NewNodeMask=BasalLogical,&
         Variables=InterpVar1)
    CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)
    CALL ParallelActive(.FALSE.)

    DEALLOCATE(InterpDim, BasalLogical, IceMeshBasePerm, InterpVar1Copy % Perm,&
               InterpVar1Copy % Values)!, InterpVar2Copy % Perm,&
               !InterpVar2Copy % Values, InterpVar3Copy % Perm,&
               !InterpVar3Copy % Values)
    NULLIFY(HydroSolver, WorkVar, NVP, NPP, InterpVar1)!, InterpVar2, InterpVar3)

  END SUBROUTINE
! *****************************************************************************/
!> Works out weights on hydro mesh at beginning of simulation
   SUBROUTINE HydroWeightsSolver( Model,Solver,Timestep,TransientSimulation )
!******************************************************************************
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh,materials,BCs,etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear equation solver options
!
!  REAL(KIND=dp) :: Timestep
!     INPUT: Timestep size for time dependent simulations
!
!******************************************************************************
     USE Differentials
     USE MaterialModels
     USE DefUtils
     USE InterpVarToVar
!------------------------------------------------------------------------------
     IMPLICIT NONE
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!    External variables
!------------------------------------------------------------------------------
     TYPE(Model_t)  :: Model
     TYPE(Solver_t), TARGET :: Solver
     LOGICAL :: TransientSimulation
     REAL(KIND=dp) :: Timestep
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Variable_t), POINTER :: WorkVar
!------------------------------------------------------------------------------    
    CALL CalculateNodalWeights(Solver, .FALSE.,VarName='HydroWeights')
    WorkVar => VariableGet(Solver % Mesh % Variables, "HydroWeights", ThisOnly=.TRUE., UnfoundFatal=.TRUE.)
    IF(ParEnv % PEs > 1) CALL ParallelSumVector(Solver % Matrix, WorkVar % Values)
    NULLIFY(WorkVar)
END SUBROUTINE
! *****************************************************************************/
!> Works out weights on ice mesh
   SUBROUTINE IceWeightsSolver( Model,Solver,Timestep,TransientSimulation )
!******************************************************************************
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh,materials,BCs,etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear equation solver options
!
!  REAL(KIND=dp) :: Timestep
!     INPUT: Timestep size for time dependent simulations
!
!******************************************************************************
     USE Differentials
     USE MaterialModels
     USE DefUtils
     USE InterpVarToVar
!------------------------------------------------------------------------------
     IMPLICIT NONE
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!    External variables
!------------------------------------------------------------------------------
     TYPE(Model_t)  :: Model
     TYPE(Solver_t), TARGET :: Solver
     LOGICAL :: TransientSimulation
     REAL(KIND=dp) :: Timestep
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Variable_t), POINTER :: WorkVar
!------------------------------------------------------------------------------    
    CALL CalculateNodalWeights(Solver, .TRUE.,VarName='IceWeights')
    WorkVar => VariableGet(Solver % Mesh % Variables, "IceWeights", ThisOnly=.TRUE., UnfoundFatal=.TRUE.)
    IF(ParEnv % PEs > 1) CALL ParallelSumVector(Solver % Matrix, WorkVar % Values)
    NULLIFY(WorkVar)
END SUBROUTINE
