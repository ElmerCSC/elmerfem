!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
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
!
!/******************************************************************************
! *
! *  Authors:   Eelis Takala(Trafotek Oy) and Juha Ruokolainen(CSC)
! *  Emails:    eelis.takala@trafotek.fi and Juha.Ruokolainen@csc.fi
! *  Web:       http://www.trafotek.fi and http://www.csc.fi/elmer
! *  Addresses: Trafotek Oy
! *             Kaarinantie 700
! *             Turku
! *
! *             and
! *
! *             CSC - IT Center for Science Ltd.
! *             Keilaranta 14
! *             02101 Espoo, Finland 
! *
! *  Original Date: October 2015
! *
! *****************************************************************************/
 
!------------------------------------------------------------------------------
!> Initialization for the primary solver: CurrentSource
!------------------------------------------------------------------------------
SUBROUTINE CircuitsAndDynamics_init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: TransientSimulation !< Steady state or transient simulation
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: RotMachine, Found
  
  Params => Solver % Values

  ! This is only created if no variable present!
  CALL ListAddNewString( Params,'Variable','-global ckt')

  ! The circuit variable is no ordinary variable!
  CALL ListAddNewLogical( Params,'Variable Output',.FALSE.)
  CALL ListAddNewInteger( Params,'Time Derivative Order',1)
  CALL ListAddNewLogical( Params,'No Matrix',.TRUE.)

  ! When we introduce the variables in this way the variables are created
  ! so that they exist when the proper simulation cycle starts.
  ! This also keeps the command file cleaner. 
  ! If "Rotor Angle" is created in Simulation section no need to do it here!
  IF( .NOT. ListCheckPresent( Model % Simulation,'Rotor Angle') ) THEN
    RotMachine = ListGetLogical( Params,'Rotating Machine',Found )
    IF(.NOT. Found ) THEN
      RotMachine = ListGetLogicalAnyBC(Model,'Rotational Projector') .OR. &
          ListGetLogicalAnyBC(Model,'Anti Rotational Projector') 
    END IF
    
    IF( RotMachine ) THEN
      CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
          '-global Rotor Angle' )
      IF( TransientSimulation ) THEN
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
            '-global Rotor Velo' )
      END IF
    END IF
  END IF
    
  Solver % Values => Params
!------------------------------------------------------------------------------
END SUBROUTINE CircuitsAndDynamics_init
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
SUBROUTINE CircuitsAndDynamics( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE CircuitUtils
  USE CircuitsMod
  USE CircMatInitMod
  USE MGDynMaterialUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: TransientSimulation !< Steady state or transient simulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  LOGICAL :: First=.TRUE.
  TYPE(Solver_t), POINTER :: Asolver => Null()
  INTEGER :: p, n, istat, max_element_dofs, i, j
  TYPE(Mesh_t), POINTER :: Mesh  
  TYPE(Matrix_t), POINTER :: CM
  INTEGER, POINTER :: n_Circuits => Null()
  TYPE(Circuit_t), POINTER :: Circuits(:)  
  REAL(KIND=dp), ALLOCATABLE :: Crt(:)     
  TYPE(Variable_t), POINTER :: LagrangeVar
  INTEGER :: Tstep=-1
  LOGICAL :: Parallel
  REAL(KIND=dp), POINTER :: px(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: MultName
  LOGICAL :: Found
  CHARACTER(LEN=MAX_NAME_LEN) :: sname
  CHARACTER(*), PARAMETER :: Caller = 'CircuitsAndDynamics'

  SAVE First, Tstep, Parallel, Crt, MultName
  
!------------------------------------------------------------------------------
  
  IF (First) THEN
    IF( TransientSimulation ) THEN
      CALL Info(Caller,'Initializing electric circuits for transient simulation',Level=6)
    ELSE
      CALL Info(Caller,'Initializing electric circuits for steady state simulation',Level=6)
    END IF
    
    First = .FALSE.
    
    Parallel = Solver % Parallel
    IF(Parallel) THEN
      CALL Info(Caller,'Assuming parallel electric circuits',Level=12)
    ELSE
      CALL Info(Caller,'Assuming serial electric circuits',Level=12)
    END IF
      
    Model % HarmonicCircuits = .FALSE.
    CALL AddComponentsToBodyLists()
    
    ALLOCATE( Model % Circuit_tot_n, Model % n_Circuits, STAT=istat )
    IF ( istat /= 0 ) THEN
      CALL Fatal(Caller, 'Memory allocation error.' )
    END IF

    n_Circuits => Model % n_Circuits
    Model % Circuit_tot_n = 0

    ! Look for the real valued solver we attach the circuit equations to:
    ! -------------------------------------------------------------------
    Asolver => NULL()
    DO i=1,Model % NumberOfSolvers      
      sname = GetString(Model % Solvers(i) % Values, 'Procedure', Found)
      j = INDEX( sname,'MagnetoDynamics2D')
      IF(j>0) THEN
        IF( INDEX( sname,'MagnetoDynamics2DHarmonic') > 0 ) CYCLE
      END IF
      IF(j==0) j = INDEX( sname,'WhitneyAVSolver')
      IF( j > 0 ) THEN
        ASolver => Model % Solvers(i) 
        EXIT
      END IF
    END DO

    IF( Solver % Parallel .NEQV. Asolver % Parallel  ) THEN
      CALL Warn(Caller,'Conflicting parallel status for circuit and A solver!')
      Solver % Parallel = .TRUE.
      ASolver % Parallel = .TRUE.
    END IF
    
    IF(.NOT. ASSOCIATED(ASolver) ) THEN
      ASolver => FindSolverWithKey('Export Lagrange Multiplier')
    END IF
    CALL Info(Caller,'Circuit equations associated with solver index: '&
        //I2S(ASolver % SolverId),Level=6)
    Model % ASolver => ASolver 
       
    CALL AllocateCircuitsList() ! CurrentModel%Circuits
    Circuits => Model % Circuits

    CALL SetBoundaryAreasToValueLists() 

    DO p=1,n_Circuits
      n = GetNofCircVariables(p)

      CALL Info(Caller,'Initializing circuit '//I2S(p)//' with '//I2S(n)//' variables!',Level=6)

      CALL AllocateCircuit(p)
      
      Circuits(p) % n_comp = CountNofCircComponents(p, n)
      ALLOCATE(Circuits(p) % Components(Circuits(p) % n_comp))
      
      Circuits(p) % Harmonic = .FALSE.
      Circuits(p) % Parallel = Parallel
      
      CALL ReadCircuitVariables(p)
      CALL ReadComponents(p)
      CALL AddComponentValuesToLists(p)  ! Lists are used to communicate values to other solvers at the moment...
      CALL AddBareCircuitVariables(p)   ! these don't belong to any components
      CALL ReadCoefficientMatrices(p)
      CALL ReadPermutationVector(p)
      CALL ReadCircuitSources(p)
      CALL WriteCoeffVectorsForCircVariables(p)
    
    END DO

    CALL CheckComponentVariables()

    
    ! Create CRS matrix structures for the circuit equations:
    ! ------------------------------------------------------
    CALL Circuits_MatrixInit()
    ALLOCATE(Crt(Model % Circuit_tot_n))

    MultName = LagrangeMultiplierName(ASolver)
  END IF
  
  ! If we have angle given explicitly, do not compute it 
  IF( .NOT. ListCheckPresent( Model % Simulation,'Rotor Angle') ) THEN
    IF( TransientSimulation ) CALL SetDynamicAngle()
  END IF
      
  IF (Tstep /= GetTimestep()) THEN
    Tstep = GetTimestep()
    ! Circuit variable values from previous timestep:
    ! -----------------------------------------------
    Crt = 0._dp

    LagrangeVar => VariableGet( Solver % Mesh % Variables, MultName )

    IF(ASSOCIATED(LagrangeVar)) THEN
      n = SIZE( LagrangeVar % Values )

      ! Debugging stuff activated only when "Max Output Level" >= 25
      IF( InfoActive( 25 ) ) THEN
        CALL VectorValuesRange(LagrangeVar % Values,n,TRIM(LagrangeVar % Name))       
      END IF
      
      IF( n < Model % Circuit_tot_n ) THEN
        CALL Fatal(Caller,'Lagrange multiplier is too small for circuits!')
      END IF
      
      ! We want to associate the variable of this solver to the LagrangeVar so that the library routines take
      ! care of the evolution of LagrangeVar for PrevValues and for parallel timestepping.
      IF(.NOT. ASSOCIATED( LagrangeVar, Solver % Variable ) ) THEN
        CALL Info(Caller,'Associating circuit variable to Lagrange values!',Level=8)
        Solver % Variable => LagrangeVar
      END IF
        
      IF( .NOT. ASSOCIATED( LagrangeVar % PrevValues ) ) THEN
        CALL Info(Caller,'Add PrevValues to Lagrange multiplier!',Level=8)
        ALLOCATE( LagrangeVar % PrevValues(n,1) )
      END IF

      ! Rotate solution here, as InitializeTimestep() doesn't do anything,  with 'no matrix' solvers...
      LagrangeVar % PrevValues(:,1) = LagrangeVar % Values 
        
      Crt = LagrangeVar % PrevValues(1:Model%Circuit_tot_n,1)
    END IF
    
    CALL Circuits_ToMeshVariable(Solver,crt) 
  END IF

  max_element_dofs = Model % Mesh % MaxElementDOFs
  Circuits => Model % Circuits
  n_Circuits => Model % n_Circuits
  CM => Model % CircuitMatrix
  
  ! Initialize Circuit matrix:
  ! -----------------------------
  IF(.NOT.ASSOCIATED(CM)) RETURN

  CM % RHS = 0._dp
  IF(ASSOCIATED(CM % Values)) CM % Values = 0._dp
  
  ! Write Circuit equations:
  ! ------------------------
  DO p = 1,n_Circuits
    CALL AddBasicCircuitEquations(p,Crt,dt)
    CALL AddComponentEquationsAndCouplings(p, max_element_dofs,dt,Crt)
  END DO
       
  IF(.NOT. ASSOCIATED( CM ) ) THEN
    Asolver %  Matrix % AddMatrix => NULL()
    RETURN
  END IF
 
  Asolver %  Matrix % AddMatrix => CM  
  IF(  CM % FORMAT == MATRIX_LIST ) CALL List_toCRSMatrix(CM)

  IF( InfoActive(20) ) THEN
    px => CM % Values
    CALL VectorValuesRange(px,SIZE(px),'CircuitMatrix',.TRUE.)
    px => CM % rhs
    CALL VectorValuesRange(px,SIZE(px),'CircuitRhs',.TRUE.)
  END IF

  IF( LIstGetLogical( Solver % Values,'Save Circuit Matrix',Found ) ) THEN
    CALL SaveLinearSystem( Solver, CM,'circuit',ASolver % Matrix % NumberOfRows)
  END IF
    
  
  IF(CM % NumberOfRows<=0)  THEN
    CALL FreeMatrix(CM)
    Asolver % Matrix % AddMatrix => NULL()
  END IF

  CALL Info(Caller,'Finished assembly of circuit matrix',Level=12)
  
CONTAINS

    
    
!------------------------------------------------------------------------------
   SUBROUTINE AddBasicCircuitEquations(p,Crt,dt)
!------------------------------------------------------------------------------
     IMPLICIT NONE
    INTEGER :: p
    REAL(KIND=dp) :: Crt(:)    
    REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------     
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(CircuitVariable_t), POINTER :: Cvar
    TYPE(ValueList_t), POINTER :: BF, Params
    TYPE(Matrix_t), POINTER :: CM
    INTEGER :: i, nm, RowId, ColId, j
    REAL(KIND=dp) :: vphi
    LOGICAL :: Found
    
    Circuit => CurrentModel % Circuits(p)
    nm = CurrentModel % Asolver % Matrix % NumberOfRows
    BF => CurrentModel % BodyForces(1) % Values
    CM => CurrentModel%CircuitMatrix
  
    Params => GetSolverParams()

    DO i=1,Circuit % n
      Cvar => Circuit % CircuitVariables(i)

      IF( Parallel ) THEN
        IF(Cvar % Owner /= ParEnv % myPE) CYCLE
      END IF
        
      RowId = Cvar % ValueId + nm
      
      IF( LEN_TRIM( Circuit % Source(i) ) > 0 ) THEN
        vphi = GetCReal(Params, Circuit % Source(i), Found)
        IF ( .NOT. Found .AND. ASSOCIATED(BF) ) THEN
          vphi = GetCReal(BF, Circuit % Source(i), Found)
        END IF
        IF (Found) Cvar % SourceRe(i) = vphi
      ELSE
        vphi = 0.0_dp
      END IF
      
      Cvar % SourceRe(i) = vphi
      CM % RHS(RowId) = Cvar % SourceRe(i)
        
      DO j=1,Circuit % n

        ColId = Circuit % CircuitVariables(j) % ValueId + nm

        IF ( TransientSimulation ) THEN 
          ! A d/dt(x): (x could be voltage or current):
          !--------------------------------------------
          IF(Cvar % A(j) /= 0._dp) THEN
            CALL AddToMatrixElement(CM, RowId, ColId, Cvar % A(j)/dt)
            CM % RHS(RowId) = CM % RHS(RowId) + Cvar % A(j) * Crt(ColId-nm) / dt
          END IF
        END IF  
        ! B x:
        ! ------
        IF(Cvar % B(j) /= 0._dp) THEN
          CALL AddToMatrixElement(CM, RowId, ColId, Cvar % B(j))
        END IF
      END DO
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE AddBasicCircuitEquations
!------------------------------------------------------------------------------


   
!------------------------------------------------------------------------------
   SUBROUTINE AddComponentEquationsAndCouplings(p, nn, dt, crt)
!------------------------------------------------------------------------------
    USE MGDynMaterialUtils
    IMPLICIT NONE
    INTEGER :: p, CompInd, nm, nn, nd
    TYPE(Solver_t), POINTER :: ASolver
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(Matrix_t), POINTER :: CM
    TYPE(Component_t), POINTER :: Comp
    TYPE(CircuitVariable_t), POINTER :: Cvar
    TYPE(Valuelist_t), POINTER :: CompParams
    TYPE(Element_t), POINTER :: Element
    INTEGER :: VvarId, IvarId, q, j, astat
    REAL(KIND=dp), ALLOCATABLE :: Tcoef(:,:,:)
    REAL(KIND=dp) :: RotM(3,3,nn)
    REAL(KIND=dp) :: val, dt, crt(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilType
    LOGICAL :: Found, IsActive

    ASolver => CurrentModel % Asolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('AddComponentEquationsAndCouplings','ASolver not found!')

    Circuit => CurrentModel % Circuits(p)
    
    nm = Asolver % Matrix % NumberOfRows
    CM => CurrentModel%CircuitMatrix

    ALLOCATE(Tcoef(3,3,nn), STAT=astat)
    IF (astat /= 0) THEN
      CALL Fatal('AddComponentEquationsAndCouplings','Memory allocation failed!')
    END IF
    
    DO CompInd = 1, Circuit % n_comp
      Comp => Circuit % Components(CompInd)

      Comp % Resistance = 0._dp 
      Comp % Conductance = 0._dp 

      Cvar => Comp % vvar
      vvarId = Comp % vvar % ValueId + nm
      IvarId = Comp % ivar % ValueId + nm

      CompParams => CurrentModel % Components(Comp % ComponentId) % Values
      IF (.NOT. ASSOCIATED(CompParams)) CALL Fatal ('AddComponentEquationsAndCouplings', 'Component parameters not found')
      IF (Comp % CoilType == 'stranded' .OR. Comp % ComponentType == 'resistor') THEN
        Comp % Resistance = ListGetCReal(CompParams, 'Resistance', Found)
        IF (Found) THEN
          Comp % UseCoilResistance = .TRUE.
        ELSE
          Comp % UseCoilResistance = .FALSE.
        END IF
      END IF

      IsActive = .TRUE.
      IF( Circuit % Parallel ) THEN
        IsActive = ( Cvar % Owner == ParEnv % myPE )
      END IF
      
      IF ( IsActive ) THEN
        IF (Comp % ComponentType == 'resistor') THEN
            CALL Info('AddComponentEquationsAndCouplings',&
                'Writing resistor equation, component '//i2s(CompInd), Level = 7)
            CALL AddToMatrixElement(CM, VvarId, IvarId, Comp % Resistance)
            CALL AddToMatrixElement(CM, VvarId, VvarId, -1._dp)
        ELSE
          SELECT CASE (Comp % CoilType)
          CASE('stranded')
            IF (Comp % UseCoilResistance) THEN
              CALL Info('AddComponentEquationsAndCouplings',&
                  'Using coil resistance for component '//i2s(CompInd), Level = 7)
              CALL AddToMatrixElement(CM, VvarId, IvarId, Comp % Resistance)
            ELSE
              Comp % Resistance = 0._dp
            END IF
            CALL AddToMatrixElement(CM, VvarId, VvarId, -1._dp)
          CASE('massive')
            CALL AddToMatrixElement(CM, VvarId, IvarId, -1._dp)
          CASE('foil winding')
            ! Foil Winding voltage: 
            ! V + ...added next... = 0
            ! ----------------------
            CALL AddToMatrixElement(CM, VvarId, VvarId, 1._dp)
            
            DO j = 1, Cvar % pdofs 
              ! Foil Winding voltage: 
              !  ... - Nf/Lalpha * int_0^{Lalpha}(V_0+V_1*alpha+V_2*alpha**2+...) = 0
              !          => ... - Nf * (V_0*Lalpha^0 + V_1/2*Lalpha^1 + V_2/3*Lalpha^2 + ...) = 0
              ! where V_m is the mth dof of the polynomial
              ! --------------------------------------------------------------
              val = - REAL(Comp % nofturns) / REAL(j) * Comp % coilthickness**(j-1)
              CALL AddToMatrixElement(CM, VvarId, j + VvarId, val)

              ! Circuit eqns for the pdofs:
              ! - Nf/Lalpha * I * int_0^1(Vi'(alpha)) + ...added later... = 0
              ! ----------------------------------------------------------
              CALL AddToMatrixElement(CM, j + VvarId, IvarId, val)
            END DO
          END SELECT
        END IF
      END IF
      
      IF (Comp % ComponentType == 'resistor') CYCLE

      DO q=GetNOFActive(),1,-1
        Element => GetActiveElement(q)
        IF (ElAssocToComp(Element, Comp)) THEN
          CompParams => GetComponentParams( Element )
          IF (.NOT. ASSOCIATED(CompParams)) &
              CALL Fatal ('AddComponentEquationsAndCouplings', 'Component parameters not found')

          CoilType = ListGetString(CompParams, 'Coil Type', UnfoundFatal=.TRUE.)
          
          nn = GetElementNOFNodes(Element)
          nd = GetElementNOFDOFs(Element,ASolver)
          
          IF (SIZE(Tcoef,3) /= nn) THEN
            DEALLOCATE(Tcoef)
            ALLOCATE(Tcoef(3,3,nn), STAT=astat)
            IF ( astat /= 0 ) THEN
              CALL Fatal('AddComponentEquationsAndCouplings', 'Memory allocation error!' )
            END IF
          END IF
          
          Tcoef = GetElectricConductivityTensor(Element, nn, 're', .TRUE., CoilType)
          SELECT CASE(CoilType)
          CASE ('stranded')
            CALL Add_stranded(Element,Tcoef,Comp,nn,nd,dt,CompParams)
          CASE ('massive')
            IF (.NOT. HasSupport(Element,nn)) CYCLE
            CALL Add_massive(Element,Tcoef,Comp,nn,nd,dt,crt)
          CASE ('foil winding')
            IF (.NOT. HasSupport(Element,nn)) CYCLE
            CALL Add_foil_winding(Element,Tcoef,Comp,nn,nd,dt)
          CASE DEFAULT
            CALL Fatal ('AddComponentEquationsAndCouplings', 'Non-existent Coil Type Chosen!')
          END SELECT
        END IF
      END DO
    END DO

    IF( Parallel ) THEN
      DO CompInd = 1, Circuit % n_comp
        Comp => Circuit % Components(CompInd)
        Comp % Resistance = ParallelReduction(Comp % Resistance)
        Comp % Conductance = ParallelReduction(Comp % Conductance)
      END DO
    END IF
      
    DEALLOCATE(Tcoef)
!------------------------------------------------------------------------------
   END SUBROUTINE AddComponentEquationsAndCouplings
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE Add_stranded(Element,Tcoef,Comp,nn,nd,dt,CompParams)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: Tcoef(3,3,nn),dt
    TYPE(Valuelist_t), POINTER :: CompParams
    TYPE(Component_t) :: Comp
    INTEGER :: nn, nd, nm

    INTEGER :: Indexes(nd),VvarId,IvarId
    TYPE(Solver_t), POINTER :: ASolver
    INTEGER, POINTER :: PS(:)
    TYPE(Matrix_t), POINTER :: CM
    TYPE(Nodes_t), SAVE :: Nodes
    REAL(KIND=dp) :: Basis(nd), DetJ, x,POT(nd),pPOT(nd),ppPOT(nd),tscl
    REAL(KIND=dp) :: dBasisdx(nd,3), wBase(nn), w(3)
    REAL(KIND=dp) :: localC, val, circ_eq_coeff, localR !, localL
    INTEGER :: j,t
    LOGICAL :: stat

    TYPE(GaussIntegrationPoints_t) :: IP
    LOGICAL :: CSymmetry, First=.TRUE., InitHandle=.TRUE., &
               CoilUseWvec=.FALSE., CoilUseWvec0=.FALSE.,Found,Found2
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilWVecVarname, CoilType

    REAL(KIND=dp) :: WBasis(nd,3), RotWBasis(nd,3)
    INTEGER :: dim, ncdofs,q
    TYPE(VariableHandle_t), SAVE :: Wvec_h
    TYPE(Variable_t), POINTER, SAVE :: Wpot
    
    LOGICAL :: PiolaVersion = .FALSE.
    
    SAVE CSymmetry, dim, First, InitHandle

    ASolver => CurrentModel % Asolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('Add_stranded','ASolver not found!')

    IF (First) THEN
      First = .FALSE.
      CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )
      dim = CoordinateSystemDimension()

      ! There has been somewhat different philosophies in how to create the scalar and vector fields
      ! that span the current densities. This is an effort to enable the components to all use
      ! the current density computed by the CoilSolver without writing any additional keywords to the
      ! component sections. The idea is that the circuit then only has current densities defined
      ! by the CoilSolver. If this is not desired then also no such keywords should be used in the Solver
      ! section of this module. 
      !------------------------------------------------------------------------------------------------
      CoilUseWvec0 = GetLogical(CurrentModel % Solver % Values, 'Coil Use W Vector', Found2 ) 
      DO i=1,Model % NumberOfComponents
        CoilType = ListGetString(CurrentModel % Components(i) % Values, 'Coil Type',Found)
        IF(.NOT. Found) CYCLE
        IF(CoilType == 'stranded') THEN  ! massive, foil winding
          CoilWVecVarName = GetString(CurrentModel % Components(i) % Values,'W Vector Variable Name', Found)
          IF(Found) EXIT
        END IF
      END DO
      IF(.NOT. Found) THEN
        CoilWVecVarName = GetString(CurrentModel % Solver % Values,'W Vector Variable Name', Found)
        IF(.NOT. Found) THEN
          IF( GetLogical(CurrentModel % Solver % Values,'Use Nodal CoilCurrent',Found ) ) &
              CoilWVecVarname = 'CoilCurrent'
        END IF
         IF(.NOT. Found) THEN
          IF( GetLogical(CurrentModel % Solver % Values,'Use Elemental CoilCurrent',Found ) ) &
              CoilWVecVarname = 'CoilCurrent e'
        END IF
        IF(Found) CALL Info('Add_stranded','Setting coil current to: '//TRIM(CoilWVecVarname),Level=6)
        ! If we did not find w vector named in any component it is fair to assume that it is globally used!
        IF(.NOT. Found2) CoilUseWvec0 = Found 
      END IF
      IF(.NOT. Found) CoilWVecVarname = 'W Vector E'
      CALL ListInitElementVariable(Wvec_h, CoilWVecVarname)

      PiolaVersion = GetLogical( ASolver % Values, 'Use Piola Transform', Found )

      CALL GetWPotentialVar(Wpot)
    END IF

    PS => Asolver % Variable % Perm

    CM => CurrentModel % CircuitMatrix
    nm = CurrentModel % Asolver % Matrix % NumberOfRows
    
    CALL GetElementNodes(Nodes)
    nd = GetElementDOFs(Indexes,Element,ASolver)

    CALL GetLocalSolution(pPOT,UElement=Element,USolver=ASolver,tstep=-1)

    IF(Solver % Order<2.OR.GetTimeStep()<=2) THEN 
      tscl=1.0_dp
    ELSE
      tscl=1.5_dp
      CALL GetLocalSolution(ppPOT,UElement=Element,USolver=ASolver,tstep=-2)
      pPot = 2*pPOT - 0.5_dp*ppPOT
    END IF
    
    ncdofs=nd
    IF (dim == 3) THEN
      ! We can choose the base per component.
      CoilUseWvec = GetLogical(CompParams, 'Coil Use W Vector', Found)
      IF (.NOT. Found) CoilUseWvec = CoilUseWvec0
    
      IF (.NOT. CoilUseWvec) THEN
        !CALL GetLocalSolution(Wbase, 'w')
        ! when W Potential solver is used, 'w' is not enough.
        CALL GetLocalSolution( Wbase,UElement=Element,UVariable=Wpot, Found=Found)
      END IF

      ncdofs=nd-nn
    END IF

    VvarId = Comp % vvar % ValueId + nm
    IvarId = Comp % ivar % ValueId + nm

    ! Numerical integration:
    ! ----------------------
    IF(PiolaVersion) THEN
      IP = GaussPoints(Element, EdgeBasis=.TRUE., EdgeBasisDegree=1, PReferenceElement=PiolaVersion)
    ELSE
      IP = GaussPoints(Element)
    END IF
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      
      circ_eq_coeff = 1._dp
      SELECT CASE(dim)
      CASE(2)
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                  IP % W(t), detJ, Basis,dBasisdx )

        w = [0._dp, 0._dp, 1._dp]
        IF( CSymmetry ) THEN
          x = SUM( Basis(1:nn) * Nodes % x(1:nn) )
          detJ = detJ * x
        END IF
        circ_eq_coeff = GetCircuitModelDepth()
      CASE(3)

        IF (PiolaVersion) THEN
          stat = EdgeElementInfo( Element, Nodes, IP % U(t), IP % V(t), IP % W(t), &
               DetF = DetJ, Basis = Basis, dBasisdx=dBasisdx, EdgeBasis = WBasis, RotBasis = RotWBasis, &
               BasisDegree=1, ApplyPiolaTransform = .TRUE.)
        ELSE
          stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                    IP % W(t), detJ, Basis, dBasisdx )
          CALL GetEdgeBasis(Element,WBasis,RotWBasis,Basis,dBasisdx)
        END IF

        IF (CoilUseWvec) THEN
          w = ListGetElementVectorSolution( Wvec_h, Basis, Element, dofs = dim )
        ELSE
          w = -MATMUL(WBase(1:nn), dBasisdx(1:nn,:))
        END IF
      END SELECT

      localC = SUM(Tcoef(1,1,1:nn) * Basis(1:nn))
      
      IF (.NOT. Comp % UseCoilResistance) THEN
        ! I * R, where 
        ! R = (1/sigma * js,js):
        ! ----------------------
        localR = Comp % N_j **2 * IP % s(t)*detJ*SUM(w*w)/localC*circ_eq_coeff / Comp % VoltageFactor
        Comp % Resistance = Comp % Resistance + localR
      
        CALL AddToMatrixElement(CM, VvarId, IvarId, localR)
      END IF
      
      DO j=1,ncdofs
        q=j
        IF (dim == 3) q=q+nn
        IF (Comp % N_j/=0._dp) THEN
          ! ( d/dt a,w )        

          IF ( TransientSimulation ) THEN 
            IF (dim == 2) val = Comp % N_j * IP % s(t)*detJ*Basis(j)*circ_eq_coeff/dt*w(3)
            IF (dim == 3) val = Comp % N_j * IP % s(t)*detJ*SUM(WBasis(j,:)*w)/dt
            val = val / Comp % VoltageFactor

!           localL = val
!           Comp % Inductance = Comp % Inductance + localL

            CALL AddToMatrixElement(CM, VvarId, PS(Indexes(q)), tscl * val)
            CM % RHS(vvarid) = CM % RHS(vvarid) + pPOT(q) * val
         END IF
          
          ! source: 
          ! (J, rot a'), where
          ! J = w*I, thus I*(w, rot a'):
          ! ----------------------------         
          IF (dim == 2) val = -Comp % N_j*IP % s(t)*detJ*Basis(j)*w(3)
          IF (dim == 3) val = -Comp % N_j*IP % s(t)*detJ*SUM(WBasis(j,:)*w)

          val = val * Comp % SymmetryCoeff

          CALL AddToMatrixElement(CM,PS(Indexes(q)), IvarId, val)
        END IF
      END DO
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE Add_stranded
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE Add_massive(Element,Tcoef,Comp,nn,nd,dt,crt)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t) :: Element
    REAL(KIND=dp) :: Tcoef(3,3,nn),dt, crt(:)
    TYPE(Component_t) :: Comp

    TYPE(Solver_t), POINTER :: ASolver
    INTEGER, POINTER :: PS(:)
    TYPE(Matrix_t), POINTER :: CM
    REAL(KIND=dp) :: Basis(nd), DetJ, x,POT(nd),pPOT(nd),ppPOT(nd),pVel(nd),pAcc(nd),pAcc2(nd),  &
                    tscl, ascl, bscl, alpha, beta, gamma, delta, prevV, param_a, param_b
    REAL(KIND=dp) :: dBasisdx(nd,3)
    REAL(KIND=dp) :: LondonLambda(nn)
    REAL(KIND=dp) :: LondonLambda_ip, Permittivity(nn), LocalP
    REAL(KIND=dp) :: localC, val, circ_eq_coeff, grads_coeff, localConductance !, localL
    INTEGER :: nn, nd, j, t, nm, Indexes(nd), &
               VvarId, dim
    LOGICAL :: stat
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    LOGICAL :: CSymmetry, First=.TRUE.
    LOGICAL :: LondonEquations
    TYPE(ValueList_t), POINTER :: Material

    REAL(KIND=dp) :: wBase(nn), gradv(3), WBasis(nd,3), RotWBasis(nd,3)
    INTEGER :: ncdofs,q
    TYPE(Variable_t), POINTER, SAVE :: Wpot
    
    SAVE CSymmetry, dim, First

    IF (First) THEN
      First = .FALSE.
      CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )
      dim = CoordinateSystemDimension()

      CALL GetWPotentialVar(Wpot)
    END IF

    ASolver => CurrentModel % Asolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('Add_massive','ASolver not found!')
    PS => Asolver % Variable % Perm

    CM => CurrentModel % CircuitMatrix
    nm = CurrentModel % Asolver % Matrix % NumberOfRows

    vvarId = Comp % vvar % ValueId + nm

    Material => GetMaterial( Element )

    CALL GetElementNodes(Nodes)
    nd = GetElementDOFs(Indexes,Element,ASolver)
    IF (ASolver % TimeOrder==2) THEN
      CALL GetLocalSolution(pPot,UElement=Element,USolver=ASolver,tstep=-3)
      CALL GetLocalSolution(pVel,UElement=Element,USolver=ASolver,tstep=-4)
      CALL GetLocalSolution(pAcc,UElement=Element,USolver=ASolver,tstep=-5)
      CALL GetLocalSolution(pAcc2,UElement=Element,USolver=ASolver,tstep=-7)

      param_a = Asolver % alpha
      param_b = Asolver % beta
      IF(param_b>=0) THEN ! this is rho_inf for Generalized-alpha
        alpha  = (2*param_b-1)/(param_b+1)
        beta = 1/(param_b+1)**2
        gamma = (3-param_b)/(2*(param_b+1))
        delta = param_b/(param_b+1)
       ELSE !bossak
        alpha = param_a
        beta = (1-alpha)**2/4
        gamma = 0.5d0 - alpha
        delta = 0
        pAcc2(1:nd) = pAcc(1:nd)
      END IF
      prevV = Crt(vvarId-nm)
      CALL GetPermittivity(Material, Permittivity, nn)
    ELSE
      CALL GetLocalSolution(pPOT,UElement=Element,USolver=ASolver,tstep=-1)
      IF(Solver % Order<2.OR.GetTimeStep()<=2) THEN 
        tscl=1.0_dp
      ELSE
        tscl=1.5_dp
        CALL GetLocalSolution(ppPOT,UElement=Element,USolver=ASolver,tstep=-2)
        pPot = (2*pPOT - 0.5_dp*ppPOT)
      END IF
    END IF

    ncdofs=nd
    IF (dim == 3) THEN
      !CALL GetLocalSolution(Wbase, 'w')      
      CALL GetLocalSolution( Wbase,UElement=Element,UVariable=Wpot, Found=Found)
      ncdofs=nd-nn
    END IF

    LondonLambda(:) = GetReal( Material, 'London Lambda', LondonEquations, Element)

    ! Numerical integration:
    ! ----------------------
    IP = GaussPoints(Element)
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                IP % W(t), detJ, Basis,dBasisdx )
      
      grads_coeff = -1._dp
      circ_eq_coeff = 1._dp
      SELECT CASE(dim)
      CASE(2)
        IF( CSymmetry ) THEN
          x = SUM( Basis(1:nn) * Nodes % x(1:nn) )
          detJ = detJ * x
          grads_coeff = grads_coeff/x
        END IF
        circ_eq_coeff = GetCircuitModelDepth()
        grads_coeff = grads_coeff/circ_eq_coeff
      CASE(3)
        CALL GetEdgeBasis(Element,WBasis,RotWBasis,Basis,dBasisdx)
        gradv = MATMUL( WBase(1:nn), dBasisdx(1:nn,:))
      END SELECT

      localC = SUM(Tcoef(1,1,1:nn) * Basis(1:nn))
      localP = SUM(Permittivity(1:nn) * Basis(1:nn))

      ! computing the source term Vi(sigma grad v0, grad si):
      ! ------------------------------------------------
      IF(dim==2) val = IP % s(t)*detJ*localC*grads_coeff**2*circ_eq_coeff
      IF(dim==3) val = IP % s(t)*detJ*localC*SUM(gradv*gradv)
      val = val * Comp % VoltageFactor

      localConductance = ABS(val)
      Comp % Conductance = Comp % Conductance + localConductance

      CALL AddToMatrixElement(CM, vvarId, vvarId, val)

      IF(Asolver % TimeOrder==2) THEN
        IF(dim==2) val = IP % s(t)*detJ*localP*grads_coeff**2*circ_eq_coeff
        IF(dim==3) val = IP % s(t)*detJ*localP*SUM(gradv*gradv)

        CALL AddToMatrixElement(CM, vvarId, vvarId, val/dt)
        CM % RHS(vvarId) = CM % RHS(vvarId) + val*prevV/dt
      END IF

      IF ( LondonEquations ) THEN
        LondonLambda_ip = SUM( Basis(1:nn) * LondonLambda(1:nn) )

        IF(dim==2) val = IP % s(t)*detJ/LondonLambda_ip*grads_coeff**2*circ_eq_coeff
        val = val * Comp % VoltageFactor
        ! Phi (beta grad phi_0, grad phi')
        ! Here Phi takes the place of Vi
        ! -----------------------------------
        CALL AddToMatrixElement(CM, vvarId, vvarId, val)
      END IF

      DO j=1,ncdofs
        q=j
        IF (dim == 3) q=q+nn

        IF ( LondonEquations ) THEN
          ! Phi * ( beta * grad phi, a')
          ! where phi is the node flux scalar potential
          ! -------------------------------------------
          IF(dim==2) val = IP % s(t)*detJ/LondonLambda_ip*basis(j)*grads_coeff*circ_eq_coeff
          CALL AddToMatrixElement(CM, vvarId, PS(Indexes(q)), val)

          IF(dim==2) val = IP % s(t)*detJ/LondonLambda_ip*basis(j)*grads_coeff
          val = val * Comp % VoltageFactor
          CALL AddToMatrixElement(CM, PS(indexes(q)), vvarId, val)
        END IF

        ! computing the mass term (sigma d/dt a, grad si):
        ! ---------------------------------------------------------
        IF ( TransientSimulation ) THEN 
          IF(dim==2) val = IP % s(t)*detJ*basis(j)*grads_coeff*circ_eq_coeff
          IF(dim==3) val = IP % s(t)*detJ*SUM(Wbasis(j,:)*gradv)

          IF(Asolver % Timeorder==1) THEN
            CALL AddToMatrixElement(CM, vvarId, PS(Indexes(q)), tscl * val * localC/dt)
            CM % RHS(vvarid) = CM % RHS(vvarid) + pPOT(q) * val * localC/dt
          ELSE
            ! 1st time derivative
            CALL AddToMatrixElement(CM, vvarId, PS(Indexes(q)), &
                  val * localC * gamma/(beta*dt) )
            CM % RHS(vvarid) = CM % RHS(vvarid) + val * localC*( gamma/(beta*dt)*pPot(q) + &
               (gamma/beta-1)*pVel(q) - ((1-gamma)+gamma*(1-1/(2*beta)))*dt*pAcc(q) )

            ! 2nd time derivative
            CALL AddToMatrixElement(CM, vvarId, PS(Indexes(q)), &
                 val * localP * (1-alpha)/(1-delta)/(beta*dt**2) )
            CM % RHS(vvarid) = CM % RHS(vvarid) + val * localP*( &
                (1-alpha)/(1-delta)/(beta*dt**2)*pPot(q) + &
                  (1-alpha)/(1-delta)/(beta*dt)*pVel(q) + &
                    delta/(1-delta)*pAcc(q) - &
                      ((1-alpha)*(1-1/(2*beta))+alpha)/(1-delta)*pAcc2(q) )

            IF(dim==2) val = IP % s(t)*detJ*localC*basis(j)*grads_coeff
            IF(dim==3) val = IP % s(t)*detJ*localC*SUM(gradv*Wbasis(j,:))
            val = val * Comp % VoltageFactor
            CALL AddToMatrixElement(CM, PS(indexes(q)), vvarId, val /  dt)
            CM % RHS(PS(indexes(q))) = CM % RHS(PS(indexes(q))) + prevV * val / dt
          END IF
        END IF

        IF(dim==2) val = IP % s(t)*detJ*localC*basis(j)*grads_coeff
        IF(dim==3) val = IP % s(t)*detJ*localC*SUM(gradv*Wbasis(j,:))
        val = val * Comp % VoltageFactor
        CALL AddToMatrixElement(CM, PS(indexes(q)), vvarId, val)
      END DO
    END DO

!------------------------------------------------------------------------------
   END SUBROUTINE Add_massive
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE Add_foil_winding(Element,Tcoef,Comp,nn,nd,dt)
!------------------------------------------------------------------------------
    USE MGDynMaterialUtils
    IMPLICIT NONE
    INTEGER :: nn, nd
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: Tcoef(3,3,nn), C(3,3), val, dt
    TYPE(Component_t) :: Comp

    TYPE(Solver_t), POINTER :: ASolver
    INTEGER, POINTER :: PS(:)
    TYPE(Matrix_t), POINTER :: CM
    REAL(KIND=dp) :: Basis(nd), DetJ, localAlpha, localV, localVtest, &
                     x, circ_eq_coeff, grads_coeff,POT(nd),pPOT(nd),ppPOT(nd),tscl
    REAL(KIND=dp) :: dBasisdx(nd,3),alpha(nn)
    REAL(KIND=dp) :: localR !, localL
    INTEGER :: nm,p,j,t,Indexes(nd),vvarId,vpolord_tot, &
               vpolord, vpolordtest, dofId, dofIdtest, &
               dim
    LOGICAL :: stat
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    LOGICAL :: CSymmetry, First=.TRUE.

    REAL(KIND=dp) :: wBase(nn), gradv(3), WBasis(nd,3), RotWBasis(nd,3), &
                     RotMLoc(3,3), RotM(3,3,nn)
    INTEGER :: i,ncdofs,q
    TYPE(Variable_t), POINTER, SAVE :: Wpot
    
    SAVE CSymmetry, dim, First

    IF (First) THEN
      First = .FALSE.
      CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )
      dim = CoordinateSystemDimension()

      CALL GetWPotentialVar(Wpot)
    END IF

    ASolver => CurrentModel % Asolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('Add_foil_winding','ASolver not found!')
    PS => Asolver % Variable % Perm

    CM => CurrentModel % CircuitMatrix
    nm = CurrentModel % Asolver % Matrix % NumberOfRows

    CALL GetElementNodes(Nodes)
    nd = GetElementDOFs(Indexes,Element,ASolver)
    CALL GetLocalSolution(alpha,'Alpha')

    CALL GetLocalSolution(pPOT,UElement=Element,USolver=ASolver,tstep=-1)

    IF(Solver % Order<2.OR.GetTimeStep()<=2) THEN 
      tscl=1.0_dp
    ELSE
      tscl=1.5_dp
      CALL GetLocalSolution(ppPOT,UElement=Element,USolver=ASolver,tstep=-2)
      pPot = 2*pPOT - 0.5_dp*ppPOT
    END IF

    ncdofs=nd
    IF (dim == 3) THEN
      CALL GetLocalSolution( Wbase,UElement=Element,UVariable=Wpot, Found=Found)
      !CALL GetLocalSolution(Wbase, 'w')
      CALL GetElementRotM(Element, RotM, nn)
      ncdofs=nd-nn
    END IF

    vvarId = Comp % vvar % ValueId
    vpolord_tot = Comp % vvar % pdofs - 1

    ! Numerical integration:
    ! ----------------------
    IP = GaussPoints(Element)
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                IP % W(t), detJ, Basis,dBasisdx )

      grads_coeff = -1._dp
      circ_eq_coeff = 1._dp
      SELECT CASE(dim)
      CASE(2)
        IF( CSymmetry ) THEN
          x = SUM( Basis(1:nn) * Nodes % x(1:nn) )
          detJ = detJ * x
          grads_coeff = grads_coeff/x
        END IF
        circ_eq_coeff = GetCircuitModelDepth()
        grads_coeff = grads_coeff/circ_eq_coeff
        C(1,1) = SUM( Tcoef(3,3,1:nn) * Basis(1:nn) )
        ! I * R, where 
        ! R = (1/sigma * js,js):
        ! ----------------------
        localR = Comp % N_j **2 * IP % s(t)*detJ/C(1,1)*circ_eq_coeff/Comp % VoltageFactor
      CASE(3)
        CALL GetEdgeBasis(Element,WBasis,RotWBasis,Basis,dBasisdx)
        gradv = MATMUL( WBase(1:nn), dBasisdx(1:nn,:))
        ! Compute the conductivity tensor
        ! -------------------------------
        DO i=1,3
          DO j=1,3
            C(i,j) = SUM( Tcoef(i,j,1:nn) * Basis(1:nn) )
            RotMLoc(i,j) = SUM( RotM(i,j,1:nn) * Basis(1:nn) )
          END DO
        END DO

        ! I * R, where 
        ! R = (1/sigma * js,js):
        ! ----------------------
        localR = Comp % N_j **2 * IP % s(t)*detJ/C(3,3)/Comp % VoltageFactor
        ! Transform the conductivity tensor:
        ! ----------------------------------
        C = MATMUL(MATMUL(RotMLoc, C),TRANSPOSE(RotMLoc))
      END SELECT
      
      localAlpha = SUM(alpha(1:nn) * Basis(1:nn))
      
      ! alpha is normalized to be in [0,1] thus, 
      ! it needs to be multiplied by the thickness of the coil 
      ! to get the real alpha:
      ! ------------------------------------------------------
      localAlpha = localAlpha * Comp % coilthickness

      Comp % Resistance = Comp % Resistance + localR

      DO vpolordtest=0,vpolord_tot ! V'(alpha)
        localVtest = localAlpha**vpolordtest
        dofIdtest = vpolordtest + 1 + vvarId
        DO vpolord = 0, vpolord_tot ! V(alpha)

          localV = localAlpha**vpolord
          dofId = vpolord + 1 + vvarId
          
          ! Computing the stiff term (sigma V(alpha) grad v0, V'(alpha) grad si):
          ! ---------------------------------------------------------------------
          IF (dim == 2) val = IP % s(t)*detJ*localV*localVtest*C(1,1)*grads_coeff**2*circ_eq_coeff
          IF (dim == 3) val = IP % s(t)*detJ*localV*localVtest*SUM(MATMUL(C,gradv)*gradv)
          val = val * Comp % VoltageFactor
          CALL AddToMatrixElement(CM, dofIdtest+nm, dofId+nm, val)
        END DO


        IF ( TransientSimulation ) THEN 
          DO j=1,ncdofs
            q=j
            IF (dim == 3) q=q+nn
            ! computing the mass term (sigma * d/dt * a, V'(alpha) grad si):
            ! ---------------------------------------------------------
            IF (dim == 2) val = IP % s(t)*detJ*localVtest*C(1,1)*basis(j)*grads_coeff*circ_eq_coeff/dt
            IF (dim == 3) val = IP % s(t)*detJ*localVtest*SUM(MATMUL(C,Wbasis(j,:))*gradv)/dt
  !          localL = val
  !          Comp % Inductance = Comp % Inductance + localL
            CALL AddToMatrixElement(CM, dofIdtest+nm, PS(Indexes(q)), tscl * val)
            CM % RHS(dofIdtest+nm) = CM % RHS(dofIdtest+nm) + pPOT(q) * val
          END DO
        END IF

      END DO

      DO vpolord = 0, vpolord_tot ! V(alpha)
        localV = localAlpha**vpolord
        dofId = vpolord + 1 + vvarId

        DO j=1,ncdofs
            q=j
            IF (dim == 3) q=q+nn
            IF (dim == 2) val = IP % s(t)*detJ*localV*C(1,1)*basis(j)*grads_coeff
            IF (dim == 3) val = IP % s(t)*detJ*localV*SUM(MATMUL(C,gradv)*Wbasis(j,:))
            val = val * Comp % VoltageFactor
            CALL AddToMatrixElement(CM, PS(indexes(q)), dofId+nm, val)
        END DO
      END DO

    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE Add_foil_winding
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE GetConductivity(Element, Tcoef, nn)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER :: Element
    TYPE(Valuelist_t), POINTER :: Material
    REAL(KIND=dp) :: Tcoef(3,3,nn)
    REAL(KIND=dp), POINTER, SAVE :: Cwrk(:,:,:)
    INTEGER :: nn, i, j
    LOGICAL, SAVE :: visited = .FALSE., Found

    IF (.NOT. visited) THEN
      NULLIFY( Cwrk )
    END IF

    Tcoef = 0.0_dp
    Material => GetMaterial( Element )
    IF (.NOT. ASSOCIATED(Material)) CALL Fatal('Circuits_apply','Material not found.')

    CALL ListGetRealArray( Material, &
           'Electric Conductivity', Cwrk, nn, Element % NodeIndexes, Found )

    IF (.NOT. Found) CALL Fatal('GetConductivity', 'Electric Conductivity not found.')
    
    IF (Found) THEN
       IF ( SIZE(Cwrk,1) == 1 ) THEN
          DO i=1,3
             Tcoef( i,i,1:nn ) = Cwrk( 1,1,1:nn )
          END DO
       ELSE IF ( SIZE(Cwrk,2) == 1 ) THEN
          DO i=1,MIN(3,SIZE(Cwrk,1))
             Tcoef(i,i,1:nn) = Cwrk(i,1,1:nn)
          END DO
       ELSE
          DO i=1,MIN(3,SIZE(Cwrk,1))
             DO j=1,MIN(3,SIZE(Cwrk,2))
                Tcoef( i,j,1:nn ) = Cwrk(i,j,1:nn)
             END DO
          END DO
       END IF
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE GetConductivity
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Set the rotation angle in case moment of inertia and torque are given.
!------------------------------------------------------------------------------
  SUBROUTINE SetDynamicAngle()
    TYPE(Variable_t), POINTER :: AngVar, VeloVar
    TYPE(ValueList_t), POINTER :: Simulation
    REAL(KIND=dp) :: dt, ang, velo, ang0, velo0, imom, torq    
    INTEGER :: tStep, tStepPrev = 0
    LOGICAL :: Found
    
    SAVE ang0, velo0, imom, tStepPrev
    
    ! Variable should already exist as it was introduced in the _init section.
    AngVar => DefaultVariableGet( 'Rotor Angle' )
    IF(.NOT. ASSOCIATED( AngVar ) ) RETURN
    
    VeloVar => DefaultVariableGet( 'Rotor Velo' )
    IF(.NOT. ASSOCIATED( VeloVar ) ) THEN
      CALL Fatal('SetRotation','Variable > Rotor Velo < does not exist!')
    END IF
    
    Simulation => GetSimulation()

    IF( ListCheckPresent( Model % Simulation,'Rotor Angle') ) THEN
      CALL Info('SetRotation','Using "Rotor Velo" from simulation section',Level=6)
    ELSE      
      CALL Info('SetRotation','Using computed torque to set rotation!',Level=6)
      tStep = GetTimestep()

      imom = GetConstReal( Simulation,'Imom',Found) ! interatial moment of the motor      
      IF(.NOT. Found ) THEN
        CALL Info('SetDynamicAngle','Moment of inertia "Imom" not giving, skipping dynamics!',Level=7)
        RETURN
      END IF
      
      IF(imom < EPSILON(imom) ) THEN
        CALL Info('SetDynamicAngle','Moment of inertia "Imom" close to zero, skipping dynamics!',Level=7)
        RETURN
      END IF
      
      ! We initiate these at the start of the timestep when they still present the previous
      ! computed values. 
      IF( tStep /= tStepPrev ) THEN
        ang0 = AngVar % Values(1)
        velo0 = VeloVar % Values(1)
        tStepPrev = tStep
      END IF

      torq = GetConstReal( Simulation,'res: Air Gap Torque', Found)
      IF(.NOT. Found ) THEN
        CALL Warn('SetRotation','Without torque rotor velocity stays the same!')
      ELSE
        velo = velo0 + dt * (torq-0) / imom
        ang  = ang0  + dt * velo
      END IF
      VeloVar % Values(1) = velo
      AngVar % Values(1) = ang
    END IF
      
    CALL ListAddConstReal(Simulation,'res: Angle(rad)', ang)
    CALL ListAddConstReal(Simulation,'res: Speed(rpm)', velo/(2._dp*pi)*60)
      
  END SUBROUTINE SetDynamicAngle
  

!------------------------------------------------------------------------------
END SUBROUTINE CircuitsAndDynamics
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Initialization for the primary solver: CurrentSource
!------------------------------------------------------------------------------
SUBROUTINE CircuitsAndDynamicsHarmonic_init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: TransientSimulation !< Steady state or transient simulation
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params
  
  Params => Solver % Values

  ! This is only created if no variable present!
  CALL ListAddNewString( Params,'Variable','-global cmplxckt')
  CALL ListAddNewLogical( Params,'Variable Output',.FALSE.)
  CALL ListAddNewLogical( Params,'No Matrix',.TRUE.)

  Solver % Values => Params
    
!------------------------------------------------------------------------------
END SUBROUTINE CircuitsAndDynamicsHarmonic_init
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
SUBROUTINE CircuitsAndDynamicsHarmonic( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE CircuitUtils
  USE CircMatInitMod
  USE MGDynMaterialUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: TransientSimulation !< Steady state or transient simulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  LOGICAL :: First=.TRUE.
  TYPE(Solver_t), POINTER :: Asolver => Null()
  INTEGER :: p, n, istat, max_element_dofs, i, j
  TYPE(Mesh_t), POINTER :: Mesh  
  TYPE(Matrix_t), POINTER :: CM
  INTEGER, POINTER :: n_Circuits => Null(), circuit_tot_n => Null()
  TYPE(Circuit_t), POINTER :: Circuits(:)  
  LOGICAL :: Parallel, Found, EigenSystem
  REAL(KIND=dp), POINTER :: px(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: sname
  CHARACTER(*), PARAMETER :: Caller = 'CircuitsAndDynamicsHarmonic'
  
  SAVE First, Parallel
  
!------------------------------------------------------------------------------


  CALL DefaultStart()

  IF (First) THEN
    CALL Info(Caller,'Initializing electric circuits for harmonic simulation',Level=6)

    First = .FALSE.

    Parallel = Solver % Parallel 
    IF(Parallel) THEN
      CALL Info(Caller,'Assuming parallel electric circuits',Level=12)
    ELSE
      CALL Info(Caller,'Assuming serial electric circuits',Level=12)
    END IF
    
    Model % HarmonicCircuits = .TRUE.
    CALL AddComponentsToBodyLists()
    
    ALLOCATE( Model % Circuit_tot_n, Model % n_Circuits, STAT=istat )
    IF ( istat /= 0 ) THEN
      CALL Fatal( Caller, 'Memory allocation error.' )
    END IF

    n_Circuits => Model % n_Circuits
    Model % Circuit_tot_n = 0

    ! Look for the solver we attach the circuit equations to:
    ! -------------------------------------------------------
    Asolver => NULL()
    DO i=1,Model % NumberOfSolvers      
      sname = GetString(Model % Solvers(i) % Values, 'Procedure', Found)
      j = INDEX( sname,'MagnetoDynamics2DHarmonic')
      IF(j==0) j = INDEX( sname,'WhitneyAVHarmonicSolver')
      IF( j > 0 ) THEN
        ASolver => Model % Solvers(i) 
        EXIT
      END IF
    END DO
    
    IF(.NOT. ASSOCIATED(ASolver) ) THEN
      ASolver => FindSolverWithKey('Export Lagrange Multiplier')
    END IF
    CALL Info(Caller,'Circuit equations associated with solver index: '&
        //I2S(ASolver % SolverId),Level=6)
    Model % ASolver => ASolver 

    IF( Solver % Parallel .NEQV. Asolver % Parallel  ) THEN
      CALL Warn(Caller,'Conflicting parallel status for circuit and A solver!')
      Solver % Parallel = .TRUE.
      ASolver % Parallel = .TRUE.
    END IF
    
    CALL AllocateCircuitsList() ! CurrentModel%Circuits
    Circuits => Model % Circuits

    CALL SetBoundaryAreasToValueLists() 
    
    DO p=1,n_Circuits
      
      n = GetNofCircVariables(p)
      CALL AllocateCircuit(p)
      
      Circuits(p) % n_comp = CountNofCircComponents(p, n)
      ALLOCATE(Circuits(p) % Components(Circuits(p) % n_comp))
      
      Circuits(p) % Harmonic = .TRUE.
      Circuits(p) % Parallel = Parallel 
      
      CALL ReadCircuitVariables(p)
      CALL ReadComponents(p)
      CALL AddComponentValuesToLists(p)  ! Lists are used to communicate values to other solvers at the moment...
      CALL AddBareCircuitVariables(p)   ! these don't belong to any components
      CALL ReadCoefficientMatrices(p)
      CALL ReadPermutationVector(p)
      CALL ReadCircuitSources(p)
      CALL WriteCoeffVectorsForCircVariables(p)
    
    END DO

    ! Create CRS matrix structures for the circuit equations:
    ! ------------------------------------------------------
    CALL Circuits_MatrixInit()
  END IF
  
  EigenSystem = GetLogical( Asolver % Values, 'Eigen Analysis', Found )

  max_element_dofs = Model % Mesh % MaxElementDOFs
  Circuits => Model % Circuits
  n_Circuits => Model % n_Circuits
  CM => Model % CircuitMatrix
  
  ! Initialize Circuit matrix:
  ! -----------------------------
  IF(.NOT.ASSOCIATED(CM)) RETURN

  IF (SIZE(CM % values) <= 0) RETURN
  CM % RHS = 0._dp
  IF(ASSOCIATED(CM % Values)) CM % Values = 0._dp

  ! Write Circuit equations:
  ! ------------------------
  DO p = 1,n_Circuits
    CALL AddBasicCircuitEquations(p)
    CALL AddComponentEquationsAndCouplings(p, max_element_dofs)
  END DO
  Asolver % Matrix % AddMatrix => CM

  IF(ASSOCIATED(CM)) THEN
    IF(  CM % Format == MATRIX_LIST ) CALL List_toCRSMatrix(CM)
    IF(CM % NumberOfRows<=0)  THEN
      CALL FreeMatrix(CM)
      Asolver % Matrix % AddMatrix => Null()
    END IF
  ELSE
     ASolver % Matrix % AddMatrix => Null()
  END IF

  IF( InfoActive(20) ) THEN
    px => CM % Values
    CALL VectorValuesRange(px,SIZE(px),'CircuitMatrix',.TRUE.)
    px => CM % rhs
    CALL VectorValuesRange(px,SIZE(px),'CircuitRhs',.TRUE.)
  END IF

  IF( ListGetLogical( Solver % Values,'Save Circuit Matrix',Found ) ) THEN
    CALL SaveLinearSystem( Solver, CM,'circuit',ASolver % Matrix % NumberOfRows)
  END IF
    
  CALL DefaultFinish()

  CALL Info(Caller,'Finished assembly of circuit matrix',Level=12)

  
  CONTAINS

!------------------------------------------------------------------------------
   SUBROUTINE AddBasicCircuitEquations(p)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(CircuitVariable_t), POINTER :: Cvar
    TYPE(ValueList_t), POINTER :: BF, Params
    TYPE(Matrix_t), POINTER :: CM
    INTEGER :: p, i, nm, RowId, ColId, j
    REAL(KIND=dp) :: Omega, vphi
    COMPLEX(KIND=dp) :: cmplx_val
    LOGICAL :: Found
    COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
    
    Circuit => CurrentModel % Circuits(p)
    nm = CurrentModel % Asolver % Matrix % NumberOfRows
    BF => CurrentModel % BodyForces(1) % Values
    CM => CurrentModel%CircuitMatrix

    Omega = GetAngularFrequency()
    WRITE(Message,'(A,ES12.3)') 'Angular frequency for circuit equations: ',Omega
    CALL Info(Caller, Message, Level=6) 
    
    Params => GetSolverParams()
    
    DO i=1,Circuit % n
      Cvar => Circuit % CircuitVariables(i)

      IF(Circuit % Parallel ) THEN
        IF(Cvar % Owner /= ParEnv % myPE) CYCLE
      END IF
        
      RowId = Cvar % ValueId + nm
      
      vphi = GetCReal(Params, TRIM(Circuit % Source(i))//" re", Found)
      IF ( .NOT.Found.AND.ASSOCIATED(BF) ) THEN
        vphi = GetCReal(BF, TRIM(Circuit % Source(i))//" re", Found)
      END IF
      IF (Found) Cvar % SourceRe(i) = vphi

      !IF(Found) PRINT *,'vphi re',Found,i,vphi,TRIM(Circuit % Source(i))
       
      vphi = GetCReal(Params, TRIM(Circuit % Source(i))//" im", Found)
      IF ( .NOT.Found.AND.ASSOCIATED(BF) ) THEN
        vphi = GetCReal(BF, TRIM(Circuit % Source(i))//" im", Found)
      END IF
      IF (Found) Cvar % SourceIm(i) = vphi
      
      !IF(Found) PRINT *,'vphi im',i,vphi,TRIM(Circuit % Source(i))

      CM % RHS(RowId) = Cvar % SourceRe(i)
      CM % RHS(RowId+1) = Cvar % SourceIm(i)
        
      DO j=1,Circuit % n

        ColId = Circuit % CircuitVariables(j) % ValueId + nm

        ! im * Omega * A x: (x could be voltage or current):
        !--------------------------------------------
        IF(Cvar % A(j) /= 0._dp) THEN
          CALL AddMatrixEntry( CM, RowId, ColId, (0._dp,0._dp), im*Omega*Cvar % A(j) )
        END IF

        ! B x:
        ! ------
        IF(Cvar % B(j) /= 0._dp) THEN
          IF (Cvar % Mre(j) /= 0._dp .OR. Cvar % Mim(j) /= 0._dp) THEN
            cmplx_val = Cvar % Mre(j) + im * Cvar % Mim(j)
            cmplx_val = cmplx_val * Cvar % B(j)
          ELSE
            cmplx_val = Cvar % B(j)
          END IF
          
          CALL AddToCmplxMatrixElement(CM, RowId, ColId, REAL(cmplx_val), AIMAG(cmplx_val))
        END IF
      END DO
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE AddBasicCircuitEquations
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE AddComponentEquationsAndCouplings(p, nn)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: p, nn
    INTEGER :: CompInd, nm
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(Matrix_t), POINTER :: CM
    TYPE(Component_t), POINTER :: Comp
    TYPE(CircuitVariable_t), POINTER :: Cvar
    TYPE(Valuelist_t), POINTER :: CompParams
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp), ALLOCATABLE :: sigma_33(:), sigmaim_33(:)
    COMPLEX(KIND=dp), ALLOCATABLE :: Tcoef(:,:,:)
    INTEGER :: VvarId, IvarId, q, j, astat
    COMPLEX(KIND=dp) :: i_multiplier, cmplx_val
    COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
    REAL(KIND=dp) :: RotM(3,3,nn)
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilType
    LOGICAL :: Found

    Circuit => CurrentModel % Circuits(p)
    nm = Asolver % Matrix % NumberOfRows
    CM => CurrentModel%CircuitMatrix

    DO CompInd = 1, Circuit % n_comp
      Comp => Circuit % Components(CompInd)

      Comp % Resistance = 0._dp 
      Comp % Conductance = 0._dp 

      Cvar => Comp % vvar
      vvarId = Comp % vvar % ValueId + nm
      IvarId = Comp % ivar % ValueId + nm

      CompParams => CurrentModel % Components(Comp % ComponentId) % Values
      IF (.NOT. ASSOCIATED(CompParams)) CALL Fatal ('AddComponentEquationsAndCouplings', &
               'Component parameters not found')

      IF (Comp % CoilType == 'stranded') THEN
        Comp % Resistance = ListGetCReal(CompParams, 'Resistance', Found)
        IF (Found) THEN
          Comp % UseCoilResistance = .TRUE.
        ELSE
          Comp % UseCoilResistance = .FALSE.
        END IF
      END IF

      IF ( Cvar % Owner == ParEnv % myPE .OR. .NOT. Circuit % Parallel ) THEN
        SELECT CASE (Comp % CoilType)
        CASE('stranded')
          IF (Comp % UseCoilResistance) THEN
            CALL Info('AddComponentEquationsAndCouplings', &
                'Using coil resistance for component '//i2s(CompInd), Level = 7)
            CALL AddToCmplxMatrixElement(CM, VvarId, IvarId, Comp % Resistance, 0._dp)
          ELSE
            Comp % Resistance = 0._dp
          END IF
 
          CALL AddToCmplxMatrixElement(CM, VvarId, VvarId, -1._dp, 0._dp)
        CASE('massive')
          i_multiplier = Comp % i_multiplier_re + im * Comp % i_multiplier_im
          IF (i_multiplier /= 0_dp) THEN
            CALL AddToCmplxMatrixElement(CM, VvarId, IvarId, -REAL(i_multiplier), -AIMAG(i_multiplier))
          ELSE
            CALL AddToCmplxMatrixElement(CM, VvarId, IvarId, -1._dp, 0._dp)
          END IF
        CASE('foil winding')
          ! Foil Winding voltage: 
          ! V + ...added next... = 0
          ! ----------------------
          i_multiplier = Comp % i_multiplier_re + im * Comp % i_multiplier_im
          CALL AddToCmplxMatrixElement(CM, VvarId, VvarId, 1._dp, 0._dp)
          
          IF (i_multiplier == 0_dp) i_multiplier = 1.0_dp
          
          DO j = 1, Cvar % pdofs 
            ! Foil Winding voltage: 
            !  ... - Nf/Lalpha * int_0^{Lalpha}(V_0+V_1*alpha+V_2*alpha**2+...) = 0
            !          => ... - Nf * (V_0*Lalpha^0 + V_1/2*Lalpha^1 + V_2/3*Lalpha^2 + ...) = 0
            ! where V_m is the mth dof of the polynomial
            ! --------------------------------------------------------------
            cmplx_val = -i_multiplier * REAL(Comp % nofturns) / REAL(j) * Comp % coilthickness**(j-1)
            CALL AddToCmplxMatrixElement(CM, VvarId, 2*j + VvarId, &
                REAL(cmplx_val), AIMAG(cmplx_val))

            ! Circuit eqns for the pdofs:
            ! - Nf/Lalpha * I * int_0^1(Vi'(alpha)) + ...added later... = 0
            ! ----------------------------------------------------------
            CALL AddToCmplxMatrixElement(CM, 2*j + VvarId, IvarId, &
               REAL(cmplx_val), AIMAG(cmplx_val))
          END DO
        END SELECT
      END IF

      DO q=GetNOFActive(),1,-1
        Element => GetActiveElement(q)
        CALL AddComponentElementContributions(Element, Comp, Tcoef, &
                                              sigma_33, sigmaim_33, .False.)
      END DO

      DO q=GetNOFBoundaryElements(),1,-1
        Element => GetBoundaryElement(q)
        CALL AddComponentElementContributions(Element, Comp, Tcoef, &
                                              sigma_33, sigmaim_33, .True.)
      END DO
    END DO

    IF( Circuit % Parallel ) THEN
      DO CompInd = 1, Circuit % n_comp
        Comp => Circuit % Components(CompInd)
        Comp % Resistance = ParallelReduction(Comp % Resistance)
        Comp % Conductance = ParallelReduction(Comp % Conductance)
      END DO
    END IF
      
    IF (ALLOCATED(Tcoef)) THEN
      DEALLOCATE(Tcoef,sigma_33,sigmaim_33)
    END IF
!------------------------------------------------------------------------------
   END SUBROUTINE AddComponentEquationsAndCouplings
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE AddComponentElementContributions(Element, Comp, Tcoef, & 
                                               sigma_33, sigmaim_33, boundary)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    !------------------------------------
    TYPE(Element_t), POINTER :: Element
    TYPE(Component_t), POINTER :: Comp
    COMPLEX(KIND=dp), ALLOCATABLE :: Tcoef(:,:,:)
    REAL(KIND=dp), ALLOCATABLE :: sigma_33(:), sigmaim_33(:)
    !------------------------------------
    TYPE(Solver_t), POINTER :: ASolver
    TYPE(Valuelist_t), POINTER :: CompParams
    LOGICAL :: StrandedHomogenization
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilType
    INTEGER :: nn_elem, nd_elem
    INTEGER :: astat
    LOGICAL :: Found, FoundIm, boundary

    IF (ElAssocToComp(Element, Comp)) THEN
      ASolver => CurrentModel % Asolver
      IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('AddComponentEquationsAndCouplings','ASolver not found!')

      CompParams => GetComponentParams( Element )
      IF (.NOT. ASSOCIATED(CompParams)) CALL Fatal ('AddComponentElementContributions',&
                                                    'Component parameters not found')

      StrandedHomogenization = .FALSE.
      CoilType = GetString(CompParams, 'Coil Type', Found)
      IF (.NOT. Found) CoilType = ''
      
      nn_elem = GetElementNOFNodes(Element)
      nd_elem = GetElementNOFDOFs(Element,ASolver)

      IF (.NOT. ALLOCATED(Tcoef)) THEN
        ALLOCATE(Tcoef(3,3,nn_elem), sigma_33(nn_elem), sigmaim_33(nn_elem), STAT=astat)
        IF (astat /= 0) THEN
          CALL Fatal ('AddComponentEquationsAndCouplings','Memory allocation failed')
        END IF
      ELSE IF (SIZE(Tcoef,3) /= nn_elem) THEN
        DEALLOCATE(Tcoef, sigma_33, sigmaim_33)
        ALLOCATE(Tcoef(3,3,nn_elem),sigma_33(nn_elem), sigmaim_33(nn_elem), STAT=astat)
        IF (astat /= 0) THEN
          CALL Fatal ('AddComponentEquationsAndCouplings','Memory allocation failed')
        END IF
      END IF
      
      SELECT CASE(CoilType)
      CASE ('stranded')
        StrandedHomogenization = GetLogical(CompParams, 'Homogenization Model', Found)
        IF ( StrandedHomogenization ) THEN 
          sigma_33 = GetReal(CompParams, 'sigma 33', Found)
          IF ( .NOT. Found ) sigma_33 = 0._dp
          sigmaim_33 = GetReal(CompParams, 'sigma 33 im', FoundIm)
          IF ( .NOT. FoundIm ) sigmaim_33 = 0._dp
          IF ( .NOT. Found .AND. .NOT. FoundIm ) CALL Fatal ('LocalMatrix', 'Homogenization Model sigma 33 not found!')
          IF ( .NOT. Found .AND. .NOT. FoundIm ) CALL Fatal ('AddComponentEquationsAndCouplings', &
              'Homogenization Model Sigma 33 not found!')
          Tcoef = CMPLX(0._dp, 0._dp, KIND=dp)
          Tcoef(3,3,1:nn_elem) = CMPLX(sigma_33, sigmaim_33, KIND=dp)
        ELSE
          Tcoef = GetCMPLXElectricConductivityTensor(Element, nn_elem, .TRUE., CoilType) 
        END IF
        CALL Add_stranded(Element,Tcoef,Comp,nn_elem,nd_elem,CompParams)
      CASE ('massive')
        IF (HasSupport(Element,nn_elem)) THEN
          Tcoef = GetCMPLXElectricConductivityTensor(Element, nn_elem, .TRUE., CoilType) 
          CALL Add_massive(Element,Tcoef,Comp,nn_elem,nd_elem)
        END IF
      CASE ('foil winding')
        IF (HasSupport(Element,nn_elem)) THEN
          Tcoef = GetCMPLXElectricConductivityTensor(Element, nn_elem, .TRUE., CoilType) 
          CALL Add_foil_winding(Element,Tcoef,Comp,nn_elem,nd_elem,CompParams)
        END IF
      CASE DEFAULT
        CALL Fatal ('AddComponentEquationsAndCouplings', 'Non existent Coil Type Chosen!')
      END SELECT
    END IF
!------------------------------------------------------------------------------
   END SUBROUTINE AddComponentElementContributions
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE Add_stranded(Element,Tcoef,Comp,nn,nd,CompParams)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER :: Element
    COMPLEX(KIND=dp) :: Tcoef(3,3,nn)
    TYPE(Component_t) :: Comp
    TYPE(Valuelist_t), POINTER :: CompParams
    INTEGER :: nn, nd, nm, Indexes(nd),VvarId,IvarId

    TYPE(Solver_t), POINTER :: ASolver
    INTEGER, POINTER :: PS(:)
    TYPE(Matrix_t), POINTER :: CM
    REAL(KIND=dp) :: Omega
    TYPE(Nodes_t), SAVE :: Nodes
    REAL(KIND=dp) :: Basis(nd), DetJ, x, circ_eq_coeff
    REAL(KIND=dp) :: dBasisdx(nd,3), wBase(nn), w(3)
    COMPLEX(KIND=dp) :: localC, i_multiplier, cmplx_val
    REAL(KIND=dp) :: localR !, localL
    INTEGER :: j,t
    LOGICAL :: stat

    TYPE(GaussIntegrationPoints_t) :: IP
    COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
    LOGICAL :: CSymmetry, First=.TRUE.

    REAL(KIND=dp) :: WBasis(nd,3), RotWBasis(nd,3)
    INTEGER :: dim, ncdofs,q
    
    LOGICAL :: CoilUseWvec=.FALSE., CoilUseWvec0=.FALSE.,Found,Found2
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilWVecVarname, CoilType

    LOGICAL :: PiolaVersion = .FALSE.

    TYPE(VariableHandle_t), SAVE :: Wvec_h
    TYPE(Variable_t), POINTER, SAVE :: Wpot
    
    SAVE CSymmetry, dim, First

    IF (First) THEN
      First = .FALSE.
      CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )
      dim = CoordinateSystemDimension()

      CoilUseWvec0 = GetLogical(CurrentModel % Solver % Values, 'Coil Use W Vector', Found2 ) 
      DO i=1,Model % NumberOfComponents
        CoilType = ListGetString(CurrentModel % Components(i) % Values, 'Coil Type',Found)
        IF(.NOT. Found) CYCLE
        IF(CoilType == 'stranded') THEN  ! massive, foil winding
          CoilWVecVarName = GetString(CurrentModel % Components(i) % Values,'W Vector Variable Name', Found)
          IF(Found) EXIT
        END IF
      END DO
      IF(.NOT. Found) THEN
        CoilWVecVarName = GetString(CurrentModel % Solver % Values,'W Vector Variable Name', Found)
        IF(.NOT. Found) THEN
          IF( GetLogical(CurrentModel % Solver % Values,'Use Nodal CoilCurrent',Found ) ) &
              CoilWVecVarname = 'CoilCurrent'
        END IF
         IF(.NOT. Found) THEN
          IF( GetLogical(CurrentModel % Solver % Values,'Use Elemental CoilCurrent',Found ) ) &
              CoilWVecVarname = 'CoilCurrent e'
        END IF
        IF(Found) CALL Info('Add_stranded','Setting coil current to: '//TRIM(CoilWVecVarname),Level=6)
        ! If we did not find w vector named in any component it is fair to assume that it is globally used!
        IF(.NOT. Found2) CoilUseWvec0 = Found 
      END IF
      IF(.NOT. Found) CoilWVecVarname = 'W Vector E'
      CALL ListInitElementVariable(Wvec_h, CoilWVecVarname)

      CALL GetWPotentialVar(Wpot)
    END IF

    ASolver => CurrentModel % Asolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('Add_stranded','ASolver not found!')
    PiolaVersion = GetLogical( ASolver % Values, 'Use Piola Transform', Found )

    PS => Asolver % Variable % Perm

    CM => CurrentModel % CircuitMatrix
    nm = CurrentModel % Asolver % Matrix % NumberOfRows
    Omega = GetAngularFrequency()
    
    CALL GetElementNodes(Nodes)
    nd = GetElementDOFs(Indexes,Element,ASolver)
    
    ncdofs=nd
    IF (dim == 3) THEN
      ncdofs=nd-nn

      CoilUseWvec = GetLogical(CompParams, 'Coil Use W Vector', Found)
      IF(.NOT. Found) CoilUseWvec = CoilUseWvec0 

      IF (.NOT. CoilUseWvec) THEN
        CALL GetLocalSolution( Wbase,UElement=Element,UVariable=Wpot, Found=Found)
        !CALL GetWPotential(WBase)
      END IF
    END IF

    VvarId = Comp % vvar % ValueId + nm
    IvarId = Comp % ivar % ValueId + nm

    i_multiplier = Comp % i_multiplier_re + im * Comp % i_multiplier_im


    ! Numerical integration:
    ! ----------------------
    IF(PiolaVersion) THEN
      IP = GaussPoints(Element, EdgeBasis=.TRUE., EdgeBasisDegree=1, PReferenceElement=PiolaVersion)
    ELSE
      IP = GaussPoints(Element)
    END IF

    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                IP % W(t), detJ, Basis,dBasisdx )

 
      circ_eq_coeff = 1._dp
      SELECT CASE(dim)
      CASE(2)
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                  IP % W(t), detJ, Basis, dBasisdx )
        w = [0._dp, 0._dp, 1._dp]
        IF( CSymmetry ) THEN
          x = SUM( Basis(1:nn) * Nodes % x(1:nn) )
          detJ = detJ * x
        END IF
        circ_eq_coeff = GetCircuitModelDepth()
      CASE(3)
        IF (PiolaVersion) THEN
          stat = EdgeElementInfo( Element, Nodes, IP % U(t), IP % V(t), IP % W(t), &
               DetF = DetJ, Basis = Basis, dBasisdx=dBasisdx, EdgeBasis = WBasis, RotBasis = RotWBasis, &
               BasisDegree=1, ApplyPiolaTransform = .TRUE.)
        ELSE
          stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                    IP % W(t), detJ, Basis, dBasisdx )
          CALL GetEdgeBasis(Element,WBasis,RotWBasis,Basis,dBasisdx)
        END IF
        IF (CoilUseWvec) THEN
          w = ListGetElementVectorSolution( Wvec_h, Basis, Element, Found = Found, dofs = dim )
          IF(.NOT. Found ) THEN
            CALL Fatal('Add_stranded','Could not find coil current density!')            
          END IF
        ELSE
          w = -MATMUL(WBase(1:nn), dBasisdx(1:nn,:))
        END IF
      END SELECT

      localC = SUM(Tcoef(3,3,1:nn) * Basis(1:nn))
      
      IF (.NOT. Comp % UseCoilResistance) THEN
        ! I * R, where 
        ! R = (1/sigma * js,js):
        ! ----------------------
        localR = Comp % N_j **2 * IP % s(t)*detJ*SUM(w*w)/localC*circ_eq_coeff / Comp % VoltageFactor

        Comp % Resistance = Comp % Resistance + localR
        
        CALL AddToCmplxMatrixElement(CM, VvarId, IvarId, &
              REAL(Comp % N_j**2 * IP % s(t)*detJ*SUM(w*w)/localC*circ_eq_coeff / Comp % VoltageFactor), &
             AIMAG(Comp % N_j**2 * IP % s(t)*detJ*SUM(w*w)/localC*circ_eq_coeff / Comp % VoltageFactor))
      END IF
      
      DO j=1,ncdofs
        q=j
        IF (dim == 3) q=q+nn
        IF (Comp % N_j/=0._dp) THEN
          ! ( im * Omega a,w )
          IF (dim == 2) cmplx_val = im * Omega * Comp % N_j &
                  * IP % s(t)*detJ*Basis(j)*circ_eq_coeff*w(3)
          IF (dim == 3) cmplx_val = im * Omega * Comp % N_j &
                  * IP % s(t)*detJ*SUM(WBasis(j,:)*w)

          cmplx_val = cmplx_val / Comp % VoltageFactor
!          localL = ABS(cmplx_val)
!          Comp % Inductance = Comp % Inductance + localL

          CALL AddToCmplxMatrixElement(CM, VvarId, ReIndex(PS(Indexes(q))), &
                 REAL(cmplx_val), AIMAG(cmplx_val))
          
          IF (dim == 2) cmplx_val = -Comp % N_j*IP % s(t)*detJ*Basis(j)*w(3)
          IF (dim == 3) cmplx_val = -Comp % N_j*IP % s(t)*detJ*SUM(WBasis(j,:)*w)
          IF (i_multiplier /= 0._dp) cmplx_val = i_multiplier*cmplx_val

          cmplx_val = cmplx_val * Comp % SymmetryCoeff
          
          CALL AddToCmplxMatrixElement(CM,ReIndex(PS(Indexes(q))), IvarId, &
             REAL(cmplx_val), AIMAG(cmplx_val))

        END IF
      END DO
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE Add_stranded
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE Add_massive(Element,Tcoef,Comp,nn,nd)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t) :: Element
    COMPLEX(KIND=dp) :: Tcoef(3,3,nn)
    TYPE(Component_t) :: Comp

    TYPE(Solver_t), POINTER :: ASolver
    TYPE(ValueList_t), POINTER :: BC
    INTEGER, POINTER :: PS(:)
    TYPE(Matrix_t), POINTER :: CM
    REAL(KIND=dp) :: Omega, grads_coeff, circ_eq_coeff
    REAL(KIND=dp) :: Basis(nd), DetJ, x
    REAL(KIND=dp) :: dBasisdx(nd,3)
    REAL(KIND=dp) :: LondonLambda(nn)
    REAL(KIND=dp) :: LondonLambda_ip, val
    REAL(KIND=dp) :: SkinCond(nn), SkinMu(nn)
    REAL(KIND=dp) :: cond, mu, muVacuum, delta
    COMPLEX(KIND=dp) :: imu, invZs
    COMPLEX(KIND=dp) :: localC, cmplx_val=0._dp
    REAL(KIND=dp) :: localConductance !, localL
    INTEGER :: nn, nd, j, t, nm, Indexes(nd), &
               VvarId, dim
    LOGICAL :: stat
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
    LOGICAL :: CSymmetry, First=.TRUE.
    LOGICAL :: LondonEquations, SkinBc=.False., ElectroDynamics
    TYPE(ValueList_t), POINTER :: Material

    REAL(KIND=dp) :: wBase(nn), gradv(3), WBasis(nd,3), RotWBasis(nd,3)
    INTEGER :: ncdofs,q
    REAL(KIND=dp) :: ModelDepth
    COMPLEX(KIND=dp) :: Permittivity(nn), localP
    TYPE(Variable_t), POINTER, SAVE :: Wpot

    
    SAVE CSymmetry, dim, First

    IF (First) THEN
      First = .FALSE.
      CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )
      dim = CoordinateSystemDimension()

      CALL GetWPotentialVar(Wpot)
    END IF

    ASolver => CurrentModel % Asolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('Add_massive','ASolver not found!')
    PS => Asolver % Variable % Perm

    CM => CurrentModel % CircuitMatrix
    nm = CurrentModel % Asolver % Matrix % NumberOfRows
    Omega = GetAngularFrequency()

    SkinCond = 0._dp
    SkinMu = 0._dp

    CALL GetElementNodes(Nodes)
    nd = GetElementDOFs(Indexes,Element,ASolver)

    ncdofs=nd
    IF (dim == 3) THEN
      !CALL GetWPotential(WBase)     
      CALL GetLocalSolution( Wbase,UElement=Element,UVariable=Wpot, Found=Found)
      ncdofs=nd-nn
    END IF

    Material => GetMaterial( Element )

    ElectroDynamics = GetLogical( Asolver % values, 'Electrodynamics Model', Found)
    IF(ElectroDynamics) THEN
      CALL GetPermittivity(Material, Permittivity, nn)
    END IF

    LondonLambda(:) = GetReal( Material, 'London Lambda', LondonEquations, Element)

    BC => GetBC( Element )
    skinBc = .FALSE.
    IF ( ASSOCIATED(BC) ) THEN
      SkinCond = GetConstReal( BC, 'Layer Electric Conductivity', SkinBc)
      IF ( SkinBc ) THEN
        muVacuum = 4 * PI * 1d-7
        imu = CMPLX(0.0_dp, 1.0_dp, KIND=dp) 
        SkinMu = GetConstReal( BC, 'Layer Relative Permeability', Found)
      END IF
    END IF
  
    vvarId = Comp % vvar % ValueId + nm

    ! Numerical integration:
    ! ----------------------
    IP = GaussPoints(Element)
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                IP % W(t), detJ, Basis,dBasisdx )

      grads_coeff = -1._dp
      circ_eq_coeff = 1._dp
      SELECT CASE(dim)
      CASE(2)
        IF( CSymmetry ) THEN
          x = SUM( Basis(1:nn) * Nodes % x(1:nn) )
          detJ = detJ * x
          grads_coeff = grads_coeff/x
        END IF
        ModelDepth = GetCircuitModelDepth()
        circ_eq_coeff = ModelDepth
        grads_coeff = grads_coeff/ModelDepth
      CASE(3)
        CALL GetEdgeBasis(Element,WBasis,RotWBasis,Basis,dBasisdx)
        gradv = MATMUL( WBase(1:nn), dBasisdx(1:nn,:))
      END SELECT


      localC = SUM(Tcoef(1,1,1:nn) * Basis(1:nn))

      localP = 0._dp
      IF(ElectroDynamics) THEN
        localP = SUM( Basis(1:nn) * Permittivity(1:nn))
      END IF

      ! computing the source term Vi(sigma grad v0, grad si):
      ! ------------------------------------------------
      IF(dim==2) cmplx_val = IP % s(t)*detJ*grads_coeff**2*circ_eq_coeff * Comp % VoltageFactor

      IF(dim==3) THEN
        IF (SkinBc) THEN
        ! if SkinBC is activated:
        !   Boundary Condition: Layer Electric Conductivity
        !   Boundary Condition: Layer Relative Permeability
        !
        ! The term Vi/Z ( grad v0, grad_v0 )
        !
          cond = SUM(Basis(1:nn) * SkinCond(1:nn))
          mu  = muVacuum * SUM(Basis(1:nn) * SkinMu(1:nn))
          delta = SQRT( 2.0_dp/(cond*omega*mu))      
          invZs = (cond*delta)/(1.0_dp+imu)
          cmplx_val = IP % s(t)*detJ*invZs*SUM(gradv*gradv) * Comp % VoltageFactor
        ELSE
          cmplx_val = IP % s(t)*detJ*SUM(gradv*gradv) * Comp % VoltageFactor
        END IF
      END IF

      IF(SkinBC) THEN
        CALL AddToCmplxMatrixElement(CM, vvarId, vvarId, &
             REAL(cmplx_val), AIMAG(cmplx_val))
      ELSE
        CALL AddMatrixEntry(CM, vvarId, vvarId, &
             cmplx_val*localC, im*Omega*cmplx_val*localP )
      END IF

      localConductance = ABS(cmplx_val*localC)
      Comp % Conductance = Comp % Conductance + localConductance

      IF ( LondonEquations ) THEN
        LondonLambda_ip = SUM( Basis(1:nn) * LondonLambda(1:nn) )

        IF(dim==2) val = IP % s(t)*detJ/LondonLambda_ip*grads_coeff**2*circ_eq_coeff
        val = val * Comp % VoltageFactor
        ! Phi (beta grad phi_0, grad phi')
        ! Here Phi takes the place of Vi
        ! -----------------------------------
        CALL AddToCmplxMatrixElement(CM, vvarId, vvarId, val, 0._dp)
      END IF


      DO j=1,ncdofs
        q=j
        IF (dim == 3) q=q+nn
 
        IF ( LondonEquations ) THEN
          ! Phi * ( beta * grad phi, a')
          ! where phi is the node flux scalar potential
          ! -------------------------------------------
          IF(dim==2) val = IP % s(t)*detJ/LondonLambda_ip*basis(j)*grads_coeff*circ_eq_coeff
          CALL AddToCmplxMatrixElement(CM, vvarId, ReIndex(PS(Indexes(q))), val, 0._dp)

          IF(dim==2) val = IP % s(t)*detJ/LondonLambda_ip*basis(j)*grads_coeff
          val = val * Comp % VoltageFactor
          CALL AddToCmplxMatrixElement(CM, ReIndex(PS(indexes(q))), vvarId, val, 0._dp)
        END IF

        ! computing the mass term (sigma * im * Omega * a, grad si):
        ! ---------------------------------------------------------
        IF(dim==2) cmplx_val = IP % s(t)*detJ*Basis(j)*grads_coeff*circ_eq_coeff

        IF(dim==3) THEN
          IF (SkinBc) THEN
          ! if SkinBC is activated:
          !   Boundary Condition: Layer Electric Conductivity
          !   Boundary Condition: Layer Relative Permeability
          ! Then activate
          !  (1/Z*im*Omega* a , grad v')
          !
            cmplx_val = IP % s(t)*detJ*invZs*im*Omega*SUM(Wbasis(j,:)*gradv)
          ELSE
            cmplx_val = IP % s(t)*detJ*SUM(Wbasis(j,:)*gradv)
          END IF
        END IF

        IF(SkinBc ) THEN
          CALL AddToCmplxMatrixElement(CM, vvarId, ReIndex(PS(Indexes(q))), &
                    REAL(cmplx_val), AIMAG(cmplx_val))
        ELSE
          cmplx_val = im*Omega*cmplx_val*localC - omega**2*cmplx_val*localP
          CALL AddMatrixEntry(CM, vvarId, ReIndex(PS(Indexes(q))), &
                     (0._dp, 0._dp), cmplx_val )
        END IF


!        localL = ABS(1._dp/cmplx_val)
!        Comp % Inductance = Comp % Inductance + localL
        
        IF(dim==2) cmplx_val = IP % s(t)*detJ*basis(j)*grads_coeff * Comp % VoltageFactor

        IF(dim==3) THEN
          IF (SkinBc) THEN
            ! if SkinBC is activated:
            !   Layer Electric Conductivity = Real 58e6
            !   Layer Relative Permeability = Real 1
            !
            !  + 1/Z (grad v , a') 
            !
            cmplx_val = IP % s(t)*detJ*invZs*SUM(gradv*Wbasis(j,:)) * Comp % VoltageFactor
          ELSE
            cmplx_val = IP % s(t)*detJ*SUM(gradv*Wbasis(j,:)) * Comp % VoltageFactor 
          END IF
        END IF

        IF(SkinBc) THEN
          CALL AddToCmplxMatrixElement(CM, ReIndex(PS(indexes(q))), vvarId, &
                  REAL(cmplx_val), AIMAG(cmplx_val))
        ELSE
          CALL AddMatrixEntry(CM, ReIndex(PS(indexes(q))), vvarId, &
               cmplx_val*LocalC, im*Omega*cmplx_val*LocalP )
        END IF
      END DO
    END DO

!------------------------------------------------------------------------------
   END SUBROUTINE Add_massive
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE Add_foil_winding(Element,Tcoef,Comp,nn,nd,CompParams)
!------------------------------------------------------------------------------
    USE MGDynMaterialUtils
    IMPLICIT NONE
    INTEGER :: nn, nd
    TYPE(Element_t), POINTER :: Element
    COMPLEX(KIND=dp) :: Tcoef(3,3,nn), C(3,3), val
    TYPE(Component_t) :: Comp
    TYPE(Valuelist_t), POINTER :: CompParams

    TYPE(Solver_t), POINTER :: ASolver
    INTEGER, POINTER :: PS(:)
    TYPE(Matrix_t), POINTER :: CM
    REAL(KIND=dp) :: Basis(nd), DetJ, Omega, localAlpha, localV, localVtest, &
                     x, circ_eq_coeff, grads_coeff
    REAL(KIND=dp) :: dBasisdx(nd,3),alpha(nn)
    INTEGER :: nm,p,j,t,Indexes(nd),vvarId,vpolord_tot, &
               vpolord, vpolordtest, dofId, dofIdtest, &
               dim
    LOGICAL :: stat
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
    LOGICAL :: CSymmetry, First=.TRUE., InitHandle=.TRUE., &
               CoilUseWvec=.FALSE., CoilUseWvec0=.FALSE.,Found,Found2
    LOGICAL :: InitJHandle=.TRUE., FoilUseJvec=.FALSE.
    REAL(KIND=dp) :: localR
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilWVecVarname, CoilType
    CHARACTER(LEN=MAX_NAME_LEN) :: FoilJVecVarname
    TYPE(VariableHandle_t), SAVE :: Wvec_h
    TYPE(VariableHandle_t), SAVE :: Jvec_h

    REAL(KIND=dp) :: wBase(nn), gradv(3), WBasis(nd,3), RotWBasis(nd,3), &
                     RotMLoc(3,3), RotM(3,3,nn)
    REAL(KIND=dp) :: Jvec(3)
    INTEGER :: i,ncdofs,q
    TYPE(Variable_t), POINTER, SAVE :: Wpot

    
    SAVE CSymmetry, dim, First, InitHandle, InitJHandle

    IF( First ) THEN
      First = .FALSE.
      CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )
      dim = CoordinateSystemDimension()
      
      CoilUseWvec0 = GetLogical(CurrentModel % Solver % Values, 'Coil Use W Vector', Found2 ) 
      DO i=1,Model % NumberOfComponents
        CoilType = ListGetString(CurrentModel % Components(i) % Values, 'Coil Type',Found)
        IF(.NOT. Found) CYCLE
        IF(CoilType == 'foil winding') THEN  ! stranded, massive
          CoilWVecVarName = GetString(CurrentModel % Components(i) % Values,'W Vector Variable Name', Found)
          IF(Found) EXIT
        END IF
      END DO
      IF(.NOT. Found) THEN
        CoilWVecVarName = GetString(CurrentModel % Solver % Values,'W Vector Variable Name', Found)
        IF(.NOT. Found) THEN
          IF( GetLogical(CurrentModel % Solver % Values,'Use Nodal CoilCurrent',Found ) ) &
              CoilWVecVarname = 'CoilCurrent'
        END IF
        IF(.NOT. Found) THEN
          IF( GetLogical(CurrentModel % Solver % Values,'Use Elemental CoilCurrent',Found ) ) &
              CoilWVecVarname = 'CoilCurrent e'
        END IF
        IF(Found) CALL Info('Add_foil_winding','Setting coil current to: '//TRIM(CoilWVecVarname),Level=6)
        ! If we did not find w vector named in any component it is fair to assume that it is globally used!
        IF(.NOT. Found2) CoilUseWvec0 = Found 
      END IF
      IF(.NOT. Found) CoilWVecVarname = 'W Vector E'
      CALL ListInitElementVariable(Wvec_h, CoilWVecVarname)

      CALL GetWPotentialVar(Wpot)
    END IF

    ASolver => CurrentModel % Asolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('Add_foil_winding','ASolver not found!')
    PS => Asolver % Variable % Perm

    CM => CurrentModel % CircuitMatrix
    nm = CurrentModel % Asolver % Matrix % NumberOfRows
    Omega = GetAngularFrequency()

    CALL GetElementNodes(Nodes)
    nd = GetElementDOFs(Indexes,Element,ASolver)
    CALL GetLocalSolution(alpha,'Alpha')
    
    ncdofs=nd
    IF (dim == 3) THEN

      ! If we do not have a local flag then use the one from the solver section
      CoilUseWvec = GetLogical(CompParams, 'Coil Use W Vector', Found)
      IF (.NOT. Found) CoilUseWvec = CoilUseWvec0

      IF (.NOT. CoilUseWvec) THEN
        !CALL GetLocalSolution(Wbase, 'w')
        !CALL GetWPotential(WBase)
        CALL GetLocalSolution( Wbase,UElement=Element,UVariable=Wpot, Found=Found)
      END IF

      FoilUseJvec = GetLogical(CompParams, 'Foil Winding Use J Vector', Found)
      IF (.NOT. Found) FoilUseJvec = .FALSE.

      IF (FoilUseJvec) THEN
        IF( InitJHandle ) THEN
          FoilJVecVarname = GetString(CompParams, 'Foil J Vector Variable Name', Found)
          IF ( .NOT. Found) FoilJVecVarname = 'J Vector E'
          CALL ListInitElementVariable(Jvec_h, FoilJVecVarname)
          IF ( .NOT. ASSOCIATED(Jvec_h % Variable)) THEN
            CALL Fatal('Add_foil_winding','You are trying to use Foil J Vector for describing the &
                                    component source field but I cannot the variable')
          END IF
          InitJHandle = .FALSE.
        END IF
      END IF

      CALL GetElementRotM(Element, RotM, nn)

      ncdofs=nd-nn
    END IF

    vvarId = Comp % vvar % ValueId
    vpolord_tot = Comp % vvar % pdofs - 1

    ! Numerical integration:
    ! ----------------------
    IP = GaussPoints(Element)
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                IP % W(t), detJ, Basis,dBasisdx )

      grads_coeff = -1._dp
      circ_eq_coeff = 1._dp
      SELECT CASE(dim)
      CASE(2)
        IF( CSymmetry ) THEN
          x = SUM( Basis(1:nn) * Nodes % x(1:nn) )
          detJ = detJ * x
          grads_coeff = grads_coeff/x
        END IF
        circ_eq_coeff = GetCircuitModelDepth()
        grads_coeff = grads_coeff/circ_eq_coeff
        C(1,1) = SUM( Tcoef(3,3,1:nn) * Basis(1:nn) )
        ! I * R, where 
        ! R = (1/sigma * js,js):
        ! ----------------------
        localR = Comp % N_j **2 * IP % s(t)*detJ/C(1,1)*circ_eq_coeff / Comp % VoltageFactor
      CASE(3)
        CALL GetEdgeBasis(Element,WBasis,RotWBasis,Basis,dBasisdx)

        IF (CoilUseWvec) THEN
          gradv = ListGetElementVectorSolution( Wvec_h, Basis, Element, dofs = dim )
        ELSE
          gradv = MATMUL( WBase(1:nn), dBasisdx(1:nn,:))
        END IF

        ! Compute the conductivity tensor
        ! -------------------------------
        DO i=1,3
          DO j=1,3
            C(i,j) = SUM( Tcoef(i,j,1:nn) * Basis(1:nn) )
            RotMLoc(i,j) = SUM( RotM(i,j,1:nn) * Basis(1:nn) )
          END DO
        END DO

        ! I * R, where 
        ! R = (1/sigma * js,js):
        ! ----------------------
        localR = Comp % N_j **2 * IP % s(t)*detJ/C(3,3) / Comp % VoltageFactor

        C = MATMUL(MATMUL(RotMLoc, C),TRANSPOSE(RotMLoc))

        IF (FoilUseJvec) THEN
          Jvec = ListGetElementVectorSolution( Jvec_h, Basis, Element, dofs = dim )
        ELSE
          Jvec = MATMUL(C,gradv)
        END IF

        ! Transform the conductivity tensor:
        ! ----------------------------------
      END SELECT
      
      localAlpha = SUM(alpha(1:nn) * Basis(1:nn))
      
      ! alpha is normalized to be in [0,1] thus, 
      ! it needs to be multiplied by the thickness of the coil 
      ! to get the real alpha:
      ! ------------------------------------------------------
      localAlpha = localAlpha * Comp % coilthickness

      Comp % Resistance = Comp % Resistance + localR

      DO vpolordtest=0,vpolord_tot ! V'(alpha)
        localVtest = localAlpha**vpolordtest
        dofIdtest = 2*(vpolordtest + 1) + vvarId
        DO vpolord = 0, vpolord_tot ! V(alpha)

          localV = localAlpha**vpolord
          dofId = 2*(vpolord + 1) + vvarId

          ! Computing the stiff term (sigma V(alpha) grad v0, V'(alpha) grad si):
          ! ---------------------------------------------------------------------
          IF (dim == 2) val = IP % s(t)*detJ*localV*localVtest*C(1,1)*grads_coeff**2*circ_eq_coeff
          IF (dim == 3) val = IP % s(t)*detJ*localV*localVtest*SUM(Jvec*gradv)
          val = val * Comp % VoltageFactor

          CALL AddToCmplxMatrixElement(CM, dofIdtest+nm, dofId+nm, REAL(val), AIMAG(val))
        END DO

        DO j=1,ncdofs
          q=j
          IF (dim == 3) q=q+nn
          ! computing the mass term (sigma * im * Omega * a, V'(alpha) grad si):
          ! ---------------------------------------------------------
          IF (dim == 2) val = im * Omega * IP % s(t)*detJ*localVtest*C(1,1)*basis(j)*grads_coeff*circ_eq_coeff
          IF (dim == 3) val = im * Omega * IP % s(t)*detJ*localVtest*SUM(MATMUL(C,Wbasis(j,:))*gradv)
          CALL AddToCmplxMatrixElement(CM, dofIdtest+nm, ReIndex(PS(Indexes(q))), REAL(val), AIMAG(val) )
        END DO

      END DO

      DO vpolord = 0, vpolord_tot ! V(alpha)
        localV = localAlpha**vpolord
        dofId = 2*(vpolord + 1) + vvarId

        DO j=1,ncdofs
            q=j
            IF (dim == 3) q=q+nn
            IF (dim == 2) val = IP % s(t)*detJ*localV*C(1,1)*basis(j)*grads_coeff
            IF (dim == 3) val = IP % s(t)*detJ*localV*SUM(Jvec*Wbasis(j,:))
            val = val * Comp % VoltageFactor
            CALL AddToCmplxMatrixElement(CM, ReIndex(PS(indexes(q))), dofId+nm, REAL(val), AIMAG(val))
        END DO
      END DO

    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE Add_foil_winding
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE GetConductivity(Element, Tcoef, nn)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER :: Element
    TYPE(Valuelist_t), POINTER :: Material
    COMPLEX(KIND=dp) :: Tcoef(3,3,nn)
    REAL(KIND=dp), POINTER, SAVE :: Cwrk(:,:,:), Cwrk_im(:,:,:) 
    INTEGER :: nn, i, j
    LOGICAL, SAVE :: visited = .FALSE., Found

    IF (.NOT. visited) THEN
      NULLIFY( Cwrk, Cwrk_im )
    END IF

    Tcoef = cmplx(0.0d0,0.0d0)
    Material => GetMaterial( Element )
    IF (.NOT. ASSOCIATED(Material)) CALL Fatal('Circuits_apply','Material not found.')

    CALL ListGetRealArray( Material, &
           'Electric Conductivity', Cwrk, nn, Element % NodeIndexes, Found )

    IF (.NOT. Found) CALL Fatal('Circuits_apply', 'Electric Conductivity not found.')
    
    IF (Found) THEN
       IF ( SIZE(Cwrk,1) == 1 ) THEN
          DO i=1,3
             Tcoef( i,i,1:nn ) = Cwrk( 1,1,1:nn )
          END DO
       ELSE IF ( SIZE(Cwrk,2) == 1 ) THEN
          DO i=1,MIN(3,SIZE(Cwrk,1))
             Tcoef(i,i,1:nn) = Cwrk(i,1,1:nn)
          END DO
       ELSE
          DO i=1,MIN(3,SIZE(Cwrk,1))
             DO j=1,MIN(3,SIZE(Cwrk,2))
                Tcoef( i,j,1:nn ) = Cwrk(i,j,1:nn)
             END DO
          END DO
       END IF
    END IF

    CALL ListGetRealArray( Material, &
           'Electric Conductivity im', Cwrk_im, nn, Element % NodeIndexes, Found )

    IF (Found) THEN
       IF ( SIZE(Cwrk_im,1) == 1 ) THEN
          DO i=1,3
             Tcoef( i,i,1:nn ) = CMPLX( REAL(Tcoef( i,i,1:nn )), Cwrk_im( 1,1,1:nn ), KIND=dp)
          END DO
       ELSE IF ( SIZE(Cwrk_im,2) == 1 ) THEN
          DO i=1,MIN(3,SIZE(Cwrk_im,1))
             Tcoef(i,i,1:nn) = CMPLX( REAL(Tcoef( i,i,1:nn )), Cwrk_im( i,1,1:nn ), KIND=dp)
          END DO
       ELSE
          DO i=1,MIN(3,SIZE(Cwrk_im,1))
             DO j=1,MIN(3,SIZE(Cwrk_im,2))
                Tcoef( i,j,1:nn ) = CMPLX( REAL(Tcoef( i,j,1:nn )), Cwrk_im( i,j,1:nn ), KIND=dp)
             END DO
          END DO
       END IF
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE GetConductivity
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE AddMatrixEntry( A, row, col, val, tval )
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: A
    INTEGER :: row, col
    COMPLEX(KIND=dp) :: val, tval
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp) :: cval
    REAL(KIND=dp), POINTER :: svalues(:)
!------------------------------------------------------------------------------

    IF(EigenSystem) THEN
      CALL AddToCmplxMatrixElement(A, row, col, REAL(val), AIMAG(val))

      IF(.NOT.ASSOCIATED(A % MassValues)) THEN
        ALLOCATE(A % MassValues(SIZE(A % Values)))
        A % MassValues = 0._dp
      END IF

      svalues => A % Values
      A % Values => A % MassValues

      CALL AddToCmplxMatrixElement(A, row, col, REAL(tval), AIMAG(tval))

      A % Values => svalues
    ELSE
      cval = val+tval
      CALL AddToCmplxMatrixElement(A, row, col, REAL(cval), AIMAG(cval))
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE AddMatrixEntry
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE CircuitsAndDynamicsHarmonic
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
SUBROUTINE CircuitsOutput(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
   USE DefUtils
   USE CircuitUtils
   USE CircuitsMod
   IMPLICIT NONE
!------------------------------------------------------------------------------   
   TYPE(Model_t) :: Model
   TYPE(Solver_t) :: Solver
   REAL(KIND=dp) :: dt
   LOGICAL :: Transient
!------------------------------------------------------------------------------

   TYPE(Variable_t), POINTER :: LagrangeVar
   REAL(KIND=dp), ALLOCATABLE  :: crt(:), crtt(:)
   INTEGER :: nm
   
   TYPE(Solver_t), POINTER :: ASolver
   TYPE(Component_t), POINTER :: Comp
   TYPE(ValueList_t), POINTER :: CompParams

   CHARACTER(LEN=MAX_NAME_LEN) :: dofnumber, CompName 
   INTEGER :: i,p,jj,j
   TYPE(CircuitVariable_t), POINTER :: CVar

   TYPE(Matrix_t), POINTER :: CM    
   INTEGER, POINTER :: n_Circuits => Null(), circuit_tot_n => Null()
   TYPE(Circuit_t), POINTER :: Circuits(:)

   COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
   COMPLEX(KIND=dp) :: Current 
   REAL(KIND=dp) :: CompRealPower, p_dc_component

   LOGICAL :: Found
!------------------------------------------------------------------------------
! EEC variables
!------------------------------------------------------------------------------
  LOGICAL, SAVE :: EEC, First =.TRUE.
  LOGICAL :: EEC_lim
  REAL(KIND=dp), SAVE :: EEC_freq, EEC_time_0
  INTEGER, SAVE :: EEC_max, EEC_cnt = 0
  REAL :: TTime
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Variable_t), POINTER :: AzVar
  REAL (KIND=dp), POINTER :: Az(:)
  REAL (KIND=dp), ALLOCATABLE, SAVE :: Az0(:)
  REAL (KIND=dp), POINTER :: Acorr(:)
  CHARACTER(*), PARAMETER :: Caller = 'CircuitsOutput'
  CHARACTER(LEN=MAX_NAME_LEN), SAVE :: CktPrefix, sname
  LOGICAL :: Parallel
!------------------------------------------------------------------------------  
      
   CALL DefaultStart()

   Circuit_tot_n => Model%Circuit_tot_n
   n_Circuits => Model%n_Circuits
   CM => Model%CircuitMatrix
   Circuits => Model%Circuits

   ! Look for the solver we attach the circuit equations to:
   ! -------------------------------------------------------
   ASolver => CurrentModel % Asolver
   IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal(Caller,'ASolver not found!')
      
   nm =  Asolver % Matrix % NumberOfRows


  IF (First) THEN
    SolverParams => GetSolverParams(Solver)
    ! Reading parameter for supply frequency
    EEC_freq = GetConstReal( SolverParams, 'EEC Frequency', EEC)
    IF (EEC) THEN
      CALL Info(Caller, "Using EEC steady state forcing.", Level=4)
      WRITE( Message,'(A,4G11.4,A)') 'EEC signal frequency: ', EEC_freq, ' Hz'
      CALL Info(Caller, Message, Level=4)
                
      EEC_max = GetInteger( SolverParams, 'EEC Steps', EEC_lim)
      IF (.NOT. EEC_lim) EEC_max = 5 !Typically 5 correections is enough
      WRITE( Message,'(A,I5,A)') 'Applying ', EEC_max, ' halfperiod corrections'
      CALL Info(Caller, Message, Level=4)
      
      EEC_time_0 = 0.0
      
      ! Reserve memory for storing current MVP solution
      ALLOCATE(Az0(nm))
      
      ! Store MVP solution at t=0
      AzVar => Asolver % Variable 
      IF(ASSOCIATED( AzVar)) THEN
        Az => AzVar % Values
        Az0 = Az
      END IF
      
    END IF
    
    CktPrefix = ListGetString(SolverParams,'Scalars Prefix',Found )
    IF(.NOT. Found) CktPrefix = 'res:'

    First = .FALSE.
  END IF


  IF (EEC .AND. (EEC_cnt < EEC_max)) THEN
    
    TTime = GetTime()
    IF(TTime >= (EEC_time_0 + 0.5/EEC_freq)) THEN
      EEC_cnt = EEC_cnt + 1
      WRITE( Message,'(A,4G11.4)') 'Performing EEC #', EEC_cnt
      CALL Info(Caller, Message, Level=4)
      
      EEC_time_0 = EEC_time_0 + 0.5/EEC_freq

      AzVar => Asolver % Variable 
      
      IF(ASSOCIATED( AzVar)) THEN
        Az => AzVar % Values
        
        !calculate correction
        ALLOCATE(Acorr(nm))
        Acorr = -0.5*(Az0+Az)
        Az = Az+Acorr
        DEALLOCATE(Acorr)
        
        !Store corrected half-period solution
        Az0 = Az
      END IF
    END IF
  END IF


   ! Circuit variable values from previous timestep:
   ! -----------------------------------------------
  Parallel = ( ParEnv % PEs > 1 )
  IF( Parallel ) THEN
    IF( Solver % Mesh % SingleMesh ) Parallel = ListGetLogical( Model % Simulation,'Enforce Parallel',Found )
  END IF
    
  ALLOCATE(crt(circuit_tot_n), crtt(circuit_tot_n))
   crt = 0._dp
   crtt = 0._dp

   sname = LagrangeMultiplierName( ASolver )
   LagrangeVar => VariableGet( ASolver % Mesh % Variables,sname)
   IF(ASSOCIATED(LagrangeVar)) THEN
     CALL Info(Caller,'Initializing Lagrange multipliers of size: '&
         //I2S(SIZE(LagrangeVar % Values)),Level=8)
     IF( Parallel ) THEN
       DO i=1,circuit_tot_n 
         IF (ASSOCIATED(Model%CircuitMatrix)) THEN  
           IF( CM % RowOwner(nm+i)==Parenv%myPE) crtt(i) = LagrangeVar%Values(i)
         END IF
       END DO
       CALL MPI_ALLREDUCE(crtt,crt,circuit_tot_n, MPI_DOUBLE_PRECISION, &
                  MPI_SUM, ASolver % Matrix % Comm, j)
     ELSE
       crt(1:circuit_tot_n) = LagrangeVar % Values
     END IF
   END IF

   !IF( ListGetLogical( Solver % Values,'Store Cyclic System',Found ) ) THEN 
   !  Solver % Variable => LagrangeVar 
   !END IF
   
   ! Export circuit & dynamic variables for "SaveScalars":
   ! -----------------------------------------------------

   CALL ListAddConstReal(GetSimulation(),TRIM(CktPrefix)//' time', GetTime())

   CALL Info(Caller, 'Writing Circuit Results', Level=5) 
   DO p=1,n_Circuits
     CALL Info(Caller, 'Writing Circuit Variables for &
       Circuit '//i2s(p), Level=8) 
     CALL Info(Caller, 'There are '//i2s(Circuits(p)%n)//&
       ' Circuit Variables', Level=8)
     DO i=1,Circuits(p) % n
       Cvar => Circuits(p) % CircuitVariables(i)
       
       IF (Circuits(p) % Harmonic) THEN 
         CALL SimListAddAndOutputConstReal(&
           TRIM(Circuits(p) % names(i))//' re', crt(Cvar % ValueId), Level=10)
         CALL SimListAddAndOutputConstReal(&
           TRIM(Circuits(p) % names(i))//' im', crt(Cvar % ImValueId), Level=10)

         IF (Cvar % pdofs /= 0 ) THEN
           DO jj = 1, Cvar % pdofs
             CALL SimListAddAndOutputConstReal(&
               TRIM(Circuits(p) % names(i))&
               //'re dof '//I2S(jj), crt(Cvar % ValueId + ReIndex(jj)), Level=10)
             CALL SimListAddAndOutputConstReal(&
               TRIM(Circuits(p) % names(i))&
               //'im dof '//I2S(jj), crt(Cvar % ValueId + ImIndex(jj)), Level=10)
           END DO
         END IF
       ELSE
         CALL SimListAddAndOutputConstReal(&
           TRIM(Circuits(p) % names(i)), crt(Cvar % ValueId), Level=10)
         
         IF (Cvar % pdofs /= 0 ) THEN
           DO jj = 1, Cvar % pdofs
             CALL SimListAddAndOutputConstReal(&
               TRIM(Circuits(p) % names(i))&
               //'dof '//I2S(jj), crt(Cvar % ValueId + jj), Level=10)
           END DO
         END IF
       END IF

     END DO

     CALL Info(Caller, 'Writing Component Variables for &
       Circuit '//i2s(p), Level=8) 
     DO j = 1, SIZE(Circuits(p) % Components)
         Comp => Circuits(p) % Components(j)
         IF (Comp % Resistance < TINY(0._dp) .AND. Comp % Conductance > TINY(0._dp)) &
             Comp % Resistance = 1._dp / Comp % Conductance

         CALL SimListAddAndOutputConstReal('r_component('//&
           i2s(Comp % ComponentId)//')', Comp % Resistance, Level=8) 

         Current = 0._dp + im * 0._dp
         Current = crt(Comp % ivar % ValueId) 
         IF ( Circuits(p) % Harmonic ) Current = Current + im * crt(Comp % ivar % ImValueId) 
              
         CompParams => CurrentModel % Components (Comp % ComponentId) % Values
         IF (.NOT. ASSOCIATED(CompParams)) CALL Fatal ('CircuitsOutput', &
           'Component parameters not found!')

         p_dc_component = ABS(Current)**2._dp * Comp % Resistance
         CALL SimListAddAndOutputConstReal('p_dc_component('//i2s(Comp % ComponentId)//')',&
           p_dc_component, Level=8) 

         CompRealPower = GetConstReal( Model % Simulation, TRIM(CktPrefix)//' Power re & 
                 in Component '//i2s(Comp % ComponentId), Found)
         IF (Found .AND. ABS(Current) > TINY(CompRealPower)) THEN
           CALL SimListAddAndOutputConstReal('p_ac_component('//&
             i2s(Comp % ComponentId)//')', CompRealPower, Level=8)
           CALL SimListAddAndOutputConstReal('r_ac_component('//&
             i2s(Comp % ComponentId)//')', CompRealPower/ABS(Current)**2._dp, Level=8)
           CALL SimListAddAndOutputConstReal('AC to DC of component '&
             //i2s(Comp % ComponentId), CompRealPower/p_dc_component, Level=8)
         END IF
          
       END DO  
   END DO

   CALL Circuits_ToMeshVariable(Solver,crt)
   
   CALL DefaultFinish()


CONTAINS

  
!-------------------------------------------------------------------
  SUBROUTINE SimListAddAndOutputConstReal(VariableName, VariableValue, Level)
!-------------------------------------------------------------------
  IMPLICIT NONE
  CHARACTER(LEN=MAX_NAME_LEN) :: VarVal
  CHARACTER(LEN=*) :: VariableName
  REAL(KIND=dp) :: VariableValue
  INTEGER, OPTIONAL :: Level 
  INTEGER :: LevelVal = 3

  IF (PRESENT(Level)) LevelVal = Level
  WRITE(Message,'(A,T20,ES15.4)') TRIM(VariableName),VariableValue
  CALL Info(Caller,Message,Level=LevelVal)
  
  CALL ListAddConstReal(GetSimulation(),TRIM(CktPrefix)//' '//TRIM(VariableName), VariableValue)
!-------------------------------------------------------------------
  END SUBROUTINE SimListAddAndOutputConstReal
!-------------------------------------------------------------------

END SUBROUTINE CircuitsOutput
