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

  Params => Solver % Values

  ! When we introduce the variables in this way the variables are created
  ! so that they exist when the proper simulation cycle starts.
  ! This also keeps the command file cleaner.
  CALL ListAddString( Params,'Exported Variable 1',&
      '-global Rotor Angle')
  CALL ListAddString( Params,'Exported Variable 2',&
      '-global Rotor Velo')
  CALL ListAddLogical( Params,'No Matrix',.TRUE.)

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

  INTEGER :: p, n, istat, max_element_dofs
  TYPE(Mesh_t), POINTER :: Mesh  

  TYPE(Matrix_t), POINTER :: CM
  INTEGER, POINTER :: n_Circuits => Null()
  TYPE(Circuit_t), POINTER :: Circuits(:)
  
  REAL(KIND=dp), ALLOCATABLE, SAVE :: ip(:)     
  TYPE(Variable_t), POINTER :: LagrangeVar
  INTEGER, SAVE :: Tstep=-1
!------------------------------------------------------------------------------

  IF (First) THEN
    First = .FALSE.
    
    Model % HarmonicCircuits = .FALSE.
    CALL AddComponentsToBodyLists()
    
    ALLOCATE( Model%Circuit_tot_n, Model%n_Circuits, STAT=istat )
    IF ( istat /= 0 ) THEN
      CALL Fatal( 'CircuitsAndDynamicsHarmonic', 'Memory allocation error.' )
    END IF

    n_Circuits => Model%n_Circuits
    Model%Circuit_tot_n = 0
  
    Model % ASolver => FindSolverWithKey('Export Lagrange Multiplier', 26)
    ASolver => Model % ASolver
    
    CALL AllocateCircuitsList() ! CurrentModel%Circuits
    Circuits => Model%Circuits

    CALL SetBoundaryAreasToValueLists() 
    
    DO p=1,n_Circuits
      
      n = GetNofCircVariables(p)
      CALL AllocateCircuit(p)
      
      Circuits(p) % n_comp = CountNofCircComponents(p, n)
      ALLOCATE(Circuits(p) % Components(Circuits(p) % n_comp))
      
      Circuits(p) % Harmonic = .FALSE.
      
      CALL ReadCircuitVariables(p)
      CALL ReadComponents(p)
      CALL AddComponentValuesToLists(p)  ! Lists are used to communicate values to other solvers at the moment...
      CALL AddBareCircuitVariables(p)   ! these don't belong to any components
      CALL ReadCoefficientMatrices(p)
      CALL ReadPermutationVector(p)
      CALL ReadCircuitSources(p)
      CALL WriteCoeffVectorsForCircVariables(p)
    
    END DO

    ! Create CRS matrix strucures for the circuit equations:
    ! ------------------------------------------------------
    CALL Circuits_MatrixInit()
    ALLOCATE(ip(Model%Circuit_tot_n))
  END IF
  
  IF (Tstep /= GetTimestep()) THEN
    Tstep = GetTimestep()
    ! Circuit variable values from previous timestep:
    ! -----------------------------------------------
    ip = 0._dp
    LagrangeVar => VariableGet( Solver % Mesh % Variables,'LagrangeMultiplier')
    IF(ASSOCIATED(LagrangeVar)) THEN
      IF(SIZE(LagrangeVar % Values)>=Model%Circuit_tot_n) ip=LagrangeVar % Values(1:Model%Circuit_tot_n)
    END IF
  END IF

  max_element_dofs = Model % Mesh % MaxElementDOFs
  Circuits => Model%Circuits
  n_Circuits => Model%n_Circuits
  CM=>Model%CircuitMatrix
  
  ! Initialize Circuit matrix:
  ! -----------------------------
  IF(.NOT.ASSOCIATED(CM)) RETURN

  CM % RHS = 0._dp
  IF(ASSOCIATED(CM % Values)) CM % Values = 0._dp

  ! Write Circuit equations:
  ! ------------------------
  DO p = 1,n_Circuits
    CALL AddBasicCircuitEquations(p,ip,dt)
    CALL AddComponentEquationsAndCouplings(p, max_element_dofs,dt)
  END DO
  Asolver %  Matrix % AddMatrix => CM

  IF(ASSOCIATED(CM)) THEN
    IF(  CM % Format == MATRIX_LIST ) CALL List_toCRSMatrix(CM)
    IF(CM % NumberOfRows<=0)  THEN
      CALL FreeMatrix(CM)
      Asolver % Matrix % AddMatrix => Null()
    END IF
  ELSE
     ASolver % Matrix % AddMatrix => Null()
  END IF

  CONTAINS

!------------------------------------------------------------------------------
   SUBROUTINE AddBasicCircuitEquations(p,ip,dt)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(CircuitVariable_t), POINTER :: Cvar
    TYPE(ValueList_t), POINTER :: BF
    TYPE(Matrix_t), POINTER :: CM
    INTEGER :: p, i, nm, RowId, ColId, j
    REAL(KIND=dp) :: vphi,dt
    REAL(KIND=dp) :: ip(:)    
    LOGICAL :: Found
    
    Circuit => CurrentModel % Circuits(p)
    nm = CurrentModel % Asolver % Matrix % NumberOfRows
    BF => CurrentModel % BodyForces(1) % Values
    CM => CurrentModel%CircuitMatrix

    DO i=1,Circuit % n
      Cvar => Circuit % CircuitVariables(i)

      IF(Cvar % Owner /= ParEnv % myPE) CYCLE

      RowId = Cvar % ValueId + nm
      
      vphi=0._dp
      IF ( ASSOCIATED(BF) ) &
        vphi = GetCReal(BF, TRIM(Circuit % Source(i)), Found)
      
      IF (Found) THEN 
        Cvar % SourceRe(i) = vphi
      END IF
      
      CM % RHS(RowId) = Cvar % SourceRe(i)
        
      DO j=1,Circuit % n

        ColId = Circuit % CircuitVariables(j) % ValueId + nm


        IF ( TransientSimulation ) THEN 
          ! A d/dt(x): (x could be voltage or current):
          !--------------------------------------------
          IF(Cvar % A(j) /= 0._dp) THEN
            CALL AddToMatrixElement(CM, RowId, ColId, Cvar % A(j)/dt)
            CM % RHS(RowId) = CM % RHS(RowId) + Cvar % A(j)*ip(ColId-nm)/dt
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
   SUBROUTINE AddComponentEquationsAndCouplings(p, nn, dt)
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
    REAL(KIND=dp) :: value, dt
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilType
    LOGICAL :: Found

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

      CompParams => CurrentModel % Components(CompInd) % Values
      IF (.NOT. ASSOCIATED(CompParams)) CALL Fatal ('AddComponentEquationsAndCouplings', 'Component parameters not found')
      IF (Comp % CoilType == 'stranded') THEN
        Comp % Resistance = GetConstReal(CompParams, 'Resistance', Found)
        IF (Found) THEN
          Comp % UseCoilResistance = .TRUE.
        ELSE
          Comp % UseCoilResistance = .FALSE.
        END IF
      END IF

      IF ( Cvar % Owner == ParEnv % myPE ) THEN
        SELECT CASE (Comp % CoilType)
        CASE('stranded')
          IF (Comp % UseCoilResistance) THEN
            CALL Info('AddComponentEquationsAndCouplings', 'Using coil resistance for component '//TRIM(i2s(CompInd)), Level = 5)
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
            value = - REAL(Comp % nofturns) / REAL(j) * Comp % coilthickness**(j-1)
            CALL AddToMatrixElement(CM, VvarId, j + VvarId, value)

            ! Circuit eqns for the pdofs:
            ! - Nf/Lalpha * I * int_0^1(Vi'(alpha)) + ...added later... = 0
            ! ----------------------------------------------------------
            CALL AddToMatrixElement(CM, j + VvarId, IvarId, value)
          END DO
        END SELECT
      END IF
      
      DO q=GetNOFActive(),1,-1
        Element => GetActiveElement(q)
        IF (ElAssocToComp(Element, Comp)) THEN
          CompParams => GetComponentParams( Element )
          IF (.NOT. ASSOCIATED(CompParams)) CALL Fatal ('Circuits_apply', 'Component parameters not found')

          CoilType = GetString(CompParams, 'Coil Type', Found)
          IF (.NOT. Found) CoilType = ''
          
          nn = GetElementNOFNodes(Element)
          nd = GetElementNOFDOFs(Element,ASolver)
          !          CALL GetConductivity(Element, Tcoef, nn)
          
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
            CALL Add_stranded(Element,Tcoef,Comp,nn,nd,dt)
          CASE ('massive')
            IF (.NOT. HasSupport(Element,nn)) CYCLE
            CALL Add_massive(Element,Tcoef,Comp,nn,nd,dt)
          CASE ('foil winding')
            IF (.NOT. HasSupport(Element,nn)) CYCLE
            CALL Add_foil_winding(Element,Tcoef,Comp,nn,nd,dt)
          CASE DEFAULT
            CALL Fatal ('Circuits_apply', 'Non existent Coil Type Chosen!')
          END SELECT
        END IF
      END DO
    END DO

    DO CompInd = 1, Circuit % n_comp
      Comp => Circuit % Components(CompInd)
      Comp % Resistance = ParallelReduction(Comp % Resistance)
      Comp % Conductance = ParallelReduction(Comp % Conductance)
    END DO

    DEALLOCATE(Tcoef)
!------------------------------------------------------------------------------
   END SUBROUTINE AddComponentEquationsAndCouplings
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE Add_stranded(Element,Tcoef,Comp,nn,nd,dt)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t) :: Element
    REAL(KIND=dp) :: Tcoef(3,3,nn),dt
    TYPE(Component_t) :: Comp
    INTEGER :: nn, nd, nm, Indexes(nd),VvarId,IvarId

    TYPE(Solver_t), POINTER :: ASolver
    INTEGER, POINTER :: PS(:)
    TYPE(Matrix_t), POINTER :: CM
    TYPE(Nodes_t), SAVE :: Nodes
    REAL(KIND=dp) :: Basis(nn), DetJ, x,POT(nd),pPOT(nd),ppPOT(nd),tscl
    REAL(KIND=dp) :: dBasisdx(nn,3), wBase(nn), w(3)
    REAL(KIND=dp) :: localC, value, circ_eq_coeff, localR !, localL
    INTEGER :: j,t
    LOGICAL :: stat

    TYPE(GaussIntegrationPoints_t) :: IP
    LOGICAL :: CSymmetry, First=.TRUE.

    REAL(KIND=dp) :: WBasis(nd,3), RotWBasis(nd,3)
    INTEGER :: dim, ncdofs,q
    
    SAVE CSymmetry, dim
    
    IF (First) THEN
      First = .FALSE.
      CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )
      dim = CoordinateSystemDimension()
    END IF

    ASolver => CurrentModel % Asolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('Add_stranded','ASolver not found!')
    PS => Asolver % Variable % Perm

    CM => CurrentModel % CircuitMatrix
    nm = CurrentModel % Asolver % Matrix % NumberOfRows
    
    CALL GetElementNodes(Nodes)
    nd = GetElementDOFs(Indexes,Element,ASolver)

    CALL GetLocalSolution(pPOT,'A',Element,ASolver,tstep=-1)

    IF(Solver % Order<2.OR.GetTimeStep()<=2) THEN 
      tscl=1.0_dp
    ELSE
      tscl=1.5_dp
      CALL GetLocalSolution(ppPOT,'A',Element,ASolver,tstep=-2)
      pPot = 2*pPOT - 0.5_dp*ppPOT
    END IF
    
    ncdofs=nd
    IF (dim == 3) THEN
      CALL GetLocalSolution(Wbase, 'w')
      ncdofs=nd-nn
    END IF

    VvarId = Comp % vvar % ValueId + nm
    IvarId = Comp % ivar % ValueId + nm

    ! Numerical integration:
    ! ----------------------
    IP = GaussPoints(Element)
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                IP % W(t), detJ, Basis,dBasisdx )
      
      circ_eq_coeff = 1._dp
      SELECT CASE(dim)
      CASE(2)
        w = [0._dp, 0._dp, 1._dp]
        IF( CSymmetry ) THEN
          x = SUM( Basis(1:nn) * Nodes % x(1:nn) )
          detJ = detJ * x
        END IF
        circ_eq_coeff = GetCircuitModelDepth()
      CASE(3)
        CALL GetEdgeBasis(Element,WBasis,RotWBasis,Basis,dBasisdx)
        w = -MATMUL(WBase(1:nn), dBasisdx(1:nn,:))
      END SELECT

      localC = SUM(Tcoef(1,1,1:nn) * Basis(1:nn))
      
      IF (.NOT. Comp % UseCoilResistance) THEN
        ! I * R, where 
        ! R = (1/sigma * js,js):
        ! ----------------------
        localR = Comp % N_j **2 * IP % s(t)*detJ*SUM(w*w)/localC*circ_eq_coeff
        Comp % Resistance = Comp % Resistance + localR
      
        CALL AddToMatrixElement(CM, VvarId, IvarId, localR)
      END IF
      
      DO j=1,ncdofs
        q=j
        IF (dim == 3) q=q+nn
        IF (Comp % N_j/=0._dp) THEN
          ! ( d/dt a,w )        

          IF ( TransientSimulation ) THEN 
            IF (dim == 2) value = Comp % N_j * IP % s(t)*detJ*Basis(j)*circ_eq_coeff/dt*w(3)
            IF (dim == 3) value = Comp % N_j * IP % s(t)*detJ*SUM(WBasis(j,:)*w)/dt
 !          localL = value
!          Comp % Inductance = Comp % Inductance + localL

            CALL AddToMatrixElement(CM, VvarId, PS(Indexes(q)), tscl * value)
            CM % RHS(vvarid) = CM % RHS(vvarid) + pPOT(q) * value
         END IF
          
          ! source: 
          ! (J, rot a'), where
          ! J = w*I, thus I*(w, rot a'):
          ! ----------------------------         
          IF (dim == 2) value = -Comp % N_j*IP % s(t)*detJ*Basis(j)*w(3)
          IF (dim == 3) value = -Comp % N_j*IP % s(t)*detJ*SUM(WBasis(j,:)*w)
          CALL AddToMatrixElement(CM,PS(Indexes(q)), IvarId, value)

        END IF
      END DO
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE Add_stranded
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE Add_massive(Element,Tcoef,Comp,nn,nd,dt)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t) :: Element
    REAL(KIND=dp) :: Tcoef(3,3,nn),dt
    TYPE(Component_t) :: Comp

    TYPE(Solver_t), POINTER :: ASolver
    INTEGER, POINTER :: PS(:)
    TYPE(Matrix_t), POINTER :: CM
    REAL(KIND=dp) :: Basis(nn), DetJ, x,POT(nd),pPOT(nd),ppPOT(nd),tscl
    REAL(KIND=dp) :: dBasisdx(nn,3)
    REAL(KIND=dp) :: localC, value, circ_eq_coeff, grads_coeff, localConductance !, localL
    INTEGER :: nn, nd, j, t, nm, Indexes(nd), &
               VvarId, dim
    LOGICAL :: stat
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    LOGICAL :: CSymmetry, First=.TRUE.

    REAL(KIND=dp) :: wBase(nn), gradv(3), WBasis(nd,3), RotWBasis(nd,3)
    INTEGER :: ncdofs,q

    SAVE CSymmetry, dim

    IF (First) THEN
      First = .FALSE.
      CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )
      dim = CoordinateSystemDimension()
    END IF

    ASolver => CurrentModel % Asolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('Add_massive','ASolver not found!')
    PS => Asolver % Variable % Perm

    CM => CurrentModel % CircuitMatrix
    nm = CurrentModel % Asolver % Matrix % NumberOfRows

    CALL GetElementNodes(Nodes)
    nd = GetElementDOFs(Indexes,Element,ASolver)
    CALL GetLocalSolution(pPOT,'A',Element,ASolver,tstep=-1)

    IF(Solver % Order<2.OR.GetTimeStep()<=2) THEN 
      tscl=1.0_dp
    ELSE
      tscl=1.5_dp
      CALL GetLocalSolution(ppPOT,'A',Element,ASolver,tstep=-2)
      pPot = 2*pPOT - 0.5_dp*ppPOT
    END IF

    ncdofs=nd
    IF (dim == 3) THEN
      CALL GetLocalSolution(Wbase, 'w')
      ncdofs=nd-nn
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
        circ_eq_coeff = GetCircuitModelDepth()
        grads_coeff = grads_coeff/circ_eq_coeff
      CASE(3)
        CALL GetEdgeBasis(Element,WBasis,RotWBasis,Basis,dBasisdx)
        gradv = MATMUL( WBase(1:nn), dBasisdx(1:nn,:))
      END SELECT


      localC = SUM(Tcoef(1,1,1:nn) * Basis(1:nn))

      ! computing the source term Vi(sigma grad v0, grad si):
      ! ------------------------------------------------
      IF(dim==2) value = IP % s(t)*detJ*localC*grads_coeff**2*circ_eq_coeff
      IF(dim==3) value = IP % s(t)*detJ*localC*SUM(gradv*gradv)

      localConductance = ABS(value)
      Comp % Conductance = Comp % Conductance + localConductance

      CALL AddToMatrixElement(CM, vvarId, vvarId, value)

      DO j=1,ncdofs
        q=j
        IF (dim == 3) q=q+nn
        ! computing the mass term (sigma d/dt a, grad si):
        ! ---------------------------------------------------------
 
        IF ( TransientSimulation ) THEN 
          IF(dim==2) value = IP % s(t)*detJ*localC*basis(j)*grads_coeff*circ_eq_coeff/dt
          IF(dim==3) value = IP % s(t)*detJ*localC*SUM(Wbasis(j,:)*gradv)/dt
  !        localL = value
  !        Comp % Inductance = Comp % Inductance + localL
          CALL AddToMatrixElement(CM, vvarId, PS(Indexes(q)), tscl * value)
          CM % RHS(vvarid) = CM % RHS(vvarid) + pPOT(q) * value
        END IF

        IF(dim==2) value = IP % s(t)*detJ*localC*basis(j)*grads_coeff
        IF(dim==3) value = IP % s(t)*detJ*localC*SUM(gradv*Wbasis(j,:))
        CALL AddToMatrixElement(CM, PS(indexes(q)), vvarId, value)
      END DO
    END DO

!------------------------------------------------------------------------------
   END SUBROUTINE Add_massive
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE Add_foil_winding(Element,Tcoef,Comp,nn,nd,dt)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: nn, nd
    TYPE(Element_t) :: Element
    REAL(KIND=dp) :: Tcoef(3,3,nn), C(3,3), value, dt
    TYPE(Component_t) :: Comp

    TYPE(Solver_t), POINTER :: ASolver
    INTEGER, POINTER :: PS(:)
    TYPE(Matrix_t), POINTER :: CM
    REAL(KIND=dp) :: Basis(nn), DetJ, localAlpha, localV, localVtest, &
                     x, circ_eq_coeff, grads_coeff,POT(nd),pPOT(nd),ppPOT(nd),tscl
    REAL(KIND=dp) :: dBasisdx(nn,3),alpha(nn)
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

    SAVE CSymmetry, dim, First

    IF (First) THEN
      First = .FALSE.
      CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )
      dim = CoordinateSystemDimension()
    END IF

    ASolver => CurrentModel % Asolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('Add_foil_winding','ASolver not found!')
    PS => Asolver % Variable % Perm

    CM => CurrentModel % CircuitMatrix
    nm = CurrentModel % Asolver % Matrix % NumberOfRows

    CALL GetElementNodes(Nodes)
    nd = GetElementDOFs(Indexes,Element,ASolver)
    CALL GetLocalSolution(alpha,'Alpha')

    CALL GetLocalSolution(pPOT,'A',Element,ASolver,tstep=-1)

    IF(Solver % Order<2.OR.GetTimeStep()<=2) THEN 
      tscl=1.0_dp
    ELSE
      tscl=1.5_dp
      CALL GetLocalSolution(ppPOT,'A',Element,ASolver,tstep=-2)
      pPot = 2*pPOT - 0.5_dp*ppPOT
    END IF

    ncdofs=nd
    IF (dim == 3) THEN
      CALL GetLocalSolution(Wbase, 'w')
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

      ! I * R, where 
      ! R = (1/sigma * js,js):
      ! ----------------------
      IF (dim == 2) localR = Comp % N_j **2 * IP % s(t)*detJ/C(1,1)*circ_eq_coeff
      IF (dim == 3) localR = Comp % N_j **2 * IP % s(t)*detJ*SUM(gradv*MATMUL(C,gradv))
      Comp % Resistance = Comp % Resistance + localR

      DO vpolordtest=0,vpolord_tot ! V'(alpha)
        localVtest = localAlpha**vpolordtest
        dofIdtest = vpolordtest + 1 + vvarId
        DO vpolord = 0, vpolord_tot ! V(alpha)

          localV = localAlpha**vpolord
          dofId = vpolord + 1 + vvarId
          
          ! Computing the stiff term (sigma V(alpha) grad v0, V'(alpha) grad si):
          ! ---------------------------------------------------------------------
          IF (dim == 2) value = IP % s(t)*detJ*localV*localVtest*C(1,1)*grads_coeff**2*circ_eq_coeff
          IF (dim == 3) value = IP % s(t)*detJ*localV*localVtest*SUM(MATMUL(C,gradv)*gradv)
          CALL AddToMatrixElement(CM, dofIdtest+nm, dofId+nm, value)
        END DO


        IF ( TransientSimulation ) THEN 
          DO j=1,ncdofs
            q=j
            IF (dim == 3) q=q+nn
            ! computing the mass term (sigma * d/dt * a, V'(alpha) grad si):
            ! ---------------------------------------------------------
            IF (dim == 2) value = IP % s(t)*detJ*localVtest*C(1,1)*basis(j)*grads_coeff*circ_eq_coeff/dt
            IF (dim == 3) value = IP % s(t)*detJ*localVtest*SUM(MATMUL(C,Wbasis(j,:))*gradv)/dt
  !          localL = value
  !          Comp % Inductance = Comp % Inductance + localL
            CALL AddToMatrixElement(CM, dofIdtest+nm, PS(Indexes(q)), tscl * value)
            CM % RHS(dofIdtest+nm) = CM % RHS(dofIdtest+nm) + pPOT(q) * value
          END DO
        END IF

      END DO

      DO vpolord = 0, vpolord_tot ! V(alpha)
        localV = localAlpha**vpolord
        dofId = vpolord + 1 + vvarId

        DO j=1,ncdofs
            q=j
            IF (dim == 3) q=q+nn
            IF (dim == 2) value = IP % s(t)*detJ*localV*C(1,1)*basis(j)*grads_coeff
            IF (dim == 3) value = IP % s(t)*detJ*localV*SUM(MATMUL(C,gradv)*Wbasis(j,:))
            CALL AddToMatrixElement(CM, PS(indexes(q)), dofId+nm, value)
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
 SUBROUTINE GetElementRotM(Element,RotM,n)
!------------------------------------------------------------------------------
   IMPLICIT NONE
   TYPE(Element_t) :: Element
   INTEGER :: k, l, m, j, n
   REAL(KIND=dp) :: RotM(3,3,n)
   TYPE(Variable_t), POINTER, SAVE :: RotMvar
   LOGICAL, SAVE :: visited = .FALSE.
   INTEGER, PARAMETER :: ind1(9) = [1,1,1,2,2,2,3,3,3]
   INTEGER, PARAMETER :: ind2(9) = [1,2,3,1,2,3,1,2,3]
  
   IF(.NOT. visited) THEN
     visited = .TRUE.
     RotMvar => VariableGet( CurrentModel % Mesh % Variables, 'RotM E')
     IF(.NOT. ASSOCIATED(RotMVar)) THEN
       CALL Fatal('GetElementRotM','RotM E variable not found')
     END IF
   END IF

   RotM = 0._dp
   DO j = 1, n
     DO k=1,RotMvar % DOFs
       RotM(ind1(k),ind2(k),j) = RotMvar % Values( &
             RotMvar % DOFs*(RotMvar % Perm(Element % DGIndexes(j))-1)+k)
     END DO
   END DO

!------------------------------------------------------------------------------
 END SUBROUTINE GetElementRotM
!------------------------------------------------------------------------------


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

  ! When we introduce the variables in this way the variables are created
  ! so that they exist when the proper simulation cycle starts.
  ! This also keeps the command file cleaner.
  CALL ListAddString( Params,'Exported Variable 1',&
      '-global Rotor Angle')
  CALL ListAddString( Params,'Exported Variable 2',&
      '-global Rotor Velo')
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

  INTEGER :: p, n, istat, max_element_dofs
  TYPE(Mesh_t), POINTER :: Mesh  

  TYPE(Matrix_t), POINTER :: CM
  INTEGER, POINTER :: n_Circuits => Null(), circuit_tot_n => Null()
  TYPE(Circuit_t), POINTER :: Circuits(:)
    
!------------------------------------------------------------------------------

  CALL DefaultStart()

  IF (First) THEN
    First = .FALSE.
    
    Model % HarmonicCircuits = .TRUE.

    CALL AddComponentsToBodyLists()
    
    ALLOCATE( Model%Circuit_tot_n, Model%n_Circuits, STAT=istat )
    IF ( istat /= 0 ) THEN
      CALL Fatal( 'CircuitsAndDynamicsHarmonic', 'Memory allocation error.' )
    END IF

    n_Circuits => Model%n_Circuits
    Model%Circuit_tot_n = 0
  
    Model % ASolver => FindSolverWithKey('Export Lagrange Multiplier', 26)
    ASolver => Model % ASolver
    
    CALL AllocateCircuitsList() ! CurrentModel%Circuits
    Circuits => Model%Circuits

    CALL SetBoundaryAreasToValueLists() 
    
    DO p=1,n_Circuits
      
      n = GetNofCircVariables(p)
      CALL AllocateCircuit(p)
      
      Circuits(p) % n_comp = CountNofCircComponents(p, n)
      ALLOCATE(Circuits(p) % Components(Circuits(p) % n_comp))
      
      Circuits(p) % Harmonic = .TRUE.
      
      CALL ReadCircuitVariables(p)
      CALL ReadComponents(p)
      CALL AddComponentValuesToLists(p)  ! Lists are used to communicate values to other solvers at the moment...
      CALL AddBareCircuitVariables(p)   ! these don't belong to any components
      CALL ReadCoefficientMatrices(p)
      CALL ReadPermutationVector(p)
      CALL ReadCircuitSources(p)
      CALL WriteCoeffVectorsForCircVariables(p)
    
    END DO

    ! Create CRS matrix strucures for the circuit equations:
    ! ------------------------------------------------------
    CALL Circuits_MatrixInit()
  END IF
  
  max_element_dofs = Model % Mesh % MaxElementDOFs
  Circuits => Model%Circuits
  n_Circuits => Model%n_Circuits
  CM=>Model%CircuitMatrix
  
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
  Asolver %  Matrix % AddMatrix => CM

  IF(ASSOCIATED(CM)) THEN
    IF(  CM % Format == MATRIX_LIST ) CALL List_toCRSMatrix(CM)
    IF(CM % NumberOfRows<=0)  THEN
      CALL FreeMatrix(CM)
      Asolver % Matrix % AddMatrix => Null()
    END IF
  ELSE
     ASolver % Matrix % AddMatrix => Null()
  END IF

  CALL DefaultFinish()

  CONTAINS

!------------------------------------------------------------------------------
   SUBROUTINE AddBasicCircuitEquations(p)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(CircuitVariable_t), POINTER :: Cvar
    TYPE(ValueList_t), POINTER :: BF
    TYPE(Matrix_t), POINTER :: CM
    INTEGER :: p, i, nm, RowId, ColId, j
    REAL(KIND=dp) :: Omega, vphi
    COMPLEX(KIND=dp) :: cmplx_value
    LOGICAL :: Found
    COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
    
    Circuit => CurrentModel % Circuits(p)
    nm = CurrentModel % Asolver % Matrix % NumberOfRows
    Omega = GetAngularFrequency()
    BF => CurrentModel % BodyForces(1) % Values
    CM => CurrentModel%CircuitMatrix
    
    DO i=1,Circuit % n
      Cvar => Circuit % CircuitVariables(i)

      IF(Cvar % Owner /= ParEnv % myPE) CYCLE

      RowId = Cvar % ValueId + nm
      
      vphi=0._dp
      IF ( ASSOCIATED(BF) ) &
        vphi = GetCReal(BF, TRIM(Circuit % Source(i))//" re", Found)
      
      IF (Found) THEN 
        Cvar % SourceRe(i) = vphi
      END IF
      
      IF (ASSOCIATED(BF) ) &
        vphi = GetCReal(BF, TRIM(Circuit % Source(i))//" im", Found)

      IF (Found) THEN 
        Cvar % SourceIm(i) = vphi
      END IF
      
      CM % RHS(RowId) = Cvar % SourceRe(i)
      CM % RHS(RowId+1) = Cvar % SourceIm(i)
        
      DO j=1,Circuit % n

        ColId = Circuit % CircuitVariables(j) % ValueId + nm

        ! im * Omega * A x: (x could be voltage or current):
        !--------------------------------------------
        IF(Cvar % A(j) /= 0._dp) THEN
          CALL AddToCmplxMatrixElement(CM, RowId, ColId, 0._dp, Omega * Cvar % A(j))
        END IF
        ! B x:
        ! ------
        IF(Cvar % B(j) /= 0._dp) THEN
          
          IF (Cvar % Mre(j) /= 0._dp .OR. Cvar % Mim(j) /= 0._dp) THEN
            cmplx_value = Cvar % Mre(j) + im * Cvar % Mim(j)
            cmplx_value = cmplx_value * Cvar % B(j)
          ELSE
            cmplx_value = Cvar % B(j)
          END IF
          
          CALL AddToCmplxMatrixElement(CM, RowId, ColId, REAL(cmplx_value), AIMAG(cmplx_value))
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
    INTEGER :: CompInd, nm, nn_elem, nd_elem
    TYPE(Solver_t), POINTER :: ASolver
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(Matrix_t), POINTER :: CM
    TYPE(Component_t), POINTER :: Comp
    TYPE(CircuitVariable_t), POINTER :: Cvar
    TYPE(Valuelist_t), POINTER :: CompParams
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: Omega
    REAL(KIND=dp), ALLOCATABLE :: sigma_33(:), sigmaim_33(:)
    INTEGER :: VvarId, IvarId, q, j, astat
    COMPLEX(KIND=dp) :: i_multiplier, cmplx_value
    COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
    COMPLEX(KIND=dp), ALLOCATABLE :: Tcoef(:,:,:)
    REAL(KIND=dp) :: RotM(3,3,nn)
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilType
    LOGICAL :: Found, FoundIm, StrandedHomogenization

    ASolver => CurrentModel % Asolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('AddComponentEquationsAndCouplings','ASolver not found!')

    Circuit => CurrentModel % Circuits(p)
    nm = Asolver % Matrix % NumberOfRows
    Omega = GetAngularFrequency()
    CM => CurrentModel%CircuitMatrix

    DO CompInd = 1, Circuit % n_comp
      Comp => Circuit % Components(CompInd)

      Comp % Resistance = 0._dp 
      Comp % Conductance = 0._dp 

      Cvar => Comp % vvar
      vvarId = Comp % vvar % ValueId + nm
      IvarId = Comp % ivar % ValueId + nm

      CompParams => CurrentModel % Components(CompInd) % Values
      IF (.NOT. ASSOCIATED(CompParams)) CALL Fatal ('AddComponentEquationsAndCouplings', 'Component parameters not found')
      IF (Comp % CoilType == 'stranded') THEN
        Comp % Resistance = GetConstReal(CompParams, 'Resistance', Found)
        IF (Found) THEN
          Comp % UseCoilResistance = .TRUE.
        ELSE
          Comp % UseCoilResistance = .FALSE.
        END IF
      END IF

      IF ( Cvar % Owner == ParEnv % myPE ) THEN
        SELECT CASE (Comp % CoilType)
        CASE('stranded')
          IF (Comp % UseCoilResistance) THEN
            CALL Info('AddComponentEquationsAndCouplings', 'Using coil resistance for component '//TRIM(i2s(CompInd)), Level = 5)
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
            cmplx_value = -i_multiplier * REAL(Comp % nofturns) / REAL(j) * Comp % coilthickness**(j-1)
            CALL AddToCmplxMatrixElement(CM, VvarId, 2*j + VvarId, &
                REAL(cmplx_value), AIMAG(cmplx_value))

            ! Circuit eqns for the pdofs:
            ! - Nf/Lalpha * I * int_0^1(Vi'(alpha)) + ...added later... = 0
            ! ----------------------------------------------------------
            CALL AddToCmplxMatrixElement(CM, 2*j + VvarId, IvarId, &
               REAL(cmplx_value), AIMAG(cmplx_value))
          END DO
        END SELECT
      END IF

      DO q=GetNOFActive(),1,-1
        Element => GetActiveElement(q)
        IF (ElAssocToComp(Element, Comp)) THEN
          CompParams => GetComponentParams( Element )
          IF (.NOT. ASSOCIATED(CompParams)) CALL Fatal ('AddComponentEquationsAndCouplings',&
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
              sigma_33 = 0._dp
              sigmaim_33 = 0._dp
              sigma_33 = GetReal(CompParams, 'sigma 33', Found)
              sigmaim_33 = GetReal(CompParams, 'sigma 33 im', FoundIm)
              IF ( .NOT. Found .AND. .NOT. FoundIm ) CALL Fatal ('AddComponentEquationsAndCouplings', &
                                                                 'Homogenization Model Sigma 33 not found!')
              Tcoef = CMPLX(0._dp, 0._dp, KIND=dp)
              Tcoef(3,3,1:nn_elem) = CMPLX(sigma_33, sigmaim_33, KIND=dp)
            ELSE
              Tcoef = GetCMPLXElectricConductivityTensor(Element, nn_elem, .TRUE., CoilType) 
            END IF
            CALL Add_stranded(Element,Tcoef,Comp,nn_elem,nd_elem)
          CASE ('massive')
            IF (.NOT. HasSupport(Element,nn_elem)) CYCLE
         !   CALL GetConductivity(Element, Tcoef, nn_elem)
            Tcoef = GetCMPLXElectricConductivityTensor(Element, nn_elem, .TRUE., CoilType) 
            CALL Add_massive(Element,Tcoef,Comp,nn_elem,nd_elem)
          CASE ('foil winding')
            IF (.NOT. HasSupport(Element,nn_elem)) CYCLE
         !   CALL GetConductivity(Element, Tcoef, nn_elem)
            Tcoef = GetCMPLXElectricConductivityTensor(Element, nn_elem, .TRUE., CoilType) 
            CALL Add_foil_winding(Element,Tcoef,Comp,nn_elem,nd_elem)
          CASE DEFAULT
            CALL Fatal ('AddComponentEquationsAndCouplings', 'Non existent Coil Type Chosen!')
          END SELECT
        END IF
      END DO
    END DO

    DO CompInd = 1, Circuit % n_comp
      Comp => Circuit % Components(CompInd)
      Comp % Resistance = ParallelReduction(Comp % Resistance)
      Comp % Conductance = ParallelReduction(Comp % Conductance)
    END DO

    IF (ALLOCATED(Tcoef)) THEN
      DEALLOCATE(Tcoef,sigma_33,sigmaim_33)
    END IF
!------------------------------------------------------------------------------
   END SUBROUTINE AddComponentEquationsAndCouplings
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE Add_stranded(Element,Tcoef,Comp,nn,nd)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t) :: Element
    COMPLEX(KIND=dp) :: Tcoef(3,3,nn)
    TYPE(Component_t) :: Comp
    INTEGER :: nn, nd, nm, Indexes(nd),VvarId,IvarId

    TYPE(Solver_t), POINTER :: ASolver
    INTEGER, POINTER :: PS(:)
    TYPE(Matrix_t), POINTER :: CM
    REAL(KIND=dp) :: Omega
    TYPE(Nodes_t), SAVE :: Nodes
    REAL(KIND=dp) :: Basis(nn), DetJ, x, circ_eq_coeff
    REAL(KIND=dp) :: dBasisdx(nn,3), wBase(nn), w(3)
    COMPLEX(KIND=dp) :: localC, i_multiplier, cmplx_value
    REAL(KIND=dp) :: localR !, localL
    INTEGER :: j,t
    LOGICAL :: stat

    TYPE(GaussIntegrationPoints_t) :: IP
    COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
    LOGICAL :: CSymmetry, First=.TRUE.

    REAL(KIND=dp) :: WBasis(nd,3), RotWBasis(nd,3)
    INTEGER :: dim, ncdofs,q
    
    SAVE CSymmetry, dim

    IF (First) THEN
      First = .FALSE.
      CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )
      dim = CoordinateSystemDimension()
    END IF

    ASolver => CurrentModel % Asolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('Add_stranded','ASolver not found!')
    PS => Asolver % Variable % Perm

    CM => CurrentModel % CircuitMatrix
    nm = CurrentModel % Asolver % Matrix % NumberOfRows
    Omega = GetAngularFrequency()
    
    CALL GetElementNodes(Nodes)
    nd = GetElementDOFs(Indexes,Element,ASolver)
    
    ncdofs=nd
    IF (dim == 3) THEN
      ncdofs=nd-nn
      CALL GetWPotential(WBase)
      !print *, "W Potential", Wbase
    END IF

    VvarId = Comp % vvar % ValueId + nm
    IvarId = Comp % ivar % ValueId + nm

    i_multiplier = Comp % i_multiplier_re + im * Comp % i_multiplier_im


    ! Numerical integration:
    ! ----------------------
    IP = GaussPoints(Element)
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                IP % W(t), detJ, Basis,dBasisdx )

      circ_eq_coeff = 1._dp
      SELECT CASE(dim)
      CASE(2)
        w = [0._dp, 0._dp, 1._dp]
        IF( CSymmetry ) THEN
          x = SUM( Basis(1:nn) * Nodes % x(1:nn) )
          detJ = detJ * x
        END IF
        circ_eq_coeff = GetCircuitModelDepth()
      CASE(3)
        CALL GetEdgeBasis(Element,WBasis,RotWBasis,Basis,dBasisdx)
        w = -MATMUL(WBase(1:nn), dBasisdx(1:nn,:))
        !print *, "W Pot norm:", SQRT(SUM(w**2._dp))
      END SELECT

      localC = SUM(Tcoef(3,3,1:nn) * Basis(1:nn))
      
      IF (.NOT. Comp % UseCoilResistance) THEN
        ! I * R, where 
        ! R = (1/sigma * js,js):
        ! ----------------------
        localR = Comp % N_j **2 * IP % s(t)*detJ*SUM(w*w)/localC*circ_eq_coeff
        Comp % Resistance = Comp % Resistance + localR
        
        CALL AddToCmplxMatrixElement(CM, VvarId, IvarId, &
              REAL(Comp % N_j**2 * IP % s(t)*detJ*SUM(w*w)/localC*circ_eq_coeff), &
             AIMAG(Comp % N_j**2 * IP % s(t)*detJ*SUM(w*w)/localC*circ_eq_coeff))
      END IF
      
      DO j=1,ncdofs
        q=j
        IF (dim == 3) q=q+nn
        IF (Comp % N_j/=0._dp) THEN
          ! ( im * Omega a,w )
          IF (dim == 2) cmplx_value = im * Omega * Comp % N_j &
                  * IP % s(t)*detJ*Basis(j)*circ_eq_coeff*w(3)
          IF (dim == 3) cmplx_value = im * Omega * Comp % N_j &
                  * IP % s(t)*detJ*SUM(WBasis(j,:)*w)

!          localL = ABS(cmplx_value)
!          Comp % Inductance = Comp % Inductance + localL

          CALL AddToCmplxMatrixElement(CM, VvarId, ReIndex(PS(Indexes(q))), &
                 REAL(cmplx_value), AIMAG(cmplx_value))
          
          IF (dim == 2) cmplx_value = -Comp % N_j*IP % s(t)*detJ*Basis(j)*w(3)
          IF (dim == 3) cmplx_value = -Comp % N_j*IP % s(t)*detJ*SUM(WBasis(j,:)*w)
          IF (i_multiplier /= 0._dp) cmplx_value = i_multiplier*cmplx_value
          
          CALL AddToCmplxMatrixElement(CM,ReIndex(PS(Indexes(q))), IvarId, &
             REAL(cmplx_value), AIMAG(cmplx_value))

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
    INTEGER, POINTER :: PS(:)
    TYPE(Matrix_t), POINTER :: CM
    REAL(KIND=dp) :: Omega, grads_coeff, circ_eq_coeff
    REAL(KIND=dp) :: Basis(nn), DetJ, x
    REAL(KIND=dp) :: dBasisdx(nn,3)
    COMPLEX(KIND=dp) :: localC, cmplx_value
    REAL(KIND=dp) :: localConductance !, localL
    INTEGER :: nn, nd, j, t, nm, Indexes(nd), &
               VvarId, dim
    LOGICAL :: stat
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
    LOGICAL :: CSymmetry, First=.TRUE.

    REAL(KIND=dp) :: wBase(nn), gradv(3), WBasis(nd,3), RotWBasis(nd,3)
    INTEGER :: ncdofs,q
    REAL(KIND=dp) :: ModelDepth

    SAVE CSymmetry, dim

    IF (First) THEN
      First = .FALSE.
      CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )
      dim = CoordinateSystemDimension()
    END IF

    ASolver => CurrentModel % Asolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('Add_massive','ASolver not found!')
    PS => Asolver % Variable % Perm

    CM => CurrentModel % CircuitMatrix
    nm = CurrentModel % Asolver % Matrix % NumberOfRows
    Omega = GetAngularFrequency()

    CALL GetElementNodes(Nodes)
    nd = GetElementDOFs(Indexes,Element,ASolver)

    ncdofs=nd
    IF (dim == 3) THEN
      CALL GetWPotential(WBase)
      ncdofs=nd-nn
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

      ! computing the source term Vi(sigma grad v0, grad si):
      ! ------------------------------------------------
      IF(dim==2) cmplx_value = IP % s(t)*detJ*localC*grads_coeff**2*circ_eq_coeff
      IF(dim==3) cmplx_value = IP % s(t)*detJ*localC*SUM(gradv*gradv)
      localConductance = ABS(cmplx_value)
      Comp % Conductance = Comp % Conductance + localConductance
      CALL AddToCmplxMatrixElement(CM, vvarId, vvarId, &
              REAL(cmplx_value), AIMAG(cmplx_value))

      DO j=1,ncdofs
        q=j
        IF (dim == 3) q=q+nn
        ! computing the mass term (sigma * im * Omega * a, grad si):
        ! ---------------------------------------------------------
        IF(dim==2) cmplx_value = im * Omega * IP % s(t)*detJ*localC*basis(j)*grads_coeff*circ_eq_coeff
        IF(dim==3) cmplx_value = im * Omega * IP % s(t)*detJ*localC*SUM(Wbasis(j,:)*gradv)
        CALL AddToCmplxMatrixElement(CM, vvarId, ReIndex(PS(Indexes(q))), &
               REAL(cmplx_value), AIMAG(cmplx_value))

!        localL = ABS(1._dp/cmplx_value)
!        Comp % Inductance = Comp % Inductance + localL
        
        IF(dim==2) cmplx_value = IP % s(t)*detJ*localC*basis(j)*grads_coeff
        IF(dim==3) cmplx_value = IP % s(t)*detJ*localC*SUM(gradv*Wbasis(j,:))
        CALL AddToCmplxMatrixElement(CM, ReIndex(PS(indexes(q))), vvarId, &
                REAL(cmplx_value), AIMAG(cmplx_value))
      END DO
    END DO

!------------------------------------------------------------------------------
   END SUBROUTINE Add_massive
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE Add_foil_winding(Element,Tcoef,Comp,nn,nd)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: nn, nd
    TYPE(Element_t) :: Element
    COMPLEX(KIND=dp) :: Tcoef(3,3,nn), C(3,3), value
    TYPE(Component_t) :: Comp

    TYPE(Solver_t), POINTER :: ASolver
    INTEGER, POINTER :: PS(:)
    TYPE(Matrix_t), POINTER :: CM
    REAL(KIND=dp) :: Basis(nn), DetJ, Omega, localAlpha, localV, localVtest, &
                     x, circ_eq_coeff, grads_coeff
    REAL(KIND=dp) :: dBasisdx(nn,3),alpha(nn)
    INTEGER :: nm,p,j,t,Indexes(nd),vvarId,vpolord_tot, &
               vpolord, vpolordtest, dofId, dofIdtest, &
               dim
    LOGICAL :: stat
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
    LOGICAL :: CSymmetry, First=.TRUE.
    REAL(KIND=dp) :: localR

    REAL(KIND=dp) :: wBase(nn), gradv(3), WBasis(nd,3), RotWBasis(nd,3), &
                     RotMLoc(3,3), RotM(3,3,nn)
    INTEGER :: i,ncdofs,q

    SAVE CSymmetry, dim, First

    IF (First) THEN
      First = .FALSE.
      CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )
      dim = CoordinateSystemDimension()
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
      !CALL GetLocalSolution(Wbase, 'w')
      CALL GetWPotential(WBase)
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

      ! I * R, where 
      ! R = (1/sigma * js,js):
      ! ----------------------
      IF (dim == 2) localR = Comp % N_j **2 * IP % s(t)*detJ/C(1,1)*circ_eq_coeff
      IF (dim == 3) localR = Comp % N_j **2 * IP % s(t)*detJ*SUM(gradv*MATMUL(C,gradv))
      Comp % Resistance = Comp % Resistance + localR

      DO vpolordtest=0,vpolord_tot ! V'(alpha)
        localVtest = localAlpha**vpolordtest
        dofIdtest = 2*(vpolordtest + 1) + vvarId
        DO vpolord = 0, vpolord_tot ! V(alpha)

          localV = localAlpha**vpolord
          dofId = 2*(vpolord + 1) + vvarId
          
          ! Computing the stiff term (sigma V(alpha) grad v0, V'(alpha) grad si):
          ! ---------------------------------------------------------------------
          IF (dim == 2) value = IP % s(t)*detJ*localV*localVtest*C(1,1)*grads_coeff**2*circ_eq_coeff
          IF (dim == 3) value = IP % s(t)*detJ*localV*localVtest*SUM(MATMUL(C,gradv)*gradv)

          CALL AddToCmplxMatrixElement(CM, dofIdtest+nm, dofId+nm, REAL(value), AIMAG(value))
        END DO

        DO j=1,ncdofs
          q=j
          IF (dim == 3) q=q+nn
          ! computing the mass term (sigma * im * Omega * a, V'(alpha) grad si):
          ! ---------------------------------------------------------
          IF (dim == 2) value = im * Omega * IP % s(t)*detJ*localVtest*C(1,1)*basis(j)*grads_coeff*circ_eq_coeff
          IF (dim == 3) value = im * Omega * IP % s(t)*detJ*localVtest*SUM(MATMUL(C,Wbasis(j,:))*gradv)
          CALL AddToCmplxMatrixElement(CM, dofIdtest+nm, ReIndex(PS(Indexes(q))), REAL(value), AIMAG(value) )
        END DO

      END DO

      DO vpolord = 0, vpolord_tot ! V(alpha)
        localV = localAlpha**vpolord
        dofId = 2*(vpolord + 1) + vvarId

        DO j=1,ncdofs
            q=j
            IF (dim == 3) q=q+nn
            IF (dim == 2) value = IP % s(t)*detJ*localV*C(1,1)*basis(j)*grads_coeff
            IF (dim == 3) value = IP % s(t)*detJ*localV*SUM(MATMUL(C,gradv)*Wbasis(j,:))
            CALL AddToCmplxMatrixElement(CM, ReIndex(PS(indexes(q))), dofId+nm, REAL(value), AIMAG(value))
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
 SUBROUTINE GetElementRotM(Element,RotM,n)
!------------------------------------------------------------------------------
   IMPLICIT NONE
   TYPE(Element_t) :: Element
   INTEGER :: k, l, m, j, n
   REAL(KIND=dp) :: RotM(3,3,n)
   TYPE(Variable_t), POINTER, SAVE :: RotMvar
   LOGICAL, SAVE :: visited = .FALSE.
   INTEGER, PARAMETER :: ind1(9) = [1,1,1,2,2,2,3,3,3]
   INTEGER, PARAMETER :: ind2(9) = [1,2,3,1,2,3,1,2,3]
  
   IF(.NOT. visited) THEN
     visited = .TRUE.
     RotMvar => VariableGet( CurrentModel % Mesh % Variables, 'RotM E')
     IF(.NOT. ASSOCIATED(RotMVar)) THEN
       CALL Fatal('GetElementRotM','RotM E variable not found')
     END IF
   END IF

   RotM = 0._dp
   DO j = 1, n
     DO k=1,RotMvar % DOFs
       RotM(ind1(k),ind2(k),j) = RotMvar % Values( &
             RotMvar % DOFs*(RotMvar % Perm(Element % DGIndexes(j))-1)+k)
     END DO
   END DO

!------------------------------------------------------------------------------
 END SUBROUTINE GetElementRotM
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
END SUBROUTINE CircuitsAndDynamicsHarmonic
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
SUBROUTINE CircuitsOutput(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
   USE DefUtils
   USE CircuitsMod
   IMPLICIT NONE
!------------------------------------------------------------------------------   
   TYPE(Model_t) :: Model
   TYPE(Solver_t) :: Solver
   REAL(KIND=dp) :: dt
   LOGICAL :: Transient
!------------------------------------------------------------------------------

   TYPE(Variable_t), POINTER :: LagrangeVar
   REAL(KIND=dp), ALLOCATABLE  :: ip(:), ipt(:)
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
  REAL, SAVE :: EEC_freq, EEC_time_0
  INTEGER, SAVE :: EEC_max, EEC_cnt = 0
  REAL :: TTime
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Variable_t), POINTER :: AzVar
  REAL (KIND=dp), POINTER :: Az(:)
  REAL (KIND=dp), ALLOCATABLE, SAVE :: Az0(:)
  REAL (KIND=dp), POINTER :: Acorr(:)
  
!------------------------------------------------------------------------------  
  
   CALL DefaultStart()

   Circuit_tot_n => Model%Circuit_tot_n
   n_Circuits => Model%n_Circuits
   CM => Model%CircuitMatrix
   Circuits => Model%Circuits

   ! Look for the solver we attach the circuit equations to:
   ! -------------------------------------------------------
   ASolver => CurrentModel % Asolver
   IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('CircuitsOutput','ASolver not found!')
      
   nm =  Asolver % Matrix % NumberOfRows


  IF (First) THEN
    SolverParams => GetSolverParams(Solver)
    ! Reading parameter for supply frequency
    EEC_freq = GetConstReal( SolverParams, 'EEC Frequency', EEC)
    IF (EEC) THEN
      CALL Info('CircuitsAndDynamicsEEC', "Using EEC steady state forcing.", Level=1)
	    WRITE( Message,'(A,4G11.4,A)') 'EEC signal frequency: ', EEC_freq, ' Hz'
      CALL Info('CircuitsAndDynamicsEEC', Message, Level=1)
      
          
      EEC_max = GetInteger( SolverParams, 'EEC Steps', EEC_lim)
      IF (.NOT. EEC_lim) EEC_max = 5 !Typically 5 correections is enough
      WRITE( Message,'(A,I5,A)') 'Applying ', EEC_max, ' halfperiod corrections'
      CALL Info('CircuitsAndDynamicsEEC', Message, Level=1)
      
      EEC_time_0 = 0.0
      
      ! Reserve memory for storing current MVP solution
      ALLOCATE(Az0(nm))
      
      ! Store MVP solution at t=0
      AzVar => VariableGet( Asolver % Mesh % Variables, 'A')
      IF(ASSOCIATED( AzVar)) THEN
        Az => AzVar % Values
        Az0 = Az
      END IF
      
    END IF
    First = .FALSE.
  END IF


  IF (EEC .AND. (EEC_cnt .LT. EEC_max)) THEN
    
    TTime = GetTime()
    IF(TTime .GE. (EEC_time_0 + 0.5/EEC_freq)) THEN
      EEC_cnt = EEC_cnt + 1
      WRITE( Message,'(A,4G11.4)') 'Performing EEC #', EEC_cnt
      CALL Info('CircuitsAndDynamicsEEC', Message, Level=1)
      
      EEC_time_0 = EEC_time_0 + 0.5/EEC_freq

      AzVar => VariableGet( Asolver % Mesh % Variables, 'A')
      
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
   ALLOCATE(ip(circuit_tot_n), ipt(circuit_tot_n))
   ip = 0._dp
   ipt = 0._dp
   LagrangeVar => VariableGet( Solver % Mesh % Variables,'LagrangeMultiplier')
   IF(ASSOCIATED(LagrangeVar)) THEN
     IF(ParEnv % PEs>1) THEN
       DO i=1,circuit_tot_n 
         IF (ASSOCIATED(Model%CircuitMatrix)) THEN  
           IF( CM % RowOwner(nm+i)==Parenv%myPE) ipt(i) = LagrangeVar%Values(i)
         END IF
       END DO
       CALL MPI_ALLREDUCE(ipt,ip,circuit_tot_n, MPI_DOUBLE_PRECISION, &
                  MPI_SUM, ASolver % Matrix % Comm, j)
     ELSE
       ip(1:circuit_tot_n) = LagrangeVar % Values
     END IF
   END IF
    
   ! Export circuit & dynamic variables for "SaveScalars":
   ! -----------------------------------------------------

   CALL ListAddConstReal(GetSimulation(),'res: time', GetTime())

   CALL Info('CircuitsOutput', 'Writing Circuit Results', Level=3) 
   DO p=1,n_Circuits
     CALL Info('CircuitsOutput', 'Writing Circuit Variables for &
       Circuit '//TRIM(i2s(p)), Level=3) 
     CALL Info('CircuitsOutput', 'There are '//TRIM(i2s(Circuits(p)%n))//&
       ' Circuit Variables', Level=3)
     DO i=1,Circuits(p) % n
       Cvar => Circuits(p) % CircuitVariables(i)
       
       IF (Circuits(p) % Harmonic) THEN 
         CALL SimListAddAndOutputConstReal(&
           TRIM(Circuits(p) % names(i))//' re', ip(Cvar % ValueId), Level=10)
         CALL SimListAddAndOutputConstReal(&
           TRIM(Circuits(p) % names(i))//' im', ip(Cvar % ImValueId), Level=10)

         IF (Cvar % pdofs /= 0 ) THEN
           DO jj = 1, Cvar % pdofs
             write (dofnumber, "(I2)") jj
             CALL SimListAddAndOutputConstReal(&
               TRIM(Circuits(p) % names(i))&
               //'re dof '//TRIM(dofnumber), ip(Cvar % ValueId + ReIndex(jj)), Level=10)
             CALL SimListAddAndOutputConstReal(&
               TRIM(Circuits(p) % names(i))&
               //'im dof '//TRIM(dofnumber), ip(Cvar % ValueId + ImIndex(jj)), Level=10)
           END DO
         END IF
       ELSE
         CALL SimListAddAndOutputConstReal(&
           TRIM(Circuits(p) % names(i)), ip(Cvar % ValueId), Level=10)
         
         IF (Cvar % pdofs /= 0 ) THEN
           DO jj = 1, Cvar % pdofs
             write (dofnumber, "(I2)") jj
             CALL SimListAddAndOutputConstReal(&
               TRIM(Circuits(p) % names(i))&
               //'dof '//TRIM(dofnumber), ip(Cvar % ValueId + jj), Level=10)
           END DO
         END IF
       END IF

     END DO

     CALL Info('CircuitsOutput', 'Writing Component Variables for &
       Circuit '//TRIM(i2s(p)), Level=3) 
     DO j = 1, SIZE(Circuits(p) % Components)
         Comp => Circuits(p) % Components(j)
         IF (Comp % Resistance < TINY(0._dp) .AND. Comp % Conductance > TINY(0._dp)) &
             Comp % Resistance = 1._dp / Comp % Conductance

         CALL SimListAddAndOutputConstReal('r_component('//&
           TRIM(i2s(Comp % ComponentId))//')', Comp % Resistance, Level=8) 

         Current = 0._dp + im * 0._dp
         Current = ip(Comp % ivar % ValueId) 
         IF ( Circuits(p) % Harmonic ) Current = Current + im * ip(Comp % ivar % ImValueId) 
              
         CompParams => CurrentModel % Components (Comp % ComponentId) % Values
         IF (.NOT. ASSOCIATED(CompParams)) CALL Fatal ('CircuitsOutput', &
           'Component parameters not found!')

         p_dc_component = ABS(Current)**2._dp * Comp % Resistance
         CALL SimListAddAndOutputConstReal('p_dc_component('//TRIM(i2s(Comp % ComponentId))//')',&
           p_dc_component, Level=8) 

         CompRealPower = GetConstReal( Model % Simulation, 'res: Power re & 
                 in Component '//TRIM(i2s(Comp % ComponentId)), Found)
         IF (Found .AND. ABS(Current) > TINY(CompRealPower)) THEN
           CALL SimListAddAndOutputConstReal('p_ac_component('//&
             TRIM(i2s(Comp % ComponentId))//')', CompRealPower, Level=8)
           CALL SimListAddAndOutputConstReal('r_ac_component('//&
             TRIM(i2s(Comp % ComponentId))//')', CompRealPower/ABS(Current)**2._dp, Level=8)
           CALL SimListAddAndOutputConstReal('AC to DC of component '&
             //TRIM(i2s(Comp % ComponentId)), CompRealPower/p_dc_component, Level=8)
         END IF
          
       END DO  
   END DO

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

  WRITE(VarVal,'(ES15.4)') VariableValue
  CALL Info('CircuitsOutput', TRIM(VariableName)//' '//&
    TRIM(VarVal), Level=LevelVal)

  CALL ListAddConstReal(GetSimulation(),'res: '//TRIM(VariableName), VariableValue)
!-------------------------------------------------------------------
  END SUBROUTINE SimListAddAndOutputConstReal
!-------------------------------------------------------------------


END SUBROUTINE CircuitsOutput
