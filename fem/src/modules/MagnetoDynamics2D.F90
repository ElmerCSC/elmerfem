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
! *  Module for solving magnetic vector potential in Cartesian and
! *  cylindrically symmetric 2D cases. In both cases the vector potential
! *  is reduced to a single component. 
! *
! *  Authors: Juha Ruokolainen, Mika Malinen, Peter Råback
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 30.11.2012
! *
! *****************************************************************************/

!> \ingroup Solvers
!> \{

!------------------------------------------------------------------------------
SUBROUTINE MagnetoDynamics2D_Init( Model,Solver,dt,Transient ) ! {{{
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: Transient           !< Steady state or transient simulation
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: HandleAsm, Found, ElectroDynamics
  CHARACTER(*), PARAMETER :: Caller = 'MagnetoDynamics2D_Init'
 
  Params => GetSolverParams()
  CALL ListAddInteger( Params, 'Variable Dofs',1 )
  CALL ListAddNewString( Params,'Variable','Potential')
  CALL ListAddNewLogical( Params,'Apply Mortar BCs',.TRUE.)
  CALL ListAddNewLogical( Params,'Use Global Mass Matrix',.TRUE.)

  ElectroDynamics = GetLogical (Params, 'Electrodynamics Model', Found)
  IF(ElectroDynamics) THEN
    CALL ListAddInteger( Params, 'Time Derivative Order', 2)
  END IF

  HandleAsm = ListGetLogical( Params,'Handle Assembly',Found )
  
  IF( HandleAsm ) THEN
    IF( CurrentCoordinateSystem() == AxisSymmetric .OR. &
        CurrentCoordinateSystem() == CylindricSymmetric ) THEN 
      CALL Warn(Caller,'Handle assembly not yet available in axisymmetric case!')
      HandleAsm = .FALSE.
    END IF
    IF( ListGetLogicalAnyMaterial(Model, 'Zirka material') ) THEN
      CALL Warn(Caller,'Handle assembly not yet available for Zirka material!')
      HandleAsm = .FALSE.
    END IF
    IF( ListCheckPresentAnyBodyForce(Model, 'Lorentz velocity') ) THEN
      CALL Warn(Caller,'Handle assembly not yet available for "Lorentz velocity"!')
      HandleAsm = .FALSE.
    END IF
    IF( ListCheckPresentAnyBodyForce(Model, 'Angular Velocity') ) THEN
      CALL Warn(Caller,'Handle assembly not yet available for "Angular Velocity"!')
      HandleAsm = .FALSE.
    END IF
    IF(.NOT. HandleAsm ) THEN
      CALL Info(Caller,'Reverting to old bulk assembly routine!')
      CALL ListAddLogical( Params, 'Handle Assembly',.FALSE. )
    END IF
  END IF

  ! Historically a real array could be used for H-B Curve.
  ! This dirty piece of code makes things backward compatible.
  BLOCK
    INTEGER :: i
    LOGICAL :: Cubic
    TYPE(ValueList_t), POINTER :: Material
    DO i=1,Model % NumberOfMaterials
      Material => Model % Materials(i) % Values
      IF( ListCheckPresent( Material, 'H-B Curve') ) THEN
        Cubic = GetLogical( Material, 'Cubic spline for H-B curve',Found)
        CALL ListRealArrayToDepReal(Material,'H-B Curve','dummy',&
            CubicTable=Cubic) !,Monotone=.TRUE.)         
      END IF
    END DO
  END BLOCK
  
!------------------------------------------------------------------------------
END SUBROUTINE MagnetoDynamics2D_Init ! }}}
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solve the magnetic vector potential expressed in terms of a single component.
!> The solver may take into account rotating boundary conditions.
!> Also optionally compute moments and inertia. 
!------------------------------------------------------------------------------
SUBROUTINE MagnetoDynamics2D( Model,Solver,dt,Transient ) ! {{{
!------------------------------------------------------------------------------
  USE DefUtils
  USE CircuitUtils
  USE ZirkaUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: Transient           !< Steady state or transient simulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  LOGICAL :: Found
  TYPE(Element_t), POINTER :: Element
  REAL(KIND=dp) :: Norm
  INTEGER :: i,j,k,n, nb, nd, t, Active, NonlinIter, iter, tind
  TYPE(ValueList_t), POINTER :: BC
  TYPE(Mesh_t),   POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: SolverParams
  
  LOGICAL :: NewtonRaphson = .FALSE., CSymmetry, SkipDegenerate, &
      HandleAsm, MassAsm, ConstantMassInUse = .FALSE.
  LOGICAL :: SliceAverage
  TYPE(Variable_t), POINTER :: CoordVar

  REAL(KIND=dp), ALLOCATABLE, SAVE :: MassValues(:)

  TYPE(TabulatedBasisAtIp_t), POINTER, SAVE :: BasisFunctionsAtIp(:)=>NULL()
  LOGICAL, SAVE :: BasisFunctionsInUse = .FALSE.
  LOGICAL :: UseTorqueTol, UseNewtonRelax, ElectroDynamics
  REAL(KIND=dp) :: TorqueTol, TorqueErr, PrevTorque, Torque, NewtonRelax

  CHARACTER(*), PARAMETER :: Caller = 'MagnetoDynamics2D'
  
!------------------------------------------------------------------------------

  CALL Info( Caller,'------------------------------------------------', Level=4 )
  CALL Info( Caller, 'Solving equation for magnetic vector potential', Level=4 )
  CALL Info( Caller,'------------------------------------------------', Level=4 )

  CALL DefaultStart()
  
  ! Allocate some permanent storage, this is done first time only:
  ! --------------------------------------------------------------
  NULLIFY(BC)
  Mesh => GetMesh()
  SolverParams => GetSolverParams()

  ElectroDynamics = GetLogical (SolverParams, 'Electrodynamics Model', Found)

  IF( ListGetLogical( SolverParams,'Store Basis Functions',Found ) ) THEN
    CALL TabulateBasisFunctions()
  END IF
  
  CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )

  HandleAsm = ListGetLogical( SolverParams,'Handle Assembly',Found )
  IF( HandleAsm ) THEN
    CALL Info(Caller,'Performing handle version of bulk element assembly',Level=7)
  ELSE
    CALL Info(Caller,'Performing legacy version of bulk element assembly',Level=7)      
  END IF

  MassAsm = Transient
  IF( ConstantMassInUse ) MassAsm = .FALSE.
  
  NewtonRaphson = GetLogical(SolverParams, 'Newton-Raphson Iteration', Found)
  IF(GetCoupledIter()>1) NewtonRaphson = .TRUE.
  NewtonRelax = GetCReal(SolverParams,'Nonlinear System Newton Relaxation',UseNewtonRelax)
    
  TorqueTol = GetCReal(SolverParams,'Nonlinear System Torque Tolerance',UseTorqueTol)
  IF(UseTorqueTol) CALL Info(Caller,'Using additional nonlinear tolerance for torque',Level=10)
  Torque = 0.0_dp
  
  NonlinIter = GetInteger(SolverParams,'Nonlinear System Max Iterations',Found)
  IF(.NOT.Found) NonlinIter = 1

  SkipDegenerate = GetLogical(SolverParams, 'Skip Degenerate Elements',Found ) 
  
  CALL Info(Caller,'Initializing Zirka hysteresis models', Level=10)
  CALL InitHysteresis(Model, Solver)

  DO iter = 1,NonlinIter
    CALL Info(Caller,'Performing nonlinear iteration: '//I2S(iter),Level=12)

    IF(Iter > 1) NewtonRaphson=.TRUE.
    ! System assembly:
    ! ----------------

    Active = GetNOFActive()
    CALL DefaultInitialize()
    IF( ConstantMassInUse ) THEN
      Solver % Matrix % MassValues = MassValues
    END IF

    tind = 0
    !$omp parallel do private(Element,n,nd,nb,t)   
    DO t=1,active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes(Element)
      nd = GetElementNOFDOFs(Element)
      nb = GetElementNOFBDOFs(Element)
      IF( SkipDegenerate .AND. DegenerateElement( Element ) ) THEN
        CALL Info(Caller,'Skipping degenerate element:'//I2S(t),Level=12)
        CYCLE
      END IF
      IF( HandleAsm ) THEN
        CALL LocalMatrixHandles(  Element, n, nd+nb, nb )
      ELSE
        CALL LocalMatrix(Element, n, nd)
      END IF
    END DO
    !$omp end parallel do  
      
    CALL DefaultFinishBulkAssembly()
    
    Active = GetNOFBoundaryElements()
!$omp parallel do private(Element, n, nd, BC,Found, t)
    DO t=1,active
      Element => GetBoundaryElement(t)
      BC=>GetBC( Element )
      IF(.NOT.ASSOCIATED(BC)) CYCLE

      IF(GetLogical(BC,'Infinity BC',Found)) THEN
         n  = GetElementNOFNodes( Element )
         nd = GetElementNOFDOFs( Element )
         CALL LocalMatrixInfinityBC(Element, n, nd)
      ELSE IF(GetLogical(BC,'Air Gap',Found)) THEN
         n  = GetElementNOFNodes( Element )
         nd = GetElementNOFDOFs( Element )
         CALL LocalMatrixAirGapBC(Element, BC, n, nd)
      END IF
    END DO
!$omp end parallel do

    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
    
    CALL SetMagneticFluxDensityBC()
    CALL DefaultDirichletBCs()

    IF( ListGetLogical( SolverParams,'Constant Mass Matrix',Found ) ) THEN
      IF( .NOT. ConstantMassInUse ) THEN
        ALLOCATE( MassValues( SIZE( Solver % Matrix % MassValues ) ) )
        MassValues = Solver % Matrix % MassValues 
        ConstantMassInUse = .TRUE.
        MassAsm = .FALSE.
      END IF
    END IF
    
    Norm = DefaultSolve()

    IF( UseTorqueTol ) THEN
      PrevTorque = Torque 
      CALL CalculateLumpedTransient(Torque)
      IF( iter < 2 ) THEN
        ! Cannot have torque tolerance with just one iteration
        CYCLE
      ELSE
        TorqueErr = 2 * ABS(PrevTorque-Torque) / (ABS(PrevTorque)+ABS(Torque))
        IF( TorqueErr > TorqueTol ) THEN
          WRITE(Message,'(A,ES12.3)') 'Torque error at iteration '//I2S(iter)//':',TorqueErr
          CALL Info(Caller,Message,Level=6)
        END IF
      END IF
    END IF

    
    CALL Info(Caller,'Convergence status: '//I2S(Solver % Variable % NonlinConverged),Level=12)
    IF( DefaultConverged() ) THEN
      CALL Info(Caller,'System has converged to tolerances after '//I2S(iter)//' iterations!',Level=12)
      IF( UseTorqueTol ) THEN
        IF( TorqueErr > TorqueTol ) THEN
          CALL Info(Caller,'Nonlinear system tolerance ok after '&
              //I2S(iter)//' but torque still wobbly!',Level=7)
          CYCLE
        END IF
      END IF
      EXIT
    END IF
  END DO
  
  ! For cylindrical symmetry the model lumping has not been implemented
  IF( .NOT. CSymmetry ) THEN
    IF(.NOT. UseTorqueTol ) THEN
      CALL CalculateLumpedTransient()
    END IF
  END IF

  CALL DriveHysteresis(model, solver)

  ! This updates coordinates when using ElmerPost for visualization
  CoordVar => VariableGet(Mesh % Variables,'Coordinates')
  IF(ASSOCIATED(CoordVar)) THEN
    DO i=1,Mesh % NumberOfNodes
      j = 3*(CoordVar % Perm(i)-1)
      CoordVar % Values(j+1) = Mesh % Nodes % x(i)
      CoordVar % Values(j+2) = Mesh % Nodes % y(i)
      CoordVar % Values(j+3) = Mesh % Nodes % z(i)
    END DO
  END IF

  CALL DefaultFinish()

  CALL Info(Caller,'All done',Level=8)
  
CONTAINS


  !> Tabulate basis functions and their weights so that we do not need to compute them in the assembly
  !> process. This assumes that the geometry is not changing. This could be later moved to library
  !> but for now we use a local implementation. 
  !---------------------------------------------------------------------------------------------------
  SUBROUTINE TabulateBasisFunctions()

    INTEGER :: i, t, n, tind
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t), SAVE :: Nodes
    LOGICAL :: stat
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
    REAL(KIND=dp) :: detJ, Weight
    INTEGER :: Phase

    IF( BasisFunctionsInUse ) RETURN
    
    n = Mesh % MaxElementDofs
    ALLOCATE(Basis(n), dBasisdx(n,3))
    
    DO Phase = 0,1
      
      tind = 0
      
      DO i=1,GetNOFActive()
        Element => GetActiveElement(i)
        
        IP = GaussPointsAdapt( Element )      
        CALL GetElementNodes( Nodes, UElement=Element )
        n  = GetElementNOFNodes(Element)     
        
        DO t=1,IP % n
          tind = tind + 1
          
          stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )
          Weight = IP % s(t) * DetJ

          ! Only populate the table the 2nd round
          IF( Phase == 1 ) THEN
            ALLOCATE(BasisFunctionsAtIp(tind) % Basis(n))
            BasisFunctionsAtIp(tind) % Basis(1:n) = Basis(1:n)
            ALLOCATE(BasisFunctionsAtIp(tind) % dBasisdx(n,3))
            BasisFunctionsAtIp(tind) % dBasisdx(1:n,1:3) = dBasisdx(1:n,1:3)
            BasisFunctionsAtIp(tind) % Weight = Weight          
          END IF
        END DO
      END DO

      ! Allocate for the basis functions when the size has been computed
      IF(Phase == 0) THEN
        ALLOCATE( BasisFunctionsAtIp(tind) )
      END IF
    END DO
      
    DEALLOCATE(Basis, dBasisdx)    

    BasisFunctionsInUse = .TRUE.
        
    CALL Info(Caller,'Number of tabulated basis functions:'//I2S(tind),Level=5)
    
  END SUBROUTINE TabulateBasisFunctions
  

!------------------------------------------------------------------------------
! This is monolithic lumping routine that has been optimized for speed.
! This way we need to evaluate the Basis functions only once for each element.
! It is assumed that inertial moment is computed the 1st time only and it
! stays constant.
!------------------------------------------------------------------------------
 SUBROUTINE CalculateLumpedTransient(Torque)
!------------------------------------------------------------------------------
   REAL(KIND=dp), OPTIONAL :: Torque

   REAL(KIND=dp) :: torq,TorqArea,IMoment,IA, &
       rinner,router,rmean,rdiff,ctorq,detJ,Weight,&
       Bp,Br,Bx,By,x,y,r,rho
   REAL(KIND=dp), ALLOCATABLE :: a(:),u(:),POT(:),dPOT(:), &
       pPot(:),Density(:)
   REAL(KIND=dp), POINTER :: Basis(:), dBasisdx(:,:)
   LOGICAL, ALLOCATABLE :: TorqueElem(:)
   INTEGER :: i,bfid,n,nd,nbf,PrevComm,NoSlices,NoTimes
   LOGICAL :: Found, Stat
   TYPE(ValueList_t),POINTER::Params
   TYPE(GaussIntegrationPoints_t) :: IP
   TYPE(Nodes_t) :: Nodes
   LOGICAL :: CalcTorque, CalcPot, CalcInert
   LOGICAL :: ThisTorque, ThisPot, ThisInert, Parallel, HaveRange
   LOGICAL :: Visited = .FALSE.
   
   SAVE Visited, Nodes, Basis, dBasisdx, a, u, POT, dPOT, pPot, Density, Ctorq, TorqueElem
   
!------------------------------------------------------------------------------

   CALL Info(Caller,'Calculating lumped parameters',Level=8)

   NoSlices = MAX(1,ListGetInteger( Model % Simulation,'Number Of Slices', SliceAverage ) )
   NoTimes = ListGetInteger( Model % Simulation,'Number Of Times', Found )
   IF( NoTimes > 1 ) THEN
     PrevComm = ParEnv % ActiveComm
     ParEnv % ActiveComm = ParallelSlicesComm() 
   END IF
      
   ! Define whether we have something to compute
   rinner = ListGetCRealAnyBody( Model,'r inner',CalcTorque )
   IF( CalcTorque ) THEN
     router = ListGetCRealAnyBody( Model,'r outer')
     rmean = (rinner+router)/2
     rdiff = (router-rinner)
     HaveRange = .TRUE.     
   ELSE
     rmean = ListGetConstReal( CurrentModel % Simulation,'Rotor Radius',CalcTorque)
     rdiff = ListGetConstReal( CurrentModel % Simulation,'Rotor Air Gap Width',Found)
     IF(.NOT. Found ) rdiff = 1.0e-3 * rmean
     HaveRange = .FALSE.
   END IF
   
   CalcPot = ListGetLogicalAnyBodyForce( Model,'Calculate Potential' )
   CalcInert = CalcTorque .AND. .NOT. Visited 

   IF( PRESENT(Torque) .AND. .NOT. CalcTorque ) THEN
     CALL Fatal(Caller,'Torque tolerance requested, but torque not computed!')
   END IF
        
   IF(.NOT. (CalcTorque .OR. CalcPot .OR. CalcInert) ) RETURN

   Parallel = ( ParEnv % PEs > 1 )
   
   nbf = Model % NumberOfBodyForces
   IF(.NOT. Visited ) THEN
     n = Model % Mesh % MaxElementDofs
     ALLOCATE( a(nbf), u(nbf), POT(n), dPOT(n), pPot(n), Density(n), &
         Basis(n), dBasisdx(n,3) )
   END IF

   IF( CalcTorque ) THEN
     torq = 0._dp
     TorqArea = 0._dp
   END IF
   IF( CalcInert ) THEN
     IMoment = 0._dp
     IA = 0.0_dp
   END IF
   IF( CalcPot ) THEN   
     U=0._dp
     a=0._dp
   END IF

   tind = 0

   IF(.NOT. Visited .AND. CalcTorque ) THEN
     ALLOCATE( TorqueElem( GetNOFActive() ) )
     TorqueElem = .FALSE.
     
     DO i=1,GetNOFActive()
       Element => GetActiveElement(i)
     
       ThisTorque = .FALSE.
       
       n  = GetElementNOFNodes(Element)     
       CALL GetElementNodes( Nodes, Element )

       IF( HaveRange ) THEN       
         ! We are given range in classical Arkkio style. 
         ! Check how the center lies with respect to the range.
         x = SUM(Nodes % x(1:n))/n
         y = SUM(Nodes % y(1:n))/n
         r = SQRT(x**2+y**2)
         IF (r >= rinner .AND. r <= router) THEN
           TorqueElem(i) = .TRUE.
         END IF
       ELSE
         ! We are not given a range. Just take any element
         ! which has even one node at the given radius. 
         DO j=1,n
           x = Nodes % x(j)
           y = Nodes % y(j)
           r = SQRT(x**2+y**2)
           IF( ABS(r-rmean) < rdiff / 2 ) THEN
             TorqueElem(i) = .TRUE.
             EXIT
           END IF
         END DO
       END IF
     END DO
              
     i = COUNT( TorqueElem )
     CALL Info(Caller,'Number of elements to compute torque: '//I2S(i))
   END IF

   
   DO i=1,GetNOFActive()
     Element => GetActiveElement(i)
     
     ThisTorque = .FALSE.
     ThisPot = .FALSE.
     ThisInert = .FALSE.
     
     IF( CalcPot ) THEN
       Params => GetBodyForce(Element)
       IF(ASSOCIATED(Params)) THEN         
         ThisPot = GetLogical(Params,'Calculate Potential',Found)
         IF( ThisPot ) THEN
           bfid = GetBodyForceId(Element)           
           CALL GetLocalSolution(POT, UElement=Element)
           CALL GetLocalSolution(pPOT,tstep=-1,UElement=Element)
           IF(Solver % Order<2.OR.GetTimeStep()<=2) THEN 
             dPot = (POT - pPOT)/dt
           ELSE
             dPot = 1.5_dp*POT - 2*pPOT
             CALL GetLocalSolution(pPOT,tstep=-2,UElement=Element)
             dPot = (dPOT + 0.5_dp*pPOT)/dt
           END IF
         END IF
       END IF
     END IF
     
     IF( CalcInert ) THEN
       Params=>GetBodyParams(Element)
       IF(ASSOCIATED(Params)) THEN
         ThisInert = GetLogical(Params,'Calculate Inertial Moment',Found)
       END IF
       Density(1:n) = GetReal(GetMaterial(),'Density',Found,Element)
     END IF
     
     IF( CalcTorque ) THEN
       ThisTorque = TorqueElem(i)
       IF(ThisTorque .AND. .NOT. ThisPot ) THEN
         CALL GetLocalSolution(POT, UElement=Element)
       END IF
     END IF

     ! Only treat the element if we have something to compute
     IF( .NOT. (ThisPot .OR. ThisInert .OR. ThisTorque ) ) THEN
       ! If nothing to compute still update the counter for IP points
       IF( BasisFunctionsInUse ) THEN
         IP = GaussPoints(Element)
         tind = tind + IP % n
       END IF
       CYCLE
     END IF
       
     nd = GetElementNOFDOFs(Element)
     n  = GetElementNOFNodes(Element)     
     CALL GetElementNodes( Nodes, Element )
     
     ! Numerical integration:
     !-----------------------
     IP = GaussPoints(Element)
     
     DO t=1,IP % n
       
       ! Basis function values & derivatives at the integration point:
       !--------------------------------------------------------------
       IF( BasisFunctionsInUse ) THEN      
         tind = tind + 1
         Basis => BasisFunctionsAtIp(tind) % Basis
         dBasisdx => BasisFunctionsAtIp(tind) % dBasisdx        
         Weight = BasisFunctionsAtIp(tind) % Weight 
       ELSE IF( ThisTorque ) THEN
         ! Only torque needs the derivatives of basis function
         stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
             IP % W(t), detJ, Basis, dBasisdx )
         weight = IP % s(t) * detJ
       ELSE
         stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
             IP % W(t), detJ, Basis )
         weight = IP % s(t) * detJ
       END IF

       ! Coordinates of the integration point
       x = SUM(Nodes % x(1:n)*Basis(1:n))
       y = SUM(Nodes % y(1:n)*Basis(1:n))
       r = SQRT(x**2+y**2)
       
       IF(ThisPot ) THEN
         A(bfid) = A(bfid) + Weight
         U(bfid) = U(bfid) + Weight * SUM(dPot(1:nd)*Basis(1:nd))
       END IF
       
       IF( ThisTorque ) THEN                      
         Bx =  SUM(POT(1:nd)*dBasisdx(1:nd,2))
         By = -SUM(POT(1:nd)*dBasisdx(1:nd,1))
         Br =  x/r*Bx + y/r*By
         Bp = -y/r*Bx + x/r*By
         Torq = Torq + Weight * r*Br*Bp / (PI*4.0d-7*rdiff)
         TorqArea = TorqArea + Weight
       END IF
       
       IF( ThisInert ) THEN
         IF( r < rmean ) THEN
           rho = SUM( density(1:n) * Basis(1:n) ) 
           IF( rho > EPSILON( rho ) ) THEN
             IA = IA + Weight
             U = U + Weight * r * rho
           END IF
         END IF
       END IF
     END DO
   END DO
     

   ! Finally perform parallel reduction if needed, and
   ! store the results for saving by SaveScalars.
   !-------------------------------------------------------------------------   
   IF( CalcPot ) THEN
     IF( ParEnv % PEs > 1 ) THEN
       DO i=1,nbf
         a(i) = ParallelReduction(a(i)) / NoSlices
         u(i) = ParallelReduction(u(i)) / NoSlices
       END DO
     END IF     
     DO i=1,nbf
       IF(a(i)>0) THEN
         CALL ListAddConstReal(Model % Simulation,'res: Potential / bodyforce ' &
             //i2s(i),u(i)/a(i))
         CALL ListAddConstReal(Model % Simulation,'res: area / bodyforce ' &
             //i2s(i),a(i))
       END IF
     END DO
   END IF
   
   IF( CalcTorque ) THEN   
     ! Arkkios formula assumes that rinner and router are nicely aligned with elements.
     ! This may not the case, so the 1st time we make a geomeric correction. 
     IF(.NOT. Visited ) THEN
       WRITE(Message,'(A,ES15.4)') 'Air gap initial torque:', Torq
       CALL Info(Caller,Message,Level=6)

       TorqArea = ParallelReduction(TorqArea) / NoSlices       
       IF (TorqArea > EPSILON(TorqArea) ) THEN
         Ctorq = 2 * PI * rmean * rdiff / TorqArea

         WRITE(Message,'(A,F8.4)') 'Air gap correction initial:', cTorq
         CALL Info(Caller,Message,Level=4)
         
         ! The correction factor also corrects for the number of periods.
         ! We don't want that - so let us take back that and the torque
         ! can be compared to inertial moment of the sector still. 
         i = ListGetInteger( CurrentModel % Simulation,'Rotor Periods',Found )         
         IF( Parallel ) i = ParallelReduction( i, 2 ) 
         IF( i > 1 ) THEN
           WRITE(Message,'(A,I0)') 'Air gap correction rotor periods: ',i
           CALL Info(Caller,Message,Level=4)
           Ctorq = Ctorq / i 
         END IF         
       ELSE
         Ctorq = 1.0_dp
       END IF
       
       WRITE(Message,'(A,F8.4)') 'Air gap correction:', cTorq
       CALL Info(Caller,Message,Level=4)
       !CALL ListAddConstReal(Model % Simulation,'res: air gap correction', cTorq)
     END IF
       
     Torq = Ctorq * Torq

     IF( SliceAverage ) THEN
       ! Save slice torque even for one slice since then the output for scalars is the same
       ! for any number of slices.       
       WRITE(Message,'(A,ES15.4)') 'Air gap torque for slice'//I2S(ParEnv % MyPe)//':', Torq
       CALL Info(Caller,Message,Level=5)
       CALL ListAddConstReal(Model % Simulation,'res: air gap torque for slice', Torq)
     END IF
       
     ! But the averaging makes sense only for more than one slice
     Torq = ParallelReduction(Torq) / NoSlices
     WRITE(Message,'(A,ES15.4)') 'Air gap torque:', Torq
     CALL Info(Caller,Message,Level=5)
     CALL ListAddConstReal(Model % Simulation,'res: air gap torque', Torq)

     IF(PRESENT(Torque)) Torque = Torq
   END IF

   IF( CalcInert ) THEN
     IF( Parallel ) THEN
       IMoment = ParallelReduction(IMoment) / NoSlices
       IA = ParallelReduction(IA) / NoSlices
     END IF

     WRITE(Message,'(A,ES15.4)') 'Inertial volume:', IA
     CALL Info(Caller,Message,Level=7)

     WRITE(Message,'(A,ES15.4)') 'Inertial moment:', Imoment
     CALL Info(Caller,Message,Level=7)
       
     CALL ListAddConstReal(Model % Simulation,'res: inertial volume', IA)
     CALL ListAddConstReal(Model % Simulation,'res: inertial moment', IMoment)
   END IF
     
   Visited = .TRUE.
   
   ! Revert the communicatior back to original
   IF( NoTimes > 1 ) ParEnv % ActiveComm = PrevComm
   
!------------------------------------------------------------------------------
 END SUBROUTINE CalculateLumpedTransient
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Old style local matrix. 
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(Element, n, nd)
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: Material, BodyForce
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(ValueList_t), POINTER :: CompParams

    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    REAL(KIND=dp) :: MASS(nd,nd), DAMP(nd,nd), STIFF(nd,nd), FORCE(nd), &
        LOAD(nd),R(2,2,n),LondonLambda(nd), C(n), mu,muder,Babs,POT(nd), &
        JAC(nd,nd),Agrad(3),C_ip,M(2,n),M_ip(2),x, y,&
        Lorentz_velo(3,nd), Velo(3), omega_velo
    REAL(KIND=dp) :: LondonLambda_ip, P_ip, Permittivity(nd)
    REAL(KIND=dp) :: Bt(nd,2), Ht(nd,2)
    REAL(KIND=dp) :: nu_tensor(2,2)
    REAL(KIND=dp) :: B_ip(2), Alocal, H_ip(2)

    INTEGER :: i,p,q,t

    LOGICAL :: HBcurve, WithVelocity, WithAngularVelocity, Found, Stat
    LOGICAL :: CoilBody, StrandedCoil    

    CHARACTER(LEN=MAX_NAME_LEN) :: CoilType

    ! Zirka related
    LOGICAL :: Zirka
    LOGICAL :: LondonEquations = .TRUE.
    TYPE(Variable_t), POINTER :: hystvar
    TYPE(GlobalHysteresisModel_t), pointer :: zirkamodel

!$omp threadprivate(Nodes)
    
!------------------------------------------------------------------------------

    IF( UseLocalMatrixCopy( Solver, Element % ElementIndex ) ) GOTO 10

    CALL GetElementNodes( Nodes,Element )
    STIFF = 0._dp
    JAC  = 0._dp
    FORCE = 0._dp
    IF(Transient) THEN
      MASS = 0._dp; DAMP=0._dp
    END IF

    Material => GetMaterial(Element)

    HBCurve = ListCheckPresent(Material,'H-B Curve')
    Zirka = ListGetLogical(Material, 'Zirka material', Zirka)

    IF (zirka) THEN
      CALL GetLocalSolution(POT,UElement=Element,USolver=Solver)
      zirkamodel => GetZirkaPointer(Material)
      hystvar => GetZirkaVariable(Material)
    ELSE IF(HBcurve) THEN
      CALL GetLocalSolution(POT,UElement=Element,USolver=Solver)
    ELSE
      CALL GetReluctivity(Material,R,n,Element)
    END IF

    C = GetReal( Material, 'Electric Conductivity', Found, Element)

    M(1,:) = GetReal( Material, 'Magnetization 1', Found, Element)
    M(2,:) = GetReal( Material, 'Magnetization 2', Found, Element)

    Load = 0.0d0

    WithVelocity = .FALSE.
    WithAngularVelocity = .FALSE.
    
    BodyForce => GetBodyForce(Element)
    IF ( ASSOCIATED(BodyForce) ) THEN
      Load(1:n) = GetReal(BodyForce, 'Current Density', Found, Element)
      CALL GetRealVector(BodyForce, Lorentz_velo, 'Lorentz velocity', WithVelocity)
      omega_velo = ListGetCReal(BodyForce, 'Angular velocity', WithAngularVelocity) 
    END IF

    CoilBody = .FALSE.
    StrandedCoil = .FALSE.
    LondonEquations = .TRUE.
    CompParams => GetComponentParams( Element )
    IF (ASSOCIATED(CompParams)) THEN
      CoilType = GetString(CompParams, 'Coil Type', Found)
      IF (Found) THEN
        CoilBody = .TRUE.
        SELECT CASE (CoilType)
        CASE ('stranded')
          StrandedCoil = .TRUE.
        CASE ('massive')
          LondonEquations = ListGetLogical(CompParams, 'London Equations', LondonEquations)
        CASE ('foil winding')
!          CALL GetElementRotM(Element, RotM, n)
        CASE DEFAULT
          CALL Fatal (Caller, 'Non existent Coil Type Chosen!')
        END SELECT
      END IF
    END IF

    IF ( LondonEquations ) THEN
      ! Lambda = m/(n_s * e^2):
      ! -------------------
      LondonLambda(:) = GetReal( Material, 'London Lambda', LondonEquations, Element)
    END IF

    Permittivity(1:n) = GetReal( Material, 'Permittivity', Found )
    
    !Numerical integration:
    !----------------------
    IP = GaussPoints(Element)
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )

      IF ( CSymmetry ) THEN
        x = SUM( Basis(1:n) * Nodes % x(1:n) )
        detJ = detJ * x
      END IF

      ! The source term at the integration point:
      !------------------------------------------
      LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )

      nu_tensor = 0.0_dp

      IF(Zirka .OR. HBCUrve) THEN
        Agrad = 0.0_dp
        Agrad = MATMUL( POT,dBasisdx )
        Alocal = SUM( POT(1:nd) * Basis(1:nd) )
        ! Sign? This convention: \vec A = A u_z
        ! -----
        B_ip(1) = Agrad(2) 
        B_ip(2) = -Agrad(1)
        IF( CSymmetry ) THEN
          B_ip = -B_ip
          B_ip(2) = B_ip(2) + Alocal/x
        END IF
      END IF

      IF (HBcurve ) THEN
        Babs = MAX( SQRT(SUM(B_ip**2)), 1.d-8 )

        IF( NewtonRaphson ) THEN
          mu = ListGetFun( Material,'h-b curve',babs,dFdx=muder) / Babs
          muder = (muder-mu)/babs
        ELSE
          mu = ListGetFun( Material,'h-b curve',babs) / Babs
        END IF

        nu_tensor(1,1) = mu ! Mu is really nu!!! too lazy to correct now...
        nu_tensor(2,2) = mu
      ELSE IF(Zirka) THEN
        CALL GetZirkaHBAtIP(t, solver, element, hystvar, zirkamodel, B_ip, H_ip, nu_tensor)
      ELSE
        muder=0._dp
        DO p=1,2
          DO q=1,2
            nu_tensor(p,q) = SUM(Basis(1:n) * R(p,q,1:n))
          END DO
        END DO
      END IF

      C_ip = SUM( Basis(1:n) * C(1:n) )
      M_ip = MATMUL( M,Basis(1:n) )
      P_ip = SUM( Basis(1:n) * Permittivity(1:n) )

      ! Finally, the elemental matrix & vector:
      !----------------------------------------
      IF (Transient .AND. C_ip/=0._dp .AND. .NOT. StrandedCoil .OR. ElectroDynamics ) THEN
        DO p=1,nd
          DO q=1,nd
            IF(ElectroDynamics) THEN
              DAMP(p,q) = DAMP(p,q) + IP % s(t) * detJ * C_ip * Basis(q)*Basis(p)
              MASS(p,q) = MASS(p,q) + IP % s(t) * detJ * P_ip * Basis(q)*Basis(p)
            ELSE
              MASS(p,q) = MASS(p,q) + IP % s(t) * detJ * C_ip * Basis(q)*Basis(p)
            END IF
          END DO
        END DO
      END IF

      IF ( LondonEquations ) THEN
        LondonLambda_ip = SUM( Basis(1:n) * LondonLambda(1:n) )
        DO p=1,nd
          DO q=1,nd
            ! (LondonLambda a,  a'):
            ! --------------
            STIFF(p,q) = STIFF(p,q) + IP % s(t) * detJ / LondonLambda_ip * Basis(q)*Basis(p)
          END DO
        END DO
      END IF
      
      ! Is the sign correct?
      !---------------------
      Bt(1:nd,1) =  dbasisdx(1:nd,2)
      Bt(1:nd,2) = -dbasisdx(1:nd,1)
      IF ( CSymmetry ) THEN
        Bt(1:nd,1:2) = -Bt(1:nd,1:2)
        Bt(1:nd,2) = Bt(1:nd,2) + Basis(1:nd)/x
      END IF

      DO p = 1,nd
        Ht(p,:) = MATMUL(nu_tensor, Bt(p,:))
      END DO

      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + IP % s(t) * DetJ * &
             MATMUL(Ht(1:nd,:), TRANSPOSE(Bt(1:nd,:)))

      ! Csymmetry is not yet considered in the Newton linearization
      IF (HBcurve .AND. NewtonRaphson) THEN
!        DO p=1,nd
!          DO q=1,nd
!            JAC(p,q) = JAC(p,q) + IP % s(t) * DetJ * &
!              muder/babs*SUM(Agrad*dBasisdx(q,:))*SUM(Agrad*dBasisdx(p,:))
!          END DO
!        END DO
        DO p=1,nd
          DO q=1,nd
            JAC(p,q) = JAC(p,q) + IP % s(t) * DetJ * &
              muder/babs*SUM(B_ip*Bt(q,:))*SUM(B_ip*Bt(p,:))
          END DO
        END DO
      END IF

      IF (WithVelocity .OR. WithAngularVelocity ) THEN
        !
        ! Create an additional Lorentz effect so that the electric field
        ! has an added term v x curl A:
        !
        IF( WithVelocity ) THEN        
          Velo(1:2) = [ SUM(Basis(1:n)*Lorentz_velo(1,1:n)), &
              SUM(Basis(1:n)*Lorentz_velo(2,1:n)) ]
        ELSE
          x = SUM( Basis(1:n) * Nodes % x(1:n) )
          y = SUM( Basis(1:n) * Nodes % y(1:n) ) 
          
          ! Simplified omega \times r in 2D
          Velo(1) = -omega_velo * y
          Velo(2) = omega_velo * x
        END IF
          
        IF (CSymmetry) THEN
          DO p=1,nd
            STIFF(p,1:nd) = STIFF(p,1:nd) + IP % s(t) * DetJ * C_ip * Basis(p) * ( & 
                -Velo(2) * Bt(1:nd,1) + Velo(1) * Bt(1:nd,2) )
          END DO
        ELSE
          DO p=1,nd
            STIFF(p,1:nd) = STIFF(p,1:nd) + IP % s(t) * DetJ * C_ip * Basis(p) * ( & 
                Velo(2) * Bt(1:nd,1) - Velo(1) * Bt(1:nd,2) )
          END DO
        END IF
      END IF

      IF ( CSymmetry ) THEN
        FORCE(1:nd) = FORCE(1:nd) + IP % s(t) * DetJ * (LoadAtip * Basis(1:nd) - &
            M_ip(1)*dBasisdx(1:nd,2)+M_ip(2)*(dBasisdx(1:nd,1) + Basis(1:nd)/x))
      ELSE
        FORCE(1:nd) = FORCE(1:nd) + IP % s(t) * DetJ * (LoadAtip * Basis(1:nd) + &
            M_ip(1)*dBasisdx(1:nd,2)-M_ip(2)*dBasisdx(1:nd,1))
      END IF
      IF(zirka) then
        FORCE(1:nd) = FORCE(1:nd) - (H_ip(1)*Bt(1:nd,1) + H_ip(2)*Bt(1:nd,2)) * IP % s(t) * detJ
      END IF
    END DO

    IF (HBcurve .AND. NewtonRaphson) THEN
      IF(UseNewtonRelax) JAC = JAC * NewtonRelax
      STIFF = STIFF + JAC
      FORCE = FORCE + MATMUL(JAC,POT)
    END IF

    IF(Zirka) THEN
      FORCE = FORCE + MATMUL(STIFF, POT)
    END IF

    IF(Transient) THEN
      IF(ElectroDynamics) THEN
        CALL Default2ndOrderTime( MASS, DAMP, STIFF, FORCE,UElement=Element, USolver=Solver )
      ELSE
        CALL Default1stOrderTime( MASS, STIFF, FORCE,UElement=Element, USolver=Solver )
      END IF
    END IF
10  CALL DefaultUpdateEquations( STIFF, FORCE,UElement=Element, USolver=Solver)

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Assembly using handles. A little faster even for linear triangles, maybe 20%. 
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixHandles( Element, n, nd, nb )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, nd, nb
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER, SAVE :: Basis(:), dBasisdx(:,:)
    REAL(KIND=dp), ALLOCATABLE, SAVE :: MASS(:,:), DAMP(:,:), STIFF(:,:), FORCE(:), POT(:)    
    REAL(KIND=dp) :: Nu0, Nu, weight, SourceAtIp, CondAtIp, DetJ, Mu, MuDer, Babs
    LOGICAL :: Stat,Found, HBCurve
    INTEGER :: t,p,q,m,allocstat
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(ValueList_t), POINTER :: Material, PrevMaterial => NULL()
    REAL(KIND=dp) :: B_ip(2), Ht(nd,2), Bt(nd,2), Agrad(2), JAC(nd,nd), Alocal, &
            Permittivity(nd), P_ip
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilType
    LOGICAL :: StrandedCoil
    TYPE(ValueHandle_t), SAVE :: SourceCoeff_h, CondCoeff_h, PermCoeff_h, &
        RelPermCoeff_h, RelucCoeff_h, Mag1Coeff_h, Mag2Coeff_h, CoilType_h
    INTEGER :: PrevElemInd = HUGE(PrevElemInd)
    
    SAVE HBCurve, Nu0, PrevMaterial, PrevElemInd
    
    !$omp threadprivate(Basis, dBasisdx, MASS, DAMP, STIFF, FORCE, POT, &
    !$omp               Nodes, Nu0, HBCurve, PrevMaterial, &
    !$omp               SourceCoeff_h, CondCoeff_h, PermCoeff_h, RelPermCoeff_h, &
    !$omp               RelucCoeff_h, Mag1Coeff_h, Mag2Coeff_h, CoilType_h, PrevElemInd )
    
!------------------------------------------------------------------------------

    ! The elements should be in growing order. Hence we initialize if we start the list.
    IF( Element % ElementIndex < PrevElemInd ) THEN
      CALL ListInitElementKeyword( SourceCoeff_h,'Body Force','Current Density')
      CALL ListInitElementKeyword( CondCoeff_h,'Material','Electric Conductivity')
      CALL ListInitElementKeyword( PermCoeff_h,'Material','Permeability')
      CALL ListInitElementKeyword( RelPermCoeff_h,'Material','Relative Permeability')
      CALL ListInitElementKeyword( RelucCoeff_h,'Material','Reluctivity')
      CALL ListInitElementKeyword( Mag1Coeff_h,'Material','Magnetization 1')
      CALL ListInitElementKeyword( Mag2Coeff_h,'Material','Magnetization 2')
      CALL ListInitElementKeyword( CoilType_h,'Component','Coil Type')
      Found = .FALSE.
      IF( ASSOCIATED( Model % Constants ) ) THEN
        Nu0 = ListGetCReal( Model % Constants,'Permeability of Vacuum',Found)
      END IF
      IF( .NOT. Found ) Nu0 = PI * 4.0d-7
    END IF
    PrevElemInd = Element % ElementIndex
    
    ! Allocate storage if needed
    IF (.NOT. ALLOCATED(MASS)) THEN
      m = Mesh % MaxElementDofs
      ALLOCATE(MASS(m,m), DAMP(m,m), STIFF(m,m),FORCE(m), POT(m), STAT=allocstat)      
      IF (allocstat /= 0) THEN
        CALL Fatal(Caller,'Local storage allocation failed')
      END IF
      IF(.NOT. BasisFunctionsInUse ) THEN
        ALLOCATE(Basis(m), dBasisdx(m,3))
      END IF
    END IF
    
    IF( UseLocalMatrixCopy( Solver, Element % ElementIndex ) ) GOTO 20
    
    Material => GetMaterial(Element)
    IF( .NOT. ASSOCIATED( Material, PrevMaterial ) ) THEN
      PrevMaterial => Material           
      HbCurve = ListCheckPresent(Material,'H-B Curve')
    END IF

    IF(ElectroDynamics) THEN
      Permittivity(1:n) = GetReal( Material, 'Permittivity', Found )
    END IF

    StrandedCoil = .FALSE.
    CoilType = ListGetElementString(CoilType_h, Element, Found ) 
    IF( Found ) THEN
      SELECT CASE (CoilType)
      CASE ('stranded')
        StrandedCoil = .TRUE.
      CASE ('massive')
        CONTINUE
      CASE ('foil winding')
        CONTINUE
      CASE DEFAULT
        CALL Fatal (Caller, 'Non existent Coil Type Chosen!')
      END SELECT
    END IF
        
    ! Initialize
    MASS  = 0.0_dp
    DAMP  = 0.0_dp
    STIFF = 0.0_dp
    FORCE = 0.0_dp
    
    ! Integration rule
    IP = GaussPointsAdapt( Element )
      
    CALL GetElementNodes( Nodes, UElement=Element )

    IF(HBcurve) THEN
      CALL GetLocalSolution(POT,UElement=Element,USolver=Solver)
      JAC = 0.0_dp
    END IF
    
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      IF( BasisFunctionsInUse ) THEN      
        tind = tind + 1
        Basis => BasisFunctionsAtIp(tind) % Basis
        dBasisdx => BasisFunctionsAtIp(tind) % dBasisdx        
        Weight = BasisFunctionsAtIp(tind) % Weight 
      ELSE
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis, dBasisdx )
        Weight = IP % s(t) * DetJ
      END IF
        
      ! diffusion term (D*grad(u),grad(v)):
      ! -----------------------------------
      IF( HBCurve ) THEN
        Agrad(1:2) = MATMUL( POT(1:nd),dBasisdx(1:nd,1:2) )
        Alocal = SUM( POT(1:nd) * Basis(1:nd) )

        B_ip(1) = Agrad(2) 
        B_ip(2) = -Agrad(1)         
        Babs = MAX( SQRT(SUM(B_ip**2)), 1.d-8 )

        IF( NewtonRaphson ) THEN
          mu = ListGetFun( Material,'h-b curve',babs,dFdx=muder) / Babs
          muder = (muder-mu)/babs
        ELSE
          mu = ListGetFun( Material,'h-b curve',babs) / Babs
        END IF
      ELSE
        Nu = ListGetElementReal( RelPermCoeff_h, Basis, Element, Found, GaussPoint = t )
        IF( Found ) THEN
          Mu = 1.0_dp / (Nu0 * Nu)
        ELSE
          Nu = ListGetElementReal( PermCoeff_h, Basis, Element, Found, GaussPoint = t )
          IF( Found ) THEN
            Mu = 1.0_dp / Nu
          ELSE
            Mu = ListGetElementReal( RelucCoeff_h, Basis, Element, Found, GaussPoint = t )
          END IF

          IF(.NOT. Found ) THEN
            CALL Fatal(Caller,'Could not define reluctivity in any way in Body: '&
                //I2S(Element % BodyId))
          END IF
        END IF
      END IF

      P_ip = SUM( Basis(1:n) * Permittivity )

      Bt(1:nd,1) =  dbasisdx(1:nd,2)
      Bt(1:nd,2) = -dbasisdx(1:nd,1)

      ! Here isotrophy is assumed!
      Ht(1:nd,:) = mu * Bt(1:nd,:)
           
      IF ( HBCurve .AND. NewtonRaphson) THEN
        DO p=1,nd
          DO q=1,nd
            JAC(p,q) = JAC(p,q) + Weight * &
                muder/babs * SUM(B_ip(:) * Bt(q,:)) * SUM(B_ip(:)*Bt(p,:))
          END DO
        END DO
      END IF
            
      ! diffusive term: STIFF=STIFF+(a*grad(u),grad(v))   
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + Weight * &
          MATMUL(Ht(1:nd,:), TRANSPOSE(Bt(1:nd,:)))

      IF( MassAsm .AND. .NOT. StrandedCoil ) THEN
        CondAtIp = ListGetElementReal( CondCoeff_h, Basis, Element, Found )
        IF( Found ) THEN
          DO p=1,nd
            IF(ElectroDynamics) THEN
              DAMP(p,1:nd) = DAMP(p,1:nd) + Weight * CondAtIp * Basis(1:nd) * Basis(p)
              MASS(p,1:nd) = MASS(p,1:nd) + Weight * P_ip * Basis(1:nd) * Basis(p)
            ELSE
              MASS(p,1:nd) = MASS(p,1:nd) + Weight * CondAtIp * Basis(1:nd) * Basis(p)
            END IF
          END DO
        END IF
      END IF

      ! Current density source 
      SourceAtIP = ListGetElementReal( SourceCoeff_h, Basis, Element, Found ) 
      IF( Found ) THEN
        FORCE(1:nd) = FORCE(1:nd) + Weight * SourceAtIP * Basis(1:nd)
      END IF

      ! Magnetization source, weak form
      SourceAtIP = ListGetElementReal( Mag1Coeff_h, Basis, Element, Found ) 
      IF( Found ) THEN
        FORCE(1:nd) = FORCE(1:nd) + Weight * SourceAtIP * dBasisdx(1:nd,2)        
      END IF
      SourceAtIP = ListGetElementReal( Mag2Coeff_h, Basis, Element, Found ) 
      IF( Found ) THEN
        FORCE(1:nd) = FORCE(1:nd) - Weight * SourceAtIP * dBasisdx(1:nd,1)        
      END IF     
    END DO

    IF (HBcurve .AND. NewtonRaphson) THEN
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + JAC(1:nd,1:nd)
      FORCE(1:nd) = FORCE(1:nd) + MATMUL(JAC(1:nd,1:nd),POT(1:nd))
    END IF
    
    IF( MassAsm ) THEN
      IF(ElectroDynamics) THEN
        CALL DefaultUpdateDamp(DAMP,UElement=Element)
        CALL DefaultUpdateMass(MASS,UElement=Element)
      ELSE
        CALL DefaultUpdateMass(MASS,UElement=Element)
      END IF
    END IF 
    CALL CondensateP( nd-nb, nb, STIFF, FORCE )
    
20  CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element) !, VecAssembly=VecAsm)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixHandles
!------------------------------------------------------------------------------

  
!-------------------------------------------------------------------------------
! Calculates H and dHdB in 2D given B. This should be always inlined in LocalMatrix.
!-------------------------------------------------------------------------------
SUBROUTINE GetZirkaHBAtIP(i_IP, Solver, Element, HystVar, ZirkaModel, B_ip, H_ip, dHdB) ! {{{
!-------------------------------------------------------------------------------
  INTEGER, intent(in) :: i_IP
  TYPE(Solver_t) :: Solver
  TYPE(Element_t) :: Element
  TYPE(Variable_t), POINTER :: HystVar
  TYPE(GlobalHysteresisModel_t), POINTER :: ZirkaModel
  REAL(KIND=dp), INTENT(IN) :: B_ip(2)
  REAL(KIND=dp), INTENT(OUT) :: H_ip(2)
  REAL(KIND=dp), INTENT(INOUT) :: dHdB(2,2)
!-------------------------------------------------------------------------------
  INTEGER :: ipindex, n_dir, k,l
  REAL(KIND=dp) :: dH, B0(3)
!-------------------------------------------------------------------------------
  ipindex = getipindex(i_IP, usolver=solver, element=element, ipvar=hystvar)
  IF (ipindex == 0 ) RETURN

  H_ip = 0.0_dp
  DO n_dir = 1, UBOUND(zirkamodel % curves, 1)
    B0 = zirkamodel % curves(n_dir, ipindex) % B0
    ASSOCIATE(Bdir => SUM(B_ip*B0(1:2)))
      ! H_ip(1:2) = H_ip(1:2) + zirkamodel % curves(n_dir,ipindex) % &
      !     eval(sum(B_ip*B0(1:2)), cached = .true., dhdb=dH) * &
      !     B0(1:2)
      H_ip(1:2) = H_ip(1:2) + zirkamodel % curves(n_dir,ipindex) % &
          eval(Bdir, cached = .TRUE., dhdb=dH) * &
          B0(1:2)
    END ASSOCIATE
    DO k = 1,2
      DO l = 1,2
        dHdB(k,l) = dHdB(k,l) + dH*B0(k)*B0(l)
      END DO
    END DO
  END DO
  
END SUBROUTINE ! }}}
!-------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixInfinityBC(Element, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),DetJ
    LOGICAL :: Stat
    INTEGER :: i,p,q,t
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), R(2,2,n), &
            Inf_ip,Coord(3),Normal(3),mu,u,v
    TYPE(ValueList_t), POINTER :: Material
    TYPE(Element_t), POINTER :: Parent
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !$OMP THREADPRIVATE(Nodes)
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes, Element )
    STIFF = 0._dp
    FORCE = 0._dp

    Parent=>Element % BoundaryInfo % Left
    IF(.NOT.ASSOCIATED(Parent)) THEN
      Parent=>Element % BoundaryInfo % Right
    END IF
    IF(.NOT.ASSOCIATED(Parent)) RETURN

    Material => GetMaterial(Parent)
    CALL GetReluctivity(Material,R,n,Parent)

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                 IP % W(t), detJ, Basis )

      mu = SUM(Basis(1:n)*R(1,1,1:n)) ! We assume isotropic reluctivity here.

      Normal = NormalVector( Element, Nodes, u, v, .TRUE. )
      Coord(1) = SUM(Basis(1:n) * Nodes % x(1:n))
      Coord(2) = SUM(Basis(1:n) * Nodes % y(1:n))
      Coord(3) = SUM(Basis(1:n) * Nodes % z(1:n))

      IF( CSymmetry ) THEN
        detJ = detJ * Coord(1)
      END IF

      Inf_ip = mu * SUM(Coord*Normal)/SUM(Coord*Coord)

      DO p=1,nd
        DO q=1,nd
          STIFF(p,q) = STIFF(p,q) + IP % s(t)*detJ*Inf_ip*Basis(q)*Basis(p)
        END DO
      END DO
    END DO
    CALL DefaultUpdateEquations( STIFF, FORCE, UElement=Element )
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixInfinityBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixAirGapBC(Element, BC, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ
    LOGICAL :: Stat, Found
    INTEGER :: t
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), &
            mu,AirGapLength(nd), AirGapMu(nd), AirGapL
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !$OMP THREADPRIVATE(Nodes)
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes, Element )
    STIFF = 0._dp
    FORCE = 0._dp

    AirGapLength = ListGetConstReal( BC, 'Air Gap Length', UnfoundFatal = .TRUE.)

    AirGapMu = ListGetConstReal( BC, 'Air Gap Relative Permeability', Found)
    IF (.NOT. Found) AirGapMu=1.0_dp

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )

      mu = 4*pi*1d-7*SUM(Basis(1:n)*AirGapMu(1:n))
      AirGapL = SUM(Basis(1:n)*AirGapLength(1:n))

      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + IP % s(t) * DetJ * &
          AirGapL/mu*MATMUL(dBasisdx, TRANSPOSE(dBasisdx))

    END DO
    CALL DefaultUpdateEquations( STIFF, FORCE, UElement=Element )
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixAirGapBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE SetMagneticFluxDensityBC()
!------------------------------------------------------------------------------
! P. Lombard, G. Meunier, "A general purpose method for electric and magnetic 
! combined problems for 2D, axisymmetric and transient systems", IEEE Trans.
! magn. 29(2), p. 1737 - 1740, Mar 1993
! -ettaka- 
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Matrix_t), POINTER :: A
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp), POINTER :: b(:)
    INTEGER :: i, n, j, k
    TYPE(ValueList_t), POINTER :: BC
    LOGICAL :: Found
    REAL(KIND=dp) :: Bx(Solver % Mesh % MaxElementDofs), &
        By(Solver % Mesh % MaxElementDofs)
    REAL(KIND=dp) :: x, y
    INTEGER, POINTER :: Perm(:)

    Perm => Solver % Variable % Perm
    A => Solver % Matrix
    b => A % RHS
    
    DO i=1,GetNofBoundaryElements()
      Element => GetBoundaryElement(i)
      n = GetELementNofNodes()
      BC => GetBC()
      IF ( ASSOCIATED(BC)) THEN
        IF ( ListCheckPrefix( BC, 'Magnetic Flux Density') ) THEN
          Bx = 0._dp
          By = 0._dp

          Bx(1:n) = GetReal(BC, 'Magnetic Flux Density 1', Found)
          By(1:n) = GetReal(BC, 'Magnetic Flux Density 2', Found)

          DO j = 1,n
            k = Element % NodeIndexes(j)
            x = Mesh % Nodes % x(k)
            y = Mesh % Nodes % y(k)
            k = Perm(k)

            CALL UpdateDirichletDof( A, k, y * Bx(j) - x * By(j) )
          END DO 
        END IF  
      END IF  
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE SetMagneticFluxDensityBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE GetReluctivity(Material,Acoef,n,Element)
!------------------------------------------------------------------------------
    USE MGDynMaterialUtils
    TYPE(ValueList_t), POINTER :: Material
    INTEGER :: n
    REAL(KIND=dp) :: Acoef(2,2,n)
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp), SAVE :: Avacuum
    LOGICAL :: Found
    LOGICAL, SAVE :: FirstTime = .TRUE.
    !$OMP THREADPRIVATE(Avacuum, FirstTime)
!------------------------------------------------------------------------------

    IF ( FirstTime ) THEN
      Avacuum = GetConstReal( CurrentModel % Constants, &
              'Permeability of Vacuum', Found )
      IF(.NOT. Found ) Avacuum = PI * 4.0d-7
      FirstTime = .FALSE.
    END IF

    Acoef = GetTensor(Element, n, 2, 'Relative Permeability', 're', Found)

    IF ( Found ) THEN
      Acoef = Avacuum * Acoef
    ELSE
      Acoef = GetTensor(Element, n, 2, 'Permeability', 're', Found)
    END IF
    IF ( Found ) THEN
      Acoef = Get2x2TensorInverse(Acoef, n)
    ELSE
      Acoef = GetTensor(Element, n, 2, 'Reluctivity', 're', Found)
    END IF

    IF( .NOT. Found ) THEN
      CALL Warn('GetReluctivity',&
          'Could not get either > Reluctivity > or > Relative Permeability < !')
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE GetReluctivity
!------------------------------------------------------------------------------

END SUBROUTINE MagnetoDynamics2D
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE MagnetoDynamics2DHarmonic_Init0( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: Transient           !< Steady state or transient simulation
!------------------------------------------------------------------------------
  CALL ListAddNewLogical( Solver % Values, 'Apply Mortar BCs', .TRUE.)
  CALL ListAddNewLogical( Solver % Values, 'Linear System Complex', .TRUE.)
!------------------------------------------------------------------------------
END SUBROUTINE MagnetoDynamics2DHarmonic_Init0
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE MagnetoDynamics2DHarmonic_Init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: Transient           !< Steady state or transient simulation
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params

  Params => GetSolverParams(Solver)
  CALL ListAddInteger( Params, 'Variable Dofs',2 )
  CALL ListAddNewString( Params,'Variable',&
      'Potential[Potential re:1 Potential im:1]')

  CALL ListAddNewLogical( Params,'Apply Mortar BCs',.TRUE.)
  CALL ListAddNewLogical( Params,'Linear System Complex',.TRUE.)

!------------------------------------------------------------------------------
END SUBROUTINE MagnetoDynamics2DHarmonic_Init
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solve the complex magnetic vector potential having a single component.
!> The solver may take into account rotating boundary conditions.
!> Also optionally compute moments and inertia. 
!------------------------------------------------------------------------------
SUBROUTINE MagnetoDynamics2DHarmonic( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  USE CircuitUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: Transient           !< Steady state or transient simulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  LOGICAL :: Found, ElectroDynamics
  TYPE(Element_t),POINTER :: Element
  REAL(KIND=dp) :: Norm, omega
  INTEGER :: i,j,k,n, nd, t, Active, iter, NonlinIter
  TYPE(ValueList_t), POINTER :: BC
  TYPE(Mesh_t),   POINTER :: Mesh
  COMPLEX(KIND=dp), PARAMETER :: im=(0._dp,1._dp)
  LOGICAL, SAVE :: NewtonRaphson = .FALSE., CSymmetry, DoRestart, &
      SliceAverage, RestartDone = .FALSE.
  INTEGER :: TransientSolverInd
  TYPE(Variable_t), POINTER :: CoordVar, LVar
  TYPE(ValueList_t), POINTER :: Params
  CHARACTER(LEN=MAX_NAME_LEN) :: sname   
  CHARACTER(*), PARAMETER :: Caller = 'MagnetoDynamics2DHarmonic'
    
!------------------------------------------------------------------------------

  Params => GetSolverParams()
  DoRestart = ListGetLogical( Params,'Transient Restart',Found )

  IF( DoRestart ) THEN
    ! IF we do restart, do it only once!
    IF( RestartDone ) RETURN    
  END IF
        
  CALL Info( Caller,'------------------------------------------------', Level=4 )
  CALL Info( Caller,'Solving equation for magnetic vector potential', Level=4 )
  CALL Info( Caller,'------------------------------------------------', Level=4 )
  
  CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )

  ! Allocate some permanent storage, this is done first time only:
  ! --------------------------------------------------------------
  Mesh => GetMesh()
  NULLIFY(BC)

  IF(GetCoupledIter() > 1) NewtonRaphson=.TRUE.

  ! Check whether we have also transient solver present.
  ! If we do, then add namespace for the material parameters.
  ! Note that these checks assume hard-coded subroutine names in this module.
  TransientSolverInd = 0
  DO i=1,Model % NumberOfSolvers      
    sname = GetString(Model % Solvers(i) % Values, 'Procedure', Found)
    j = INDEX( sname,'MagnetoDynamics2DHarmonic')
    IF( j > 0 ) CYCLE
    k = INDEX( sname,'MagnetoDynamics2D')
    IF( k > 0 ) THEN
      TransientSolverInd = i
      EXIT
    END IF
  END DO

    
  IF( TransientSolverInd > 0 ) THEN
    CALL Info(Caller,'Transient solver index found: '//I2S(i),Level=8)
    CALL ListPushNameSpace('harmonic:')
  ELSE IF( DoRestart ) THEN
    CALL Fatal(Caller,'Could not find transient solver for restart!')
  END IF

  Omega = GetAngularFrequency()
 
  ElectroDynamics = GetLogical( GetSolverParams(), 'Electrodynamics Model', Found)
  NonlinIter = GetInteger(Params,'Nonlinear system max iterations',Found)
  IF(.NOT.Found) NonlinIter = 1

  CALL DefaultStart()
  
  DO iter=1,NonlinIter

    IF(Iter>1) NewtonRaphson=.TRUE.
    ! System assembly:
    ! ----------------
    Active = GetNOFActive()
    CALL DefaultInitialize()
!$omp parallel do private(Element,n,nd)
    DO t=1,active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes(Element)
      nd = GetElementNOFDOFs(Element)
      CALL LocalMatrix(Element, n, nd)
    END DO
!$omp end parallel do
    CALL DefaultFinishBulkAssembly()

    Active = GetNOFBoundaryElements()
!$omp parallel do private(Element, n, nd, BC, Found)
    DO t=1,active
      Element => GetBoundaryElement(t)
      BC=>GetBC(Element)
      IF(.NOT.ASSOCIATED(BC)) CYCLE
      
      n  = GetElementNOFNodes(Element)
      nd = GetElementNOFDOFs(Element)
      
      IF(GetLogical(BC,'Infinity BC',Found)) THEN
        CALL LocalMatrixInfinityBC(  Element, n, nd )
      ELSE IF(GetLogical(BC,'Air Gap',Found)) THEN
        CALL LocalMatrixAirGapBC(Element, BC, n, nd)
      ELSE IF( ListCheckPresent( BC,'Layer Electric Conductivity' ) ) THEN
        CALL LocalMatrixSkinBC(Element, BC, n, nd)
      END IF
    END DO
!$omp end parallel do

    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
    
    CALL SetMagneticFluxDensityBC()
    CALL DefaultDirichletBCs()
    Norm = DefaultSolve()

    IF( DefaultConverged() ) EXIT
  END DO
  
  IF(.NOT. CSymmetry ) THEN
    CALL CalculateLumpedHarmonic()
  END IF

  CoordVar => VariableGet(Mesh % Variables,'Coordinates')
  IF(ASSOCIATED(CoordVar)) THEN
    DO i=1,Mesh % NumberOfNodes
      j = 3*(CoordVar % Perm(i)-1)
      CoordVar % Values(j+1) = Mesh % Nodes % x(i)
      CoordVar % Values(j+2) = Mesh % Nodes % y(i)
      CoordVar % Values(j+3) = Mesh % Nodes % z(i)
    END DO
  END IF
   
  CALL DefaultFinish()
  
  ! Perform restart if continuing to transient real-values combination. 
  IF( DoRestart ) THEN
    LVar => Model % Solvers(TransientSolverInd) % Variable 
    IF( ASSOCIATED( LVar ) ) THEN         
      LVar % Values = Solver % Variable % Values(1::2)
      LVar % PrevValues(:,1) = LVar % Values
    END IF

    sname = LagrangeMultiplierName( Solver )
    Lvar => VariableGet( Mesh % Variables, sname, ThisOnly = .TRUE. )
    IF ( ASSOCIATED(Lvar) ) THEN
      CALL Info(Caller,&
          'Size of Lagrange Multiplier: '//I2S(SIZE(LVar % Values)),Level=8)
      DO i=1,SIZE( LVar % Values ) / 2
        Lvar % Values(i) = Lvar % Values(2*(i-1)+1)
      END DO
    ELSE
      CALL Info(Caller,'Could not find Lagrange Multiplier for restart: '//TRIM(sname))
    END IF
    CALL Info(Caller,'Harmonic solution provided as initial guess for transient system!')
    RestartDone = .TRUE.
  END IF

  IF( TransientSolverInd > 0 ) THEN
    CALL ListPopNamespace()
  END IF

  
CONTAINS


!------------------------------------------------------------------------------
! This is monolithic lumping routine that has been optimized for speed.
! This way we need to evaluate the Basis functions only once for each element.
! This is modified version of the transient solver, this solver did not really
! need the minor optimizations.
!------------------------------------------------------------------------------
 SUBROUTINE CalculateLumpedHarmonic()
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: torq,TorqArea,IMoment,IA,Omega, &
       rinner,router,rmean,rdiff,ctorq,detJ,Weight,x,y,r,rho
   REAL(KIND=dp), ALLOCATABLE :: a(:),POT(:,:),Density(:)
   COMPLEX(KIND=dp), ALLOCATABLE :: POTC(:),U(:)
   REAL(KIND=dp), POINTER :: Basis(:), dBasisdx(:,:)
   LOGICAL, ALLOCATABLE :: TorqueElem(:)
   INTEGER :: i,bfid,n,nd,nbf,NoSlices,NoTimes,PrevComm
   LOGICAL :: Found, Stat
   TYPE(ValueList_t),POINTER::Params
   TYPE(GaussIntegrationPoints_t) :: IP
   TYPE(Nodes_t) :: Nodes
   LOGICAL :: CalcTorque, Calcpot, CalcInert
   LOGICAL :: ThisTorque, ThisPot, ThisInert, Parallel, HaveRange
   LOGICAL :: Visited = .FALSE.
   
   SAVE Visited, Nodes, Basis, dBasisdx, a, u, POT, POTC, Density, Ctorq, TorqueElem
   
!------------------------------------------------------------------------------

   CALL Info(Caller,'Calculating lumped parameters',Level=8)
   
   NoSlices = MAX(1,ListGetInteger( Model % Simulation,'Number Of Slices', SliceAverage ) )
   NoTimes = ListGetInteger( Model % Simulation,'Number Of Times', Found )
   IF( NoTimes > 1 ) THEN
     PrevComm = ParEnv % ActiveComm
     ParEnv % ActiveComm = ParallelSlicesComm() 
   END IF
     
   ! Define whether we have something to compute
   rinner = ListGetCRealAnyBody( Model,'r inner',CalcTorque )
   IF( CalcTorque ) THEN
     router = ListGetCRealAnyBody( Model,'r outer')
     rmean = (rinner+router)/2
     rdiff = (router-rinner)
     HaveRange = .TRUE.     
   ELSE
     rmean = ListGetConstReal( CurrentModel % Simulation,'Rotor Radius',CalcTorque)
     rmean = ParallelReduction( rmean, 2 )
     IF(.NOT. CalcTorque .AND. rmean > EPSILON(rmean) ) THEN
       CALL ListAddConstReal( CurrentModel % Simulation,'Rotor Radius',rmean)
       CalcTorque = .TRUE.
     END IF
     rdiff = ListGetConstReal( CurrentModel % Simulation,'Rotor Air Gap Width',Found)
     IF(.NOT. Found ) rdiff = 1.0e-3 * rmean
     HaveRange = .FALSE.
   END IF
   
   CalcPot = ListGetLogicalAnyBodyForce( Model,'Calculate Potential' )
   CalcInert = CalcTorque .AND. .NOT. Visited 
   
   IF(.NOT. (CalcTorque .OR. CalcPot .OR. CalcInert) ) RETURN

   Parallel = ( ParEnv % PEs > 1 ) 
   
   nbf = Model % NumberOfBodyForces
   IF(.NOT. Visited ) THEN
     n = Model % Mesh % MaxElementDofs
     ALLOCATE( a(nbf), u(nbf), POT(2,n), POTC(n), Density(n), &
         Basis(n), dBasisdx(n,3) )
   END IF

   IF( CalcTorque ) THEN
     torq = 0._dp
     TorqArea = 0._dp
   END IF
   IF( CalcInert ) THEN
     IMoment = 0._dp
     IA = 0.0_dp
   END IF
   IF( CalcPot ) THEN   
     U=0._dp
     a=0._dp
   END IF


   IF(.NOT. Visited .AND. CalcTorque ) THEN
     ALLOCATE( TorqueElem( GetNOFActive() ) )
     TorqueElem = .FALSE.
     
     DO i=1,GetNOFActive()
       Element => GetActiveElement(i)
     
       ThisTorque = .FALSE.
       
       n  = GetElementNOFNodes(Element)     
       CALL GetElementNodes( Nodes, Element )

       IF( HaveRange ) THEN       
         ! We are given range in classical Arkkio style. 
         ! Check how the center lies with respect to the range.
         x = SUM(Nodes % x(1:n))/n
         y = SUM(Nodes % y(1:n))/n
         r = SQRT(x**2+y**2)
         IF (r >= rinner .AND. r <= router) THEN
           TorqueElem(i) = .TRUE.
         END IF
       ELSE
         ! We are not given a range. Just take any element
         ! which has even one node at the given radius. 
         DO j=1,n
           x = Nodes % x(j)
           y = Nodes % y(j)
           r = SQRT(x**2+y**2)
           IF( ABS(r-rmean) < rdiff / 2 ) THEN
             TorqueElem(i) = .TRUE.
             EXIT
           END IF
         END DO
       END IF
     END DO
              
     i = COUNT( TorqueElem )
     CALL Info(Caller,'Number of elements to compute torque: '//I2S(i))
   END IF

   
   DO i=1,GetNOFActive()
     Element => GetActiveElement(i)

     nd = GetElementNOFDOFs(Element)
     n  = GetElementNOFNodes(Element)     
    
     ThisTorque = .FALSE.
     ThisPot = .FALSE.
     ThisInert = .FALSE.
     
     IF( CalcPot ) THEN
       Params => GetBodyForce(Element)
       IF(ASSOCIATED(Params)) THEN         
         ThisPot = GetLogical(Params,'Calculate Potential',Found)
         IF( ThisPot ) THEN
           bfid = GetBodyForceId(Element)           
         END IF
       END IF
     END IF
     
     IF( CalcInert ) THEN
       Params=>GetBodyParams(Element)
       IF(ASSOCIATED(Params)) THEN
         ThisInert = GetLogical(Params,'Calculate Inertial Moment',Found)
       END IF
       Density(1:n) = GetReal(GetMaterial(),'Density',Found,Element)
     END IF
     
     IF( CalcTorque ) ThisTorque = TorqueElem(i)

     ! Only treat the element if we have something to compute
     IF( .NOT. (ThisPot .OR. ThisInert .OR. ThisTorque ) ) CYCLE
     
     IF( ThisTorque .OR. ThisPot ) THEN
       CALL GetLocalSolution(POT, UElement=Element)
       POTC(1:nd) = POT(1,1:nd)+im*POT(2,1:nd)
     END IF
            
     CALL GetElementNodes( Nodes, Element )
     
     ! Numerical integration:
     !-----------------------
     IP = GaussPoints(Element)
     
     DO t=1,IP % n
       
       ! Basis function values & derivatives at the integration point:
       !--------------------------------------------------------------
       IF( ThisTorque ) THEN
         ! Only torque needs the derivatives of basis function
         stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
             IP % W(t), detJ, Basis, dBasisdx )
         weight = IP % s(t) * detJ
       ELSE
         stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
             IP % W(t), detJ, Basis )
         weight = IP % s(t) * detJ
       END IF

       ! Coordinates of the integration point
       x = SUM(Nodes % x(1:n)*Basis(1:n))
       y = SUM(Nodes % y(1:n)*Basis(1:n))
       r = SQRT(x**2+y**2)
       
       IF(ThisPot ) THEN
         Omega = GetAngularFrequency(UElement=Element)
         A(bfid) = A(bfid) + Weight
         U(bfid) = U(bfid) + Weight * im * Omega * SUM(POTC(1:nd)*Basis(1:nd))
       END IF
       
       IF( ThisTorque ) THEN                      
         BLOCK
           REAL(KIND=dp) :: BrRe,BpRe,BrIm,BpIm
           COMPLEX(KIND=dp) :: Bp,Br,Bx,By
           
           Bx =  SUM(POTC(1:nd)*dBasisdx(1:nd,2))
           By = -SUM(POTC(1:nd)*dBasisdx(1:nd,1))
           Br =  x/r*Bx + y/r*By
           Bp = -y/r*Bx + x/r*By
           
           BrRe = REAL( Br ); BrIm = AIMAG( Br )
           BpRe = REAL( Bp ); BpIm = AIMAG( Bp )

           Torq = Torq + weight * r * (BrRe*BpRe+BrIm*BpIm)/(2*PI*4.0d-7*rdiff)
           TorqArea = TorqArea + Weight
         END BLOCK
       END IF
       
       IF( ThisInert ) THEN
         IF( r < rmean ) THEN
           rho = SUM( density(1:n) * Basis(1:n) ) 
           IF( rho > EPSILON( rho ) ) THEN
             IA = IA + Weight
             U = U + Weight * r * rho
           END IF
         END IF
       END IF
     END DO
   END DO
     

   ! Finally perform parallel reduction if needed, and
   ! store the results for saving by SaveScalars.
   !-------------------------------------------------------------------------   
   IF( CalcPot ) THEN
     IF( Parallel ) THEN
       DO i=1,nbf
         a(i) = ParallelReduction(a(i)) / NoSlices
         u(i) = ParallelReduction(u(i)) / NoSlices
       END DO
     END IF     
     DO i=1,nbf
       IF(a(i)>0) THEN
         CALL ListAddConstReal(Model % Simulation,'res: Potential re / bodyforce ' &
             //i2s(i),REAL(u(i))/a(i))
         CALL ListAddConstReal(Model % Simulation,'res: Potential im / bodyforce ' &
             //i2s(i),AIMAG(u(i))/a(i))
         CALL ListAddConstReal(Model % Simulation,'res: area / bodyforce ' &
             //i2s(i),a(i)) 
       END IF
     END DO
   END IF
   
   IF( CalcTorque ) THEN   
     ! Arkkios formula assumes that rinner and router are nicely aligned with elements.
     ! This may not the case, so the 1st time we make a geometric correction. 
     IF(.NOT. Visited ) THEN
       WRITE(Message,'(A,ES15.4)') 'Air gap initial torque:', Torq
       CALL Info(Caller,Message,Level=6)

       TorqArea = ParallelReduction(TorqArea) / NoSlices
       IF (TorqArea > EPSILON(TorqArea) ) THEN
         Ctorq = 2 * PI * rmean * rdiff / TorqArea

         WRITE(Message,'(A,F8.4)') 'Air gap correction initial:', cTorq
         CALL Info(Caller,Message,Level=4)
         
         ! The correction factor also corrects for the number of periods.
         ! We don't want that - so let us take back that and the torque
         ! can be compared to inertial moment of the sector still. 
         i = ListGetInteger( CurrentModel % Simulation,'Rotor Periods',Found )
         IF( Parallel ) i = ParallelReduction( i, 2 ) 
         IF( i > 1 ) THEN
           WRITE(Message,'(A,I0)') 'Air gap correction rotor periods: ',i
           CALL Info(Caller,Message,Level=4)
           Ctorq = Ctorq / i 
         END IF         
       ELSE
         Ctorq = 1.0_dp
       END IF
       
       WRITE(Message,'(A,F8.4)') 'Air gap correction:', cTorq
       CALL Info(Caller,Message,Level=4)
       !CALL ListAddConstReal(Model % Simulation,'res: air gap correction', cTorq)
     END IF
       
     Torq = Ctorq * Torq

     IF( SliceAverage ) THEN
       ! Save slice torque even for one slice since then the output for scalars is the same
       ! for any number of slices.       
       WRITE(Message,'(A,ES15.4)') 'Air gap torque for slice'//I2S(ParEnv % MyPe)//':', Torq
       CALL Info(Caller,Message,Level=5)
       CALL ListAddConstReal(Model % Simulation,'res: air gap torque for slice', Torq)
     END IF
       
     ! But the averaging makes sense only for more than one slice
     IF(Parallel) Torq = ParallelReduction(Torq) / NoSlices

     WRITE(Message,'(A,ES15.4)') 'Air gap torque:', Torq
     CALL Info(Caller,Message,Level=5)
     CALL ListAddConstReal(Model % Simulation,'res: air gap torque', Torq)
   END IF
   
   IF( CalcInert ) THEN
     IF( Parallel ) THEN
       IMoment = ParallelReduction(IMoment) / NoSlices
       IA = ParallelReduction(IA) / NoSlices
     END IF

     WRITE(Message,'(A,ES15.4)') 'Inertial volume:', IA
     CALL Info(Caller,Message,Level=7)

     WRITE(Message,'(A,ES15.4)') 'Inertial moment:', Imoment
     CALL Info(Caller,Message,Level=7)
       
     CALL ListAddConstReal(Model % Simulation,'res: inertial volume', IA)
     CALL ListAddConstReal(Model % Simulation,'res: inertial moment', IMoment)
   END IF
     
   Visited = .TRUE.
  
   ! Revert the communicatior back to original
   IF( NoTimes > 1 ) ParEnv % ActiveComm = PrevComm
      
!------------------------------------------------------------------------------
 END SUBROUTINE CalculateLumpedHarmonic
!------------------------------------------------------------------------------


  
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  Element, n, nd)
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: Material,  BodyForce
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(ValueList_t), POINTER :: CompParams

    COMPLEX(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LoadAtIp,&
        JAC(nd,nd),Agrad(3),Load(n),M(2,n),M_ip(2),POTC(nd), C(n), C_ip
    COMPLEX(KIND=dp) :: nu_tensor(2,2)
    COMPLEX(KIND=dp) :: R(2,2,n)
    COMPLEX(KIND=dp) :: Bt(nd,2)
    COMPLEX(KIND=dp) :: Ht(nd,2) 
    COMPLEX(KIND=dp) :: B_ip(2), Alocal
    COMPLEX(KIND=dp) :: FR

    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,x,y
    REAL(KIND=dp) :: POT(2,nd),Babs,mu,muder,Omega
    REAL(KIND=dp) :: nu_11(nd), nuim_11(nd), nu_22(nd), nuim_22(nd)
    REAL(KIND=dp) :: nu_val, nuim_val
    REAL(KIND=dp) :: foilthickness, coilthickness, nofturns, skindepth, mu0 
    REAL(KIND=dp) :: Lorentz_velo(3,nd), Velo(3), omega_velo
    REAL(KIND=dp) :: LondonLambda_ip, P_ip
    REAL(KIND=dp) :: LondonLambda(nd), Permittivity(nd)

    INTEGER :: i,p,q,t

    LOGICAL :: HBcurve, Found, Stat, StrandedHomogenization
    LOGICAL :: CoilBody    
    LOGICAL :: InPlaneProximity = .FALSE., WithVelocity, WithAngularVelocity
    LOGICAL :: FoundIm, StrandedCoil
    LOGICAL :: LondonEquations = .TRUE.
    
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilType

    !$omp threadprivate(Nodes,InPlaneProximity)
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes,Element )
    STIFF = 0._dp
    JAC  = 0._dp
    FORCE = 0._dp
    IF(Transient) MASS = 0._dp

    Material => GetMaterial(Element)

    Omega = GetAngularFrequency(UElement=Element)
    
    InPlaneProximity = .FALSE.
    
    CoilBody = .FALSE.
    CompParams => GetComponentParams( Element )
    StrandedHomogenization = .FALSE.
    StrandedCoil = .FALSE.
    LondonEquations = .TRUE.
    IF (ASSOCIATED(CompParams)) THEN
      CoilType = GetString(CompParams, 'Coil Type', Found)
      IF (Found) THEN
        CoilBody = .TRUE.
        SELECT CASE (CoilType)
        CASE ('stranded')
          StrandedCoil = .TRUE.
          StrandedHomogenization = GetLogical(CompParams, 'Homogenization Model', Found)

          IF ( StrandedHomogenization ) THEN
            nu_11 = 0._dp
            nuim_11 = 0._dp
            nu_11 = GetReal(CompParams, 'nu 11', Found)
            nuim_11 = GetReal(CompParams, 'nu 11 im', FoundIm)
            IF ( .NOT. Found .AND. .NOT. FoundIm ) CALL Fatal ('LocalMatrix', 'Homogenization Model nu 11 not found!')
            nu_22 = 0._dp
            nuim_22 = 0._dp
            nu_22 = GetReal(CompParams, 'nu 22', Found)
            nuim_22 = GetReal(CompParams, 'nu 22 im', FoundIm)
            IF ( .NOT. Found .AND. .NOT. FoundIm ) CALL Fatal ('LocalMatrix', 'Homogenization Model nu 22 not found!')
          END IF

        CASE ('massive')
          LondonEquations = ListGetLogical(CompParams, 'London Equations', LondonEquations)
        CASE ('foil winding')
  !         CALL GetElementRotM(Element, RotM, n)
          InPlaneProximity = GetLogical(CompParams, 'Foil In Plane Proximity', Found)
          IF (InPlaneProximity) THEN
            coilthickness = GetConstReal(CompParams, 'Coil Thickness', Found)
            IF (.NOT. Found ) Call Fatal('LocalMatrix', 'Coil Thickness not found!')
            nofturns = GetConstReal(CompParams, 'Number Of Turns', Found)
            IF (.NOT. Found ) Call Fatal('LocalMatrix', 'Number of Turns not found!')
            foilthickness = coilthickness/nofturns
          END IF
        CASE DEFAULT
          CALL Fatal (Caller, 'Non existent Coil Type Chosen!')
        END SELECT
      END IF
    END IF

    IF ( LondonEquations ) THEN
      ! Lambda = m/(n_s * e^2):
      ! -------------------
      LondonLambda(:) = GetReal( Material, 'London Lambda', LondonEquations, Element)
    END IF

    HBCurve = ListCheckPresent(Material,'H-B Curve')
    IF(HBcurve) THEN
      CALL GetLocalSolution(POT,UElement=Element)
      POTC=POT(1,:)+im*POT(2,:)
    ELSE IF (.NOT. StrandedHomogenization) THEN 
      CALL GetReluctivity(Material,R,n,Element)
    END IF
 
    C = GetReal( Material, 'Electric Conductivity', Found, Element)
    C = C + im * GetReal( Material, 'Electric Conductivity im', Found, Element)

    M(1,:) = GetReal( Material, 'Magnetization 1', Found, Element)
    M(1,:) = M(1,:) + im*GetReal( Material, 'Magnetization 1 im', Found, Element)

    M(2,:) = GetReal( Material, 'Magnetization 2', Found, Element)
    M(2,:) = M(2,:) + im*GetReal( Material, 'Magnetization 2 im', Found, Element)

    IF(ElectroDynamics) THEN 
      Permittivity(1:n) = GetReal(Material, 'Permittivity', Found)
    END IF

    Load = 0.0d0
    WithVelocity = .FALSE.
    WithAngularVelocity = .FALSE.
    
    BodyForce => GetBodyForce(Element)
    IF ( ASSOCIATED(BodyForce) ) THEN
      Load(1:n) = GetReal( BodyForce, 'Current Density', Found, Element )
      Load(1:n) = Load(1:n) + im*GetReal( BodyForce, 'Current Density im', Found, Element )
      CALL GetRealVector(BodyForce, Lorentz_velo, 'Lorentz velocity', WithVelocity)
      omega_velo = ListGetCReal(BodyForce,'Angular velocity', WithAngularVelocity) 
    END IF

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
       IP % W(t), detJ, Basis, dBasisdx )

      IF( CSymmetry ) THEN
        x = SUM( Basis(1:n) * Nodes % x(1:n) )
        detJ = detJ * x
      END IF

      ! The source term at the integration point:
      !------------------------------------------
      LoadAtIP = SUM( LOAD(1:n)*Basis(1:n) )
      nu_tensor = 0.0_dp
      IF (HBcurve) THEN
        Agrad = 0.0_dp
        Agrad = MATMUL( POTC,dBasisdx )
        Alocal = SUM( POTC(1:n) * Basis(1:n) )
        ! Sign?
        ! -----
        B_ip(1) = -Agrad(2) 
        B_ip(2) = Agrad(1)
        IF( CSymmetry ) B_ip(2) = B_ip(2) + Alocal/x
        ! -----
        Babs = MAX(SQRT(SUM(ABS(B_ip)**2)), 1.d-8)

        IF( NewtonRaphson ) THEN
          mu = ListGetFun( Material,'h-b curve',babs,dFdx=muder) / Babs
          muder = (muder-mu)/babs
        ELSE
          mu = ListGetFun( Material,'h-b curve',babs) / Babs
        END IF

        nu_tensor(1,1) = mu ! Mu is really nu!!! too lazy to correct now...
        nu_tensor(2,2) = mu
      ELSE
        muder=0._dp
        IF (StrandedHomogenization) THEN
          nu_val = SUM( Basis(1:n) * nu_11(1:n) ) 
          nuim_val = SUM( Basis(1:n) * nuim_11(1:n) ) 
          nu_tensor(1,1) = CMPLX(nu_val, nuim_val, KIND=dp)
          nu_val = SUM( Basis(1:n) * nu_22(1:n) ) 
          nuim_val = SUM( Basis(1:n) * nuim_22(1:n) ) 
          nu_tensor(2,2) = CMPLX(nu_val, nuim_val, KIND=dp)
        ELSE 
          DO p=1,2
            DO q=1,2
              nu_tensor(p,q) = SUM(Basis(1:n) * R(p,q,1:n))
            END DO
          END DO
        END IF 
     END IF


      C_ip = SUM( Basis(1:n) * C(1:n) )
      M_ip = MATMUL( M,Basis(1:n) )
      P_ip = SUM( Basis(1:n)*Permittivity(1:n) )

      IF(.NOT. StrandedCoil ) THEN
        DO p=1,nd
          DO q=1,nd
            IF(ElectroDynamics) THEN
              STIFF(p,q) = STIFF(p,q) - &
                  IP % s(t) * detJ * omega**2 * P_ip * Basis(q)*Basis(p)
            ELSE
              STIFF(p,q) = STIFF(p,q) + &
                  IP % s(t) * detJ * im * omega * C_ip * Basis(q)*Basis(p)
            END IF
          END DO
        END DO
      END IF

      IF ( LondonEquations ) THEN
        LondonLambda_ip = SUM( Basis(1:n) * LondonLambda(1:n) )
        DO p=1,nd
          DO q=1,nd
            ! (LondonLambda a,  a'):
            ! --------------
            STIFF(p,q) = STIFF(p,q) + IP % s(t) * detJ / LondonLambda_ip * Basis(q)*Basis(p)
          END DO
        END DO
      END IF

      Bt(:,1) = -dbasisdx(:,2)
      Bt(:,2) =  dbasisdx(:,1)
      IF ( CSymmetry ) Bt(:,2) = Bt(:,2) + Basis(:)/x
      
      IF (InPlaneProximity) THEN
        FR = 0._dp + im*0._dp
        mu0 = 4d-7 * pi
        skindepth = sqrt(2._dp/(omega * C_ip * mu0))
        FR = C_ip * foilthickness * skindepth * omega * (1_dp + im)/8._dp
        FR = FR*(-im)*SIN(im*(1_dp+im)*foilthickness/skindepth)
        FR = FR/(-im * SIN(im*(1_dp+im)*foilthickness/skindepth/2._dp))**2._dp
        nu_tensor(1,1) = nu_tensor(1,1) + FR - 1._dp/mu0
        nu_tensor(2,2) = nu_tensor(2,2) + FR - 1._dp/mu0
      END IF

      DO p = 1,nd
        Ht(p,:) = MATMUL(nu_tensor, Bt(p,:))
      END DO

      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + IP % s(t) * DetJ * &
             MATMUL(Ht, TRANSPOSE(Bt))

      IF (HBcurve .AND. NewtonRaphson) THEN
!        DO p=1,nd
!          DO q=1,nd
!            JAC(p,q) = JAC(p,q) + IP % s(t) * DetJ * &
!              muder/babs*SUM(Agrad*dBasisdx(q,:))*SUM(CONJG(Agrad)*dBasisdx(p,:))
!          END DO
!        END DO
        DO p=1,nd
          DO q=1,nd
            JAC(p,q) = JAC(p,q) + IP % s(t) * DetJ * &
              muder/babs*SUM(B_ip*Bt(q,:))*SUM(CONJG(B_ip)*Bt(p,:))
          END DO
        END DO
      END IF

      IF (WithVelocity .OR. WithAngularVelocity ) THEN
        !
        ! Create an additional Lorentz effect so that the electric field
        ! has an added term v x curl A:
        !
        IF( WithVelocity ) THEN
          Velo(1:2) = [ SUM(Basis(1:n)*Lorentz_velo(1,1:n)), &
              SUM(Basis(1:n)*Lorentz_velo(2,1:n)) ]
        ELSE 
          x = SUM( Basis(1:n) * Nodes % x(1:n) )
          y = SUM( Basis(1:n) * Nodes % y(1:n) ) 

          ! Simplified omega \times r in 2D
          Velo(1) = -omega_velo * y
          Velo(2) = omega_velo * x
        END IF

        DO p=1,nd
          STIFF(p,1:nd) = STIFF(p,1:nd) + IP % s(t) * DetJ * C_ip * Basis(p) * ( & 
              -Velo(2) * Bt(1:nd,1) + Velo(1) * Bt(1:nd,2) )
        END DO
      END IF

      IF ( CSymmetry ) THEN
        FORCE(1:nd) = FORCE(1:nd) + IP % s(t) * DetJ * (LoadAtip * Basis(1:nd) - &
            M_ip(1)*dBasisdx(1:nd,2)+M_ip(2)*(dBasisdx(1:nd,1) + Basis(1:nd)/x))
      ELSE
        FORCE(1:nd) = FORCE(1:nd) + IP % s(t) * DetJ * (LoadAtip * Basis(1:nd) + &
            M_ip(1)*dBasisdx(1:nd,2)-M_ip(2)*dBasisdx(1:nd,1))
      END IF
    END DO

    IF (HBcurve .AND. NewtonRaphson) THEN
      STIFF = STIFF + JAC
      FORCE = FORCE + MATMUL(JAC,POTC)
    END IF

    IF(Transient) THEN
      CALL Default1stOrderTime( MASS, STIFF, FORCE, UElement=Element )
    END IF
    CALL DefaultUpdateEquations( STIFF, FORCE, UElement=Element )

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixInfinityBC(Element, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),DetJ
    LOGICAL :: Stat
    INTEGER :: i,p,q,t
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: Inf_ip,Coord(3),Normal(3),mu,u,v
    COMPLEX(KIND=dp) :: R(2,2,n)       
    COMPLEX(KIND=dp) :: STIFF(nd,nd), FORCE(nd)
    TYPE(ValueList_t), POINTER :: Material
    TYPE(Element_t), POINTER :: Parent
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !$OMP THREADPRIVATE(Nodes)
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes, Element )
    STIFF = 0._dp
    FORCE = 0._dp

    Parent=>Element % BoundaryInfo % Left
    IF(.NOT.ASSOCIATED(Parent)) THEN
      Parent=>Element % BoundaryInfo % Right
    END IF
    IF(.NOT.ASSOCIATED(Parent)) RETURN

    Material => GetMaterial(Parent)
    CALL GetReluctivity(Material,R,n,Parent)

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                 IP % W(t), detJ, Basis )

      mu = SUM(Basis(1:n)*R(1,1,1:n)) !We assume isotropic permeability

      Normal = NormalVector( Element, Nodes, u, v, .TRUE. )
      Coord(1) = SUM(Basis(1:n) * Nodes % x(1:n))
      Coord(2) = SUM(Basis(1:n) * Nodes % y(1:n))
      Coord(3) = SUM(Basis(1:n) * Nodes % z(1:n))
      
      IF( CSymmetry ) THEN
        detJ = detJ * Coord(1)
      END IF

      Inf_ip = mu * SUM(Coord*Normal)/SUM(Coord*Coord)

      DO p=1,nd
        DO q=1,nd
          STIFF(p,q) = STIFF(p,q) + IP % s(t)*detJ*Inf_ip*Basis(q)*Basis(p)
        END DO
      END DO
    END DO
    CALL DefaultUpdateEquations( STIFF, FORCE, UElement=Element )
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixInfinityBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixAirGapBC(Element, BC, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ
    LOGICAL :: Stat, Found
    INTEGER :: i,t
    TYPE(GaussIntegrationPoints_t) :: IP
    COMPLEX(KIND=dp) :: STIFF(nd,nd), FORCE(nd)
    REAL(KIND=dp) :: mu,x,AirGapLength(nd), AirGapMu(nd), AirGapL
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !$OMP THREADPRIVATE(Nodes)
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes, Element )
    STIFF = 0._dp
    FORCE = 0._dp

    AirGapLength=GetConstReal( BC, 'Air Gap Length', Found)
    IF (.NOT. Found) CALL Fatal('LocalMatrixAirGapBC', 'Air Gap Length not found!')

    AirGapMu=GetConstReal( BC, 'Air Gap Relative Permeability', Found)
    IF (.NOT. Found) AirGapMu=1.0_dp

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )

      IF( CSymmetry ) THEN
        x = SUM( Basis(1:n) * Nodes % x(1:n) )
        detJ = detJ * x
      END IF
      
      mu = 4*pi*1d-7*SUM(Basis(1:n)*AirGapMu(1:n))
      AirGapL = SUM(Basis(1:n)*AirGapLength(1:n))

      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + IP % s(t) * DetJ * &
          AirGapL/mu*MATMUL(dBasisdx, TRANSPOSE(dBasisdx))
      
    END DO
    CALL DefaultUpdateEquations( STIFF, FORCE, UElement=Element )
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixAirGapBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixSkinBC(Element, BC, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,x
    LOGICAL :: Stat, Found
    INTEGER :: i,p,q,t
    TYPE(GaussIntegrationPoints_t) :: IP
    COMPLEX(KIND=dp) :: STIFF(nd,nd), FORCE(nd), imu, invZs, delta
    REAL(KIND=dp) :: SkinCond(nd), Mu(nd), CondAtIp, MuAtIp, MuVacuum
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !$OMP THREADPRIVATE(Nodes)
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes, Element )
    STIFF = 0._dp
    FORCE = 0._dp

    muVacuum = 4 * PI * 1d-7
    imu = CMPLX(0.0_dp, 1.0_dp)
    
    SkinCond(1:n) = GetReal( BC,'Layer Electric Conductivity', Found)
    Mu(1:n) = GetReal( BC,'Layer Relative Permeability', Found)
      
    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )

      IF( CSymmetry ) THEN
        x = SUM( Basis(1:n) * Nodes % x(1:n) )
        detJ = detJ * x
      END IF
      
      muAtIP = muVacuum*SUM(Basis(1:n)*Mu(1:n))
      condAtIp = SUM(Basis(1:n)*SkinCond(1:n))

      delta = SQRT( 2.0_dp/(condAtIp*omega*muAtIp))      
      invZs = (condAtIp*delta)/(1.0_dp+imu)
      
      DO p=1,nd
        DO q=1,nd
          STIFF(p,q) = STIFF(p,q) + IP % s(t) * DetJ * &
              ( imu * omega * invZs ) * Basis(p) * Basis(q)
        END DO
      END DO
              
    END DO
    CALL DefaultUpdateEquations( STIFF, FORCE, UElement=Element )
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixSkinBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE SetMagneticFluxDensityBC()
!------------------------------------------------------------------------------
! P. Lombard, G. Meunier, "A general purpose method for electric and magnetic 
! combined problems for 2D, axisymmetric and transient systems", IEEE Trans.
! magn. 29(2), p. 1737 - 1740, Mar 1993
! -ettaka- 
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Matrix_t), POINTER :: A
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp), POINTER :: b(:)
    INTEGER :: i, n, j, k
    TYPE(ValueList_t), POINTER :: BC
    LOGICAL :: Found
    REAL(KIND=dp) :: Bx(Solver % Mesh % MaxElementDofs), &
                      Bxim(Solver % Mesh % MaxElementDofs), &
                      By(Solver % Mesh % MaxElementDofs), &
                      Byim(Solver % Mesh % MaxElementDofs)
    REAL(KIND=dp) :: x, y
    INTEGER, POINTER :: Perm(:)

    Perm => Solver % Variable % Perm
    A => Solver % Matrix
    b => A % RHS
    DO i=1,GetNofBoundaryElements()
      Element => GetBoundaryElement(i)
      n = GetELementNofNodes()
      BC => GetBC()
      IF ( ASSOCIATED(BC)) THEN
        IF ( ListCheckPrefix( BC, 'Magnetic Flux Density') ) THEN
          Bx = 0._dp
          Bxim = 0._dp
          By = 0._dp
          Byim = 0._dp

          Bx(1:n) = GetReal(BC, 'Magnetic Flux Density 1', Found)
          Bxim(1:n) = GetReal(BC, 'Magnetic Flux Density 1 im', Found)
          By(1:n) = GetReal(BC, 'Magnetic Flux Density 2', Found)
          Byim(1:n) = GetReal(BC, 'Magnetic Flux Density 2 im', Found)

          DO j = 1,n
            k = Element % NodeIndexes(j)
            x = Mesh % Nodes % x(k)
            y = Mesh % Nodes % y(k)
            k = Perm(k)

            CALL UpdateDirichletDof( A, 2*k-1, y * Bx(j) - x * By(j) )
            CALL UpdateDirichletDof( A, 2*k, y * Bxim(j) - x * Byim(j) )
          END DO 
        END IF  
      END IF  
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE SetMagneticFluxDensityBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE GetReluctivity(Material,Acoef,n,Element)
!------------------------------------------------------------------------------
    USE MGDynMaterialUtils
    TYPE(ValueList_t), POINTER :: Material
    INTEGER :: n
    COMPLEX(KIND=dp) :: Acoef(2,2,n)
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    LOGICAL :: Found
    REAL(KIND=dp), SAVE :: Avacuum
    LOGICAL, SAVE :: FirstTime = .TRUE.

    !$OMP THREADPRIVATE(FirstTime, Avacuum)

    IF ( FirstTime ) THEN
      Avacuum = GetConstReal( CurrentModel % Constants, &
              'Permeability of Vacuum', Found )
      IF(.NOT. Found ) Avacuum = PI * 4.0d-7
      FirstTime = .FALSE.
    END IF

    Acoef = GetCMPLXTensor(Element, n, 2, 'Relative Permeability', Found)
    
    IF ( Found ) THEN
      Acoef = Avacuum * Acoef
    ELSE
      Acoef = GetCMPLXTensor(Element, n, 2, 'Permeability', Found)
    END IF
    IF ( Found ) THEN
      Acoef = Get2x2CMPLXTensorInverse(Acoef, n)
    ELSE
      Acoef = GetCMPLXTensor(Element, n, 2, 'Reluctivity', Found)
    END IF
    
    IF( .NOT. Found ) THEN
      CALL Warn('GetReluctivity',&
          'Could not get either > Reluctivity > or > Relative Permeability < !')
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE GetReluctivity
!------------------------------------------------------------------------------

END SUBROUTINE MagnetoDynamics2DHarmonic
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Initialization for the primary solver: BSolver
!------------------------------------------------------------------------------
SUBROUTINE Bsolver_init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------

  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver  
  TYPE(Model_t) :: Model    
  REAL(KIND=dp) :: dt       
  LOGICAL :: Transient      
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(ValueList_t),POINTER :: SolverParams
  LOGICAL :: Found

  SolverParams => GetSolverParams()
  CALL ListAddNewString( SolverParams, 'Variable','-nooutput bsolver_temp' )
  IF( GetLogical( SolverParams,'Target Variable Complex',Found ) ) THEN
    CALL ListAddString( SolverParams,&
        NextFreeKeyword('Exported Variable',SolverParams),'B[B re:2 B im:2]')
  ELSE
    CALL ListAddString( SolverParams,&
        NextFreeKeyword('Exported Variable',SolverParams),'B[B:2]')
  END IF
  
  IF( ListGetLogical( SolverParams, 'Calculate Joule Heating', Found ) ) THEN
    CALL ListAddString( SolverParams, &
        NextFreeKeyword('Exported Variable',SolverParams), &
        'Joule Heating' )
    CALL ListAddString( SolverParams, &
        NextFreeKeyword('Exported Variable',SolverParams), &
        'Joule Field' )
    CALL ListAddString( SolverParams, &
        NextFreeKeyword('Exported Variable',SolverParams), &
        'Current Density[Current Density re:1 Current Density im:1]' )
  END IF

  ! The reference norm is sum of all solutions. Hence we don't really want to recompute and spoil it externally.
  CALL ListAddNewLogical( SolverParams,'Skip Compute Steady State Change',.TRUE.)
  

END SUBROUTINE Bsolver_init


!------------------------------------------------------------------------------
!> Given the vector potential compute its curl, i.e. the magnetic
!> flux density.  
!> NOTE: THIS IS OBSOLETE. It is recommended that the subroutine 
!> MagnetoDynamicsCalcFields within the module MagnetoDynamics is used for 
!> postprocessing.
!------------------------------------------------------------------------------
SUBROUTINE Bsolver( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------

  USE DefUtils
  USE CircuitUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver  !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model    !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt       !< Timestep size for time dependent simulations
  LOGICAL :: Transient      !< Steady state or transient simulation
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(ValueList_t),POINTER :: SolverParams
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName
  INTEGER :: i,j,k,dim,FluxDofs,TotDofs
  LOGICAL :: ConstantBulkMatrix, ConstantBulkMatrixInUse
  LOGICAL :: GotIt
  REAL(KIND=dp) :: Unorm, Totnorm
  REAL(KIND=dp), ALLOCATABLE, TARGET :: ForceVector(:,:)
  REAL(KIND=dp), POINTER CONTIG :: SaveRHS(:)  
  TYPE(Variable_t), POINTER :: FluxSol, HeatingSol, JouleSol, AzSol
  LOGICAL ::  CSymmetry, LossEstimation, JouleHeating, ComplexPowerCompute,&
              AverageBCompute, BodyICompute, BodyVolumesCompute = .FALSE., &
              CirCompVolumesCompute = .FALSE., HomogenizationParamCompute, &
              LorentzForceCompute = .FALSE.
  TYPE(Variable_t), POINTER :: CurrDensSol
  CHARACTER(*), PARAMETER :: Caller = 'BSolver'


  CALL Warn(Caller,'This module is obsolete! USE MagnetoDynamicsCalcFields instead')

  
  CALL Info( Caller, '-------------------------------------',Level=4 )
  CALL Info( Caller, 'Computing the magnetic field density ',Level=4 )
  CALL Info( Caller, '-------------------------------------',Level=4 )
  dim = CoordinateSystemDimension()

  CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )

!------------------------------------------------------------------------------
!  Check what needs to be computed
!------------------------------------------------------------------------------
  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
  IF ( COUNT( Solver % Variable % Perm > 0 ) <= 0 ) RETURN


  BLOCK
    INTEGER :: NoSlices, NoTimes
    NoSlices = ListGetInteger( Model % Simulation,'Number Of Slices', GotIt )
    NoTimes = ListGetInteger( Model % Simulation,'Number Of Times', GotIt )
    IF(NoSlices > 1 .OR. NoTimes > 1 ) THEN
      CALL Fatal(Caller,'BSolver cannot deal with Slices or Times, use CalcFields!')
    END IF
  END BLOCK
    
  
  SolverParams => GetSolverParams()

  VarName = GetString(GetSolverParams(),'Target Variable',GotIt)
  IF(.NOT. GotIt) VarName = 'Potential'
  AzSol => VariableGet( Solver % Mesh % Variables, VarName ) 
  IF( .NOT. ASSOCIATED( AzSol ) ) THEN
    CALL Fatal(Caller,'Target field not present: '//TRIM(VarName) )
  END IF

  FluxSol => VariableGet(Solver % Mesh % Variables, 'B')
  IF( .NOT. ASSOCIATED( FluxSol ) ) THEN
    CALL Fatal(Caller,'Solution field not present: B' )
  END IF

  IF( FluxSol % Dofs / AzSol % Dofs /= 2 ) THEN
    WRITE( Message,'(A,I0,A)') 'B field should have ',2 * AzSol % Dofs,' dofs!'
    CALL Fatal(Caller,Message)
  END IF


  FluxDofs = FluxSol % Dofs
  ConstantBulkMatrix = GetLogical( SolverParams, 'Constant Bulk Matrix', GotIt )
  ConstantBulkMatrixInUse = ConstantBulkMatrix .AND. &
      ASSOCIATED(Solver % Matrix % BulkValues)
  
  CALL DefaultInitialize(Solver, ConstantBulkMatrixInUse)
  
  TotDofs = FluxDofs
  JouleHeating = ListGetLogical( SolverParams, 'Calculate Joule Heating', GotIt )

  IF( JouleHeating ) THEN
    IF( FluxDofs /= 4 ) THEN
      CALL Fatal(Caller,'Joule heating can only be computed for complex problems!')
    ELSE
      TotDofs = TotDofs + 2
      HeatingSol => VariableGet(Solver % Mesh % Variables, 'Joule Heating')
      IF( .NOT. ASSOCIATED( HeatingSol ) ) THEN
        CALL Fatal(Caller,'Solution field not present: Joule Heating' )
      END IF
      JouleSol => VariableGet(Solver % Mesh % Variables, 'Joule Field')
      IF( .NOT. ASSOCIATED( JouleSol ) ) THEN
        CALL Fatal(Caller,'Solution field not present: Joule Field' )
      END IF
      TotDofs = TotDofs + 2
      CurrDensSol => VariableGet(Solver % Mesh % Variables, 'Current Density')
      IF( .NOT. ASSOCIATED( CurrDensSol ) ) THEN
        CALL Fatal(Caller,'Solution field not present: Current Density' )
      END IF
    END IF
  END IF

  !------------------------------------------------------------------------------
  ! In the case of time-harmonic analysis losses may be estimated in terms of B
  !------------------------------------------------------------------------------ 
  LossEstimation = GetLogical(SolverParams,'Loss Estimation',GotIt)
  IF( LossEstimation .AND. FluxDofs /= 4) THEN
    CALL Fatal( Caller, 'Real solution, loss estimation omitted' )
  END IF

  HomogenizationParamCompute = GetLogical(SolverParams, 'Calculate Homogenization Parameters', GotIt)
  IF (.NOT. GotIt ) HomogenizationParamCompute = .FALSE.
  IF( HomogenizationParamCompute.AND. FluxDofs /= 4) THEN
    CALL Fatal( Caller, 'Real solution, Calculate Homogenization Parameters omitted' )
  END IF

  ComplexPowerCompute = GetLogical(SolverParams,'Calculate Complex Power',GotIt)
  IF (.NOT. GotIt ) ComplexPowerCompute = .FALSE.
  IF( ComplexPowerCompute.AND. FluxDofs /= 4) THEN
    CALL Fatal( Caller, 'Real solution, Complex Power omitted' )
  END IF

  AverageBCompute = GetLogical(SolverParams, 'Calculate Average Magnetic Flux Density', GotIt)
  IF (.NOT. GotIt ) AverageBCompute = .FALSE.

  BodyICompute = GetLogical(SolverParams, 'Calculate Body Current', GotIt)
  IF (.NOT. GotIt ) BodyICompute = .FALSE.

  LorentzForceCompute = GetLogical(SolverParams, 'Calculate Component Lorentz Force', GotIt)
  IF (.NOT. GotIt ) LorentzForceCompute = .FALSE.

  ALLOCATE(ForceVector(SIZE(Solver % Matrix % RHS),TotDOFs))  
  ForceVector = 0.0_dp
  SaveRHS => Solver % Matrix % RHS

  CALL BulkAssembly()

  IF (ConstantBulkMatrix) THEN
    CALL DefaultFinishBulkAssembly(BulkUpdate = .NOT.ConstantBulkMatrixInUse, RHSUpdate = .FALSE.)
  ELSE
    CALL DefaultFinishBulkAssembly()
  END IF

  CALL DefaultFinishAssembly()
  CALL DefaultDirichletBCs()

  TotNorm = 0.0_dp
  DO i=1,TotDofs
    Solver % Matrix % RHS => ForceVector(:,i)
    Solver % Variable % Values = 0
    UNorm = DefaultSolve()
    TotNorm = TotNorm + Unorm**2
    IF( i <= FluxDofs ) THEN
      FluxSol % Values(i::FluxDofs) = Solver % Variable % Values
    ELSE IF( i == FluxDofs + 1 ) THEN
      JouleSol % Values = Solver % Variable % Values
    ELSE IF( i == FluxDofs + 2 ) THEN
      HeatingSol % Values = Solver % Variable % Values
    ELSE 
      CurrDensSol % Values(i-Fluxdofs-2::2) = Solver % Variable % Values
    END IF
  END DO
  DEALLOCATE( ForceVector )  

  Solver % Matrix % RHS => SaveRHS
  TotNorm = SQRT(TotNorm)
  Solver % Variable % Norm = Totnorm

!------------------------------------------------------------------------------     
  
  WRITE( Message, * ) 'Result Norm: ',TotNorm
  CALL Info( Caller, Message, Level=4 )



CONTAINS


!------------------------------------------------------------------------------
  SUBROUTINE BulkAssembly()
!------------------------------------------------------------------------------
       
    INTEGER :: elem,t,i,j,k,p,q,n,nd,BodyId
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: weight,coeff,detJ,BAtIp(8),PotAtIp(2),MuAtIp, &
        Omega,TotalHeating, DesiredHeating, HeatingCoeff
    COMPLEX(KIND=dp) :: CondAtIp
    REAL(KIND=dp) :: Freq, FreqPower, FieldPower, ComponentLoss(2), LossCoeff, &
        ValAtIp, ValAtIpim, TotalLoss, x
    LOGICAL :: Found
    TYPE(ValueList_t), POINTER :: Material

    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:,:)
    REAL(KIND=dp), ALLOCATABLE :: POT(:,:)
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
    REAL(KIND=dp), ALLOCATABLE :: Cond(:), mu(:)
    REAL(KIND=dp), ALLOCATABLE :: BodyLoss(:), BodyComplexPower(:,:), BodyCurrent(:,:), &
        CirCompComplexPower(:,:), CirCompCurrent(:,:), &
        BodyLorentzForcesRe(:,:), BodyLorentzForcesIm(:,:), &
        ComponentLorenzForcesRe(:,:), ComponentLorenzForcesIm(:,:)
    COMPLEX(KIND=dp) :: cmplx_power 
    REAL(KIND=dp), ALLOCATABLE :: BodyVolumes(:), BodyAvBim(:,:), BodyAvBre(:,:), &
        BodySkinCond(:,:), BodyProxNu(:,:), &
        CirCompVolumes(:), CirCompAvBim(:,:), CirCompAvBre(:,:), &
        CirCompSkinCond(:,:), CirCompProxNu(:,:) 
    LOGICAL, ALLOCATABLE :: BodyAverageBCompute(:)

    REAL(KIND=dp), ALLOCATABLE :: alpha(:)
    TYPE(Variable_t), POINTER :: LagrangeVar
    COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
    REAL(KIND=dp) :: localV(2), coilthickness, localAlpha, N_j
    TYPE(ValueList_t), POINTER :: CompParams
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilType, bodyNumber, XYNumber, str
    LOGICAL :: CoilBody, EddyLoss
    COMPLEX(KIND=dp) :: imag_value, imag_value2
    INTEGER :: IvarId, ReIndex, ImIndex, VvarDofs, VvarId
    REAL(KIND=DP) :: grads_coeff, nofturns
    REAL(KIND=DP) :: i_multiplier_re, i_multiplier_im, ModelDepth
    COMPLEX(KIND=dp) :: i_multiplier, Bx, By, Jz, LorentzForceDensX, &
        LorentzForceDensY

    INTEGER :: NofComponents=0, bid
    INTEGER, POINTER :: BodyIds(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: CompNumber, OutputComp

    LOGICAL :: StrandedHomogenization, FoundIm, StrandedCoil 

    REAL(KIND=dp), ALLOCATABLE :: sigma_33(:), sigmaim_33(:)
    REAL(KIND=dp), ALLOCATABLE :: CoreLossUDF(:)
    REAL(KIND=dp) :: CoreLossUDFatIp

    LOGICAL :: LaminateModelPowerCompute=.FALSE., InPlaneProximity=.FALSE.
    REAL(KIND=dp) :: LaminatePowerDensity, BMagnAtIP, Fsk, Lambda, LaminateThickness, &
        mu0=4d-7*PI, skindepth

    LOGICAL :: BertottiCompute = .FALSE., LossUDF = .FALSE.
    REAL(KIND=dp) :: BertottiLoss, BRTc1, BRTc2, BRTc3, BRTc4, BRTc5

    SAVE Nodes

    n = 2*MAX(Solver % Mesh % MaxElementDOFs,Solver % Mesh % MaxElementNodes)
    ALLOCATE( STIFF(n,n), FORCE(Totdofs,n) )
    ALLOCATE( POT(2,n), Basis(n), dBasisdx(n,3), alpha(n) )
    ALLOCATE( Cond(n), mu(n), sigma_33(n), sigmaim_33(n), CoreLossUDF(n)) 
    

    str = LagrangeMultiplierName( Azsol % Solver )
    LagrangeVar => VariableGet( Solver % Mesh % Variables, str, ThisOnly = .TRUE.)

    ModelDepth = GetCircuitModelDepth()

    IF( JouleHeating ) THEN
      Omega = GetAngularFrequency()
      TotalHeating = 0.0_dp
    END IF

    IF( LossEstimation ) THEN
      ALLOCATE( BodyLoss(Model % NumberOfBodies) )
      Freq = Omega / (2*PI)
      
      FreqPower = GetCReal( SolverParams,'Fourier Loss Frequency Exponent',Found )
      IF( .NOT. Found ) FreqPower = 2.0_dp
      
      FieldPower = GetCReal( SolverParams,'Fourier Loss Field Exponent',Found ) 
      IF( .NOT. Found ) FieldPower = 2.0_dp
      FieldPower = FieldPower / 2.0_dp

      ComponentLoss = 0.0_dp
      BodyLoss = 0.0_dp
    END IF

    IF (HomogenizationParamCompute) THEN
      NofComponents = SIZE(Model % Components)
      Omega = GetAngularFrequency()
      CALL ListAddConstReal( Model % Simulation, 'res: Angular Frequency', Omega)
      NofComponents = SIZE(Model % Components)
      ALLOCATE(BodySkinCond(2, Model % NumberOfBodies), &
          BodyProxNu(2, Model % NumberOfBodies), &
          CirCompSkinCond(2, Model % NumberOfBodies), &
          CirCompProxNu(2, Model % NumberOfBodies))
      BodySkinCond = 0.0_dp
      BodyProxNu = 0.0_dp      
      CirCompSkinCond = 0.0_dp
      CirCompProxNu = 0.0_dp
      BodyICompute = .TRUE.
      ComplexPowerCompute = .TRUE.
      AverageBCompute = .TRUE.
    END IF

    IF ( ComplexPowerCompute ) THEN
      NofComponents = SIZE(Model % Components)
      ALLOCATE( BodyComplexPower(2,Model % NumberOfBodies), &
          CirCompComplexPower(2, NofComponents ) )
      BodyComplexPower = 0.0_dp
      CirCompComplexPower = 0.0_dp
    END IF

    IF (BodyICompute) THEN
      NofComponents = SIZE(Model % Components)
      ALLOCATE(BodyCurrent(2, Model % NumberOfBodies))
      ALLOCATE(CirCompCurrent(2, Model % NumberOfBodies))
      BodyCurrent = 0.0_dp
      CirCompCurrent = 0.0_dp
    END IF

    IF (LorentzForceCompute) THEN
      NofComponents = SIZE(Model % Components)
      ALLOCATE(BodyLorentzForcesRe(2, Model % NumberOfBodies))
      ALLOCATE(BodyLorentzForcesIm(2, Model % NumberOfBodies))
      ALLOCATE(ComponentLorenzForcesRe(2, NofComponents))
      ALLOCATE(ComponentLorenzForcesIm(2, NofComponents))
      BodyLorentzForcesRe = 0.0_dp
      BodyLorentzForcesIm = 0.0_dp
      ComponentLorenzForcesRe = 0.0_dp
      ComponentLorenzForcesIm = 0.0_dp
    END IF

    IF ( AverageBCompute ) THEN
      NofComponents = SIZE(Model % Components)
      ALLOCATE( BodyAvBre(2,Model % NumberOfBodies), &
          BodyAvBim(2,Model % NumberOfBodies), &
          BodyAverageBCompute(Model % NumberOfBodies), &
          CirCompAvBre(2,NofComponents), &
          CirCompAvBim(2,NofComponents) )

      BodyAvBre = 0._dp
      BodyAvBim = 0._dp
      BodyVolumesCompute = .TRUE.        
      CirCompAvBre = 0.0_dp
      CirCompAvBim = 0.0_dp
      CirCompVolumesCompute = .TRUE.        

      DO i = 1, Model % NumberOfBodies
        BodyAverageBCompute(i) = ListGetLogical(Model % Bodies(i) % Values,&
            'Compute Average Magnetic Flux Density', Found)
        IF (.NOT. Found) BodyAverageBCompute(i) = .TRUE.
      END DO
    END IF 

    IF ( BodyVolumesCompute ) THEN
      ALLOCATE( BodyVolumes(Model % NumberOfBodies) )
      BodyVolumes = 0._dp
    END IF

    IF ( CirCompVolumesCompute ) THEN
      NofComponents = SIZE(Model % Components)
      ALLOCATE( CirCompVolumes(NofComponents) )
      CirCompVolumes = 0._dp
    END IF

    DO elem = 1,GetNOFActive()
         
      ! Element information
      ! ---------------------
      Element => GetActiveElement(elem)
      CALL GetElementNodes( Nodes )
      nd = GetElementNOFDOFs()
      n  = GetElementNOFNodes()
      
      CompParams => GetComponentParams( Element )
      StrandedHomogenization = .FALSE.
      InPlaneProximity = .FALSE.
      LaminateModelPowerCompute = .FALSE.
      StrandedCoil = .FALSE.
      CoilType = ''
      
      IF (ASSOCIATED(CompParams)) THEN    
        CoilType = GetString(CompParams, 'Coil Type', Found)
        IF (Found) CoilBody = .TRUE.
        
        SELECT CASE (CoilType)
        CASE ('stranded')
          StrandedCoil = .TRUE.
          
          IvarId = GetInteger (CompParams, 'Circuit Current Variable Id', Found)
          IF (.NOT. Found) CALL Fatal (Caller, 'Circuit Current Variable Id not found!')
 
          N_j = GetConstReal (CompParams, 'Stranded Coil N_j', Found)
          IF (.NOT. Found) CALL Fatal (Caller, 'Stranded Coil N_j not found!')
 
          !nofturns = GetConstReal(CompParams, 'Number of Turns', Found)
          !IF (.NOT. Found) CALL Fatal(Caller,'Stranded Coil: Number of Turns not found!')
          
          i_multiplier_re = GetConstReal(CompParams, 'Current Multiplier re', Found)
          i_multiplier_im = GetConstReal(CompParams, 'Current Multiplier im', Found)
          
          i_multiplier = i_multiplier_re + im * i_multiplier_im

          StrandedHomogenization = GetLogical(CompParams, 'Homogenization Model', Found)
          IF ( .NOT. Found ) StrandedHomogenization = .FALSE.

          IF ( StrandedHomogenization ) THEN 
!            nu_11 = GetReal(CompParams, 'nu 11', Found)
!            nuim_11 = GetReal(CompParams, 'nu 11 im', FoundIm)
!            IF ( .NOT. Found .AND. .NOT. FoundIm ) CALL Fatal (Caller,'Homogenization Model nu 11 not found!')
!            nu_22 = GetReal(CompParams, 'nu 22', Found)
!            nuim_22 = GetReal(CompParams, 'nu 22 im', FoundIm)
!            IF ( .NOT. Found .AND. .NOT. FoundIm ) CALL Fatal (Caller,'Homogenization Model nu 22 not found!')
            sigma_33 = GetReal(CompParams, 'sigma 33', Found)
            sigmaim_33 = GetReal(CompParams, 'sigma 33 im', FoundIm)
            IF ( .NOT. Found .AND. .NOT. FoundIm ) CALL Fatal (Caller,'Homogenization Model Sigma 33 not found!')
          END IF
 
        CASE ('massive')

          VvarId = GetInteger (CompParams, 'Circuit Voltage Variable Id', Found)
          IF (.NOT. Found) CALL Fatal (Caller, 'Circuit Voltage Variable Id not found!')

        CASE ('foil winding')
          CALL GetLocalSolution(alpha,'Alpha')

          VvarId = GetInteger (CompParams, 'Circuit Voltage Variable Id', Found)
          IF (.NOT. Found) CALL Fatal (Caller, 'Circuit Voltage Variable Id not found!')

          coilthickness = GetConstReal(CompParams, 'Coil Thickness', Found)
          IF (.NOT. Found) CALL Fatal(Caller,'Foil Winding: Coil Thickness not found!')
 
          nofturns = GetConstReal(CompParams, 'Number of Turns', Found)
          IF (.NOT. Found) CALL Fatal(Caller,'Foil Winding: Number of Turns not found!')
 
          VvarDofs = GetInteger (CompParams, 'Circuit Voltage Variable dofs', Found)
          IF (.NOT. Found) CALL Fatal (Caller, 'Circuit Voltage Variable dofs not found!')
          InPlaneProximity = GetLogical(CompParams, 'Foil In Plane Proximity', Found)
          IF (InPlaneProximity) THEN
             LaminateThickness = coilthickness/nofturns
             LaminateModelPowerCompute = .TRUE.
          END IF
 
        CASE DEFAULT
          CALL Fatal (Caller, 'Non existent Coil Type Chosen!')
        END SELECT
      END IF

      
      ! Integrate local stresses:
      ! -------------------------
      IntegStuff = GaussPoints( Element )
      STIFF  = 0.0_dp
      FORCE  = 0.0_dp

      CALL GetLocalSolution( POT, VarName )

      Material => GetMaterial()
      IF( JouleHeating ) THEN
        BodyId = GetBody() 
        Cond(1:n) = GetReal( Material, 'Electric Conductivity', Found, Element)
      END IF

      IF( LossEstimation ) THEN
        BodyId = GetBody() 
        LossCoeff = ListGetFun( Material,'Fourier Loss Coefficient',Freq,Found )
        EddyLoss = .FALSE.
        IF (.NOT. Found) EddyLoss = .TRUE.
      END IF

      BertottiCompute = .FALSE.
      BRTc1 = GetCReal( Material,'Extended Bertotti Coefficient 1',Found ) 
      IF ( Found ) THEN
        BertottiCompute = .TRUE.
        Freq = Omega / (2*PI)
        BertottiLoss = 0.0_dp
        BRTc2 = GetCReal( Material,'Extended Bertotti Coefficient 2',Found ) 
        IF (.NOT. Found) CALL Fatal (Caller,'Extended Bertotti activated, &
                    Extended Bertotti Coefficient 2 not found!')

        BRTc3 = GetCReal( Material,'Extended Bertotti Coefficient 3',Found ) 
        IF (.NOT. Found) CALL Fatal (Caller,'Extended Bertotti activated, &
                    Extended Bertotti Coefficient 3 not found!')

        BRTc4 = GetCReal( Material,'Extended Bertotti Coefficient 4',Found ) 
        IF (.NOT. Found) BRTc4 = 1.5_dp

        BRTc5 = GetCReal( Material,'Extended Bertotti Coefficient 5',Found ) 
        IF (.NOT. Found) BRTc5 = 1.5_dp
      END IF

      LossUDF = .FALSE.
      CoreLossUDF = GetReal( Material,'Core Loss User Function', LossUDF ) 
      
      IF (BodyVolumesCompute) THEN
        BodyId = GetBody()
      END IF

      IF (ComplexPowerCompute) THEN
        BodyId = GetBody()
        Material => GetMaterial()

        IF (StrandedHomogenization) CALL Fatal (Caller,'Calculate Complex Power for Stranded & 
                                                 Homogenization model is not implemented.')

        mu = GetReal( Material, 'Relative Permeability', Found)
        mu = mu * 4.d-7*PI
        IF ( .NOT. Found ) CALL Warn(Caller, 'Relative Permeability not found!')
      END IF

      IF (LorentzForceCompute) BodyId = GetBody()

      DO t=1,IntegStuff % n
        Found = ElementInfo( Element, Nodes, IntegStuff % u(t), &
            IntegStuff % v(t), IntegStuff % w(t), detJ, Basis, dBasisdx )
        
        Weight = IntegStuff % s(t) * detJ
        grads_coeff = -1._dp/GetCircuitModelDepth()
        IF( CSymmetry ) THEN
          x = SUM( Basis(1:n) * Nodes % x(1:n) )
          Weight = Weight * x
          grads_coeff = grads_coeff/x
        END IF

        IF ( .NOT. ConstantBulkMatrixInUse ) THEN
          DO p=1,nd
            DO q=1,nd
              STIFF(p,q) = STIFF(p,q) + Weight * Basis(q) * Basis(p)
            END DO
          END DO
        END IF

        ! magnetic flux density components
        ! curl in cylindrically symmetric case has different sign convention.
        IF( CSymmetry ) THEN
          BAtIp(1) = -SUM( POT(1,1:nd) * dBasisdx(1:nd,2) )
          BAtIp(2) = SUM( POT(1,1:nd) * dBasisdx(1:nd,1) ) &
              + SUM( POT(1,1:nd) * Basis(1:nd) ) / x
          IF(FluxDofs == 4) THEN
            BAtIp(3) = -SUM( POT(2,1:nd) * dBasisdx(1:nd,2) )
            BAtIp(4) = SUM( POT(2,1:nd) * dBasisdx(1:nd,1) ) &
                + SUM( POT(2,1:nd) * Basis(1:nd) ) / x
          END IF
        ELSE
          BAtIp(1) =  SUM( POT(1,1:nd) * dBasisdx(1:nd,2) )
          BAtIp(2) = -SUM( POT(1,1:nd) * dBasisdx(1:nd,1) )
          IF(FluxDofs == 4) THEN
            BAtIp(3) =  SUM( POT(2,1:nd) * dBasisdx(1:nd,2) )
            BAtIp(4) = -SUM( POT(2,1:nd) * dBasisdx(1:nd,1) )
          END IF
        END IF
  
        ! Joule heating fields
        IF( TotDofs > 4 ) THEN
          IF ( StrandedHomogenization ) THEN 
            ValAtIp = SUM(Basis(1:n) * sigma_33(1:n))
            ValAtIpim = SUM(Basis(1:n) * sigmaim_33(1:n))
          ELSE
            ValAtIp = SUM( Basis(1:n) * Cond(1:n) )
            ValAtIpim = 0._dp
          END IF
          CondAtIp = ValAtIp + im * ValAtIpim
                                                         
          IF (.NOT. StrandedCoil ) THEN
            PotAtIp(1) =   Omega * SUM(POT(2,1:nd) * Basis(1:nd))
            PotAtIp(2) = - Omega * SUM(POT(1,1:nd) * Basis(1:nd))
          ELSE
            PotAtIp = 0._dp
          END IF

          localV=0._dp
          SELECT CASE (CoilType)
          CASE ('stranded')
            imag_value = LagrangeVar % Values(IvarId) + im * LagrangeVar % Values(IvarId+1)
            IF (i_multiplier /= 0._dp) THEN
              PotAtIp(1) = PotAtIp(1)+REAL(i_multiplier * imag_value * N_j / CondAtIp)
              PotAtIp(2) = PotAtIp(2)+AIMAG(i_multiplier * imag_value * N_j / CondAtIp)
            ELSE
              PotAtIp(1) = PotAtIp(1)+REAL(imag_value * N_j / CondAtIp)
              PotAtIp(2) = PotAtIp(2)+AIMAG(imag_value * N_j / CondAtIp)
            END IF            
          CASE ('massive')
            localV(1) = localV(1) + LagrangeVar % Values(VvarId)
            localV(2) = localV(2) + LagrangeVar % Values(VvarId+1)
            PotAtIp(1) = PotAtIp(1)-grads_coeff*localV(1)
            PotAtIp(2) = PotAtIp(2)-grads_coeff*localV(2)
          CASE ('foil winding')
            localAlpha = coilthickness *SUM(alpha(1:nd) * Basis(1:nd)) 
            DO k = 1, VvarDofs-1
              Reindex = 2*k
              Imindex = Reindex+1
              localV(1) = localV(1) + LagrangeVar % Values(VvarId+Reindex) * localAlpha**(k-1)
              localV(2) = localV(2) + LagrangeVar % Values(VvarId+Imindex) * localAlpha**(k-1)
            END DO
            PotAtIp(1) = PotAtIp(1)-grads_coeff*localV(1)
            PotAtIp(2) = PotAtIp(2)-grads_coeff*localV(2)
          END SELECT

          BAtIp(5) = 0.5_dp * ( PotAtIp(1)**2 + PotAtIp(2)**2 )
          BAtIp(6) = REAL(CondAtIp * BAtIp(5))
          TotalHeating = TotalHeating + Weight * BAtIp(6)
          imag_value = CondAtIp * (PotAtIp(1) + im * PotAtIp(2))
          BAtIp(7) = REAL(imag_value)
          BAtIp(8) = AIMAG(imag_value)
          imag_value = CMPLX(BatIp(1), BatIp(3), KIND=dp)
          imag_value2 = CMPLX(BatIp(2), BatIp(4), KIND=dp)
          BMagnAtIP = SQRT(ABS(imag_value*imag_value) + ABS(imag_value2*imag_value2))
        END IF
        
        IF (LorentzForceCompute) THEN
          BodyId = GetBody()
          ! Let's compute the JxB for all the bodies and 
          ! then we sum from these for the components which are outputted.

          Bx = CMPLX(BatIp(1), BatIp(3), KIND=dp)
          By = CMPLX(BatIp(2), BatIp(4), KIND=dp)
          Jz = CMPLX(BatIp(7), BatIp(8), KIND=dp)

          IF(Jz/=0) THEN
            LorentzForceDensX = ModelDepth * Weight * By / Jz * ABS(Jz)*ABS(Jz)
            LorentzForceDensY = -ModelDepth * Weight * Bx / Jz * ABS(Jz)*ABS(Jz)
          ELSE
            LorentzForceDensX = 0
            LorentzForceDensY = 0
          END IF

          BodyLorentzForcesRe(1, BodyId) = BodyLorentzForcesRe(1, BodyId) + &
            REAL(LorentzForceDensX)
          BodyLorentzForcesRe(2, BodyId) = BodyLorentzForcesRe(2, BodyId) + &
            REAL(LorentzForceDensY) 
          BodyLorentzForcesIm(1, BodyId) = BodyLorentzForcesIm(1, BodyId) + &
            AIMAG(LorentzForceDensX)
          BodyLorentzForcesIm(2, BodyId) = BodyLorentzForcesIm(2, BodyId) + &
            AIMAG(LorentzForceDensY) 
        END IF

        IF (LaminateModelPowerCompute) THEN
          ! This assumes linear reluctivity, and real conductivity
          skindepth = sqrt(2._dp/(omega * REAL(CondAtIp) * mu0))
          Lambda = LaminateThickness/skindepth
          Fsk = 3/Lambda * (SINH(Lambda) - SIN(Lambda))/(COSH(Lambda)-COS(Lambda))
          ! This is in W/m**3
          LaminatePowerDensity = 1._dp/24._dp * REAL(CondAtIp) * &
                (LaminateThickness * Omega * BMagnAtIP)**2._dp * Fsk
          TotalHeating = TotalHeating + Weight * ModelDepth * LaminatePowerDensity
        END IF

        IF (BertottiCompute) THEN
          ! Compute Bertotti loss for core
          BertottiLoss = BRTc1*Freq*BMagnAtIP**2.+ BRTc2*(Freq*BMagnAtIP)**2.+BRTc3*Freq**BRTc4*BMagnAtIP**BRTc5
          TotalHeating = TotalHeating + BertottiLoss
          BAtIp(6) = BAtIp(6) + BertottiLoss ! unorthodox
        END IF

        IF (LossUDF) THEN
          CoreLossUDFatIp = SUM(Basis(1:n) * CoreLossUDF(1:n))
          TotalHeating = TotalHeating + CoreLossUDFatIp
          BAtIp(6) = BAtIp(6) + CoreLossUDFatIp! unorthodox
        END IF

        IF( LossEstimation ) THEN
          IF ( EddyLoss ) THEN
            BodyLoss(BodyId) = BodyLoss(BodyId) + ModelDepth * Weight * BAtIp(6)
            IF (LaminateModelPowerCompute) & 
            BodyLoss(BodyId) = BodyLoss(BodyId) + ModelDepth * Weight * LaminatePowerDensity
          ELSE
            DO i=1,2
              ValAtIP = SUM( BAtIP(2*i-1:2*i) ** 2 )
              Coeff = Weight * LossCoeff * ( Freq ** FreqPower ) * ( ValAtIp ** FieldPower )
              ComponentLoss(i) = ComponentLoss(i) + Coeff
              BodyLoss(BodyId) = BodyLoss(BodyId) + Coeff
            END DO
          END IF
        END IF

        IF (ComplexPowerCompute) THEN
          cmplx_power = 0._dp
          imag_value = CMPLX(BAtIp(7), BAtIp(8))

          MuAtIp = SUM( Basis(1:n) * mu(1:n) )

          IF ( ABS(CondAtIp) > TINY(Weight) ) THEN
            cmplx_power = cmplx_power + ModelDepth * Weight * ABS(imag_value)**2._dp / CondAtIp 
          END IF

          imag_value = CMPLX(BatIp(1), BatIp(3), KIND=dp)
          imag_value2 = CMPLX(BatIp(2), BatIp(4), KIND=dp)
          cmplx_power = cmplx_power + im * ModelDepth * Weight * Omega/MuAtIp * &
              (ABS(imag_value)**2._dp+ABS(imag_value2)**2._dp)

          IF (LaminateModelPowerCompute) cmplx_power = cmplx_power + ModelDepth * Weight * LaminatePowerDensity

          BodyComplexPower(1,BodyId)=BodyComplexPower(1,BodyId) +  REAL(cmplx_power)
          BodyComplexPower(2,BodyId)=BodyComplexPower(2,BodyId) + AIMAG(cmplx_power)
 
        END IF

        IF (BodyICompute) THEN
          BodyCurrent(1,BodyId) = BodyCurrent(1,BodyId) + Weight * BatIp(7)
          IF (Fluxdofs==4) THEN
            BodyCurrent(2,BodyId) = BodyCurrent(2,BodyId) + Weight * BatIp(8)
          END IF
        END IF

        IF (BodyVolumesCompute) BodyVolumes(BodyId) = BodyVolumes(BodyId) + Weight * ModelDepth
       
        IF (AverageBCompute) THEN
          IF (BodyAverageBCompute(BodyId)) THEN
             BodyAvBre(1,BodyId) = BodyAvBre(1,BodyId) + Weight * BAtIp(1)
             BodyAvBre(2,BodyId) = BodyAvBre(2,BodyId) + Weight * BAtIp(2)
             IF (Fluxdofs==4) THEN
               BodyAvBim(1,BodyId) = BodyAvBim(1,BodyId) + Weight * BAtIp(3)
               BodyAvBim(2,BodyId) = BodyAvBim(2,BodyId) + Weight * BAtIp(4)
             END IF
          END IF
        END IF

        DO i=1,Totdofs
          Coeff = Weight * BAtIp(i)
          FORCE(i,1:nd) = FORCE(i,1:nd) + Coeff * Basis(1:nd)
        END DO
      END DO

!------------------------------------------------------------------------------
!      Update global matrices from local matrices 
!------------------------------------------------------------------------------
      IF ( .NOT. ConstantBulkMatrixInUse ) THEN
        Solver % Matrix % Rhs => SaveRHS
        CALL DefaultUpdateEquations( STIFF, FORCE(1,1:nd) )
      END IF

      DO i=1,TotDofs
        Solver % Matrix % RHS => ForceVector(:,i)
        CALL DefaultUpdateForce( FORCE(i,1:nd) )
      END DO

    END DO

    ! Check the total heating and normalize it, if requested
    IF( JouleHeating ) THEN
      TotalHeating = 2*PI*ParallelReduction(TotalHeating) 

      WRITE(Message,'(A,ES15.4)') 'Joule Heating (W): ',TotalHeating
      CALL Info(Caller,Message)
      CALL ListAddConstReal( Model % Simulation, 'res: Joule heating',TotalHeating)
      
      DesiredHeating = ListGetConstReal( SolverParams, &
          'Desired Heating Power',Found)        
      IF( Found .AND. TotalHeating > 0.0_dp ) THEN
        HeatingCoeff = DesiredHeating / TotalHeating

        WRITE(Message,'(A,ES15.4)') 'Joule coefficient: ',HeatingCoeff
        CALL Info(Caller,Message)
        CALL ListAddConstReal( Model % Simulation, 'res: Joule coefficient',HeatingCoeff)
      
        ForceVector(:,5) = HeatingCoeff * ForceVector(:,5) 
        ForceVector(:,6) = HeatingCoeff * ForceVector(:,6) 
      END IF
    END IF

    ! Assembly of the face terms (note: this gives no contribution to RHS)
    !---------------------------------------------------------------------
    IF ( .NOT. ConstantBulkMatrixInUse ) THEN
      IF (GetLogical(GetSolverParams(),'Discontinuous Galerkin',Found)) THEN
        IF (GetLogical(GetSolverParams(),'Average Within Materials',Found)) THEN
          FORCE = 0.0d0
          CALL AddLocalFaceTerms( STIFF, FORCE(1,:) )
        END IF
      END IF
    END IF

    IF( LossEstimation ) THEN
      DO j=1,2
        ComponentLoss(j) = ParallelReduction(ComponentLoss(j)) 
      END DO      
      DO j=1,Model % NumberOfBodies
        BodyLoss(j) = ParallelReduction(BodyLoss(j)) 
      END DO
      
      TotalLoss = SUM( ComponentLoss )
      CALL ListAddConstReal( Model % Simulation,'res: fourier loss',TotalLoss )
    
      !---------------------------------------------------------------------------------
      ! Screen output for componentwise and bodywise losses 
      !--------------------------------------------------------------------------------
      WRITE( Message,'(A,ES12.3)') 'Loss for cos mode: ', ComponentLoss(1)
      CALL Info(Caller, Message, Level=6 )
      WRITE( Message,'(A,ES12.3)') 'Loss for sin mode: ', ComponentLoss(2)
      CALL Info(Caller, Message, Level=6 )
      WRITE( Message,'(A,ES12.3)') 'Total loss: ',TotalLoss
      CALL Info(Caller,Message, Level=5 )

      CALL Info(Caller,'Losses by bodies',Level=6)
      DO j=1,Model % NumberOfBodies
         IF( BodyLoss(j) < TINY( TotalLoss ) ) CYCLE
         WRITE( Message,'(A,I0,A,ES12.3)') 'Body ',j,' : ',BodyLoss(j)
         WRITE (bodyNumber, "(I0)") j
         CALL ListAddConstReal( Model % Simulation,'res: Loss in Body '//TRIM(bodyNumber)//':', BodyLoss(j) )
         CALL Info(Caller, Message, Level=6 )
      END DO

      DEALLOCATE( BodyLoss )
    END IF

    IF (LorentzForceCompute) THEN
       DO j=1,Model % NumberOfBodies
         DO i = 1, 2
           BodyLorentzForcesRe(i,j) = ParallelReduction(BodyLorentzForcesRe(i,j)) 
           BodyLorentzForcesIm(i,j) = ParallelReduction(BodyLorentzForcesIm(i,j)) 
           IF (ISNAN(BodyLorentzForcesRe(i, j))) THEN
             BodyLorentzForcesRe(i, j)=0._dp
           END IF  
           IF (ISNAN(BodyLorentzForcesIm(i, j))) THEN
             BodyLorentzForcesIm(i, j)=0._dp
           END IF  
         END DO

         WRITE( Message,'(A,I0,A,ES12.3)') 'Body ',j,' : ',BodyLorentzForcesRe(1, j)
         WRITE (bodyNumber, "(I0)") j
         CALL ListAddConstReal( Model % Simulation,'res: Lorentz Force 1 re in Body '&
              //TRIM(bodyNumber)//':', BodyLorentzForcesRe(1,j) )
         CALL Info('Lorentz Force 1 re', Message, Level=6 )
         WRITE( Message,'(A,I0,A,ES12.3)') 'Body ',j,' : ',BodyLorentzForcesRe(2, j)
         CALL ListAddConstReal( Model % Simulation,'res: Lorentz Force 2 re in Body '&
              //TRIM(bodyNumber)//':', BodyLorentzForcesRe(2,j) )
         CALL Info('Lorentz Force 2 re', Message, Level=6 )

         WRITE( Message,'(A,I0,A,ES12.3)') 'Body ',j,' : ',BodyLorentzForcesIm(1, j)
         CALL ListAddConstReal( Model % Simulation,'res: Lorentz Force 1 im in Body '&
              //TRIM(bodyNumber)//':', BodyLorentzForcesIm(1,j) )
         CALL Info('Lorentz Force 1 im', Message, Level=6 )
         WRITE( Message,'(A,I0,A,ES12.3)') 'Body ',j,' : ',BodyLorentzForcesIm(2, j)
         CALL ListAddConstReal( Model % Simulation,'res: Lorentz Force 2 im in Body '&
              //TRIM(bodyNumber)//':', BodyLorentzForcesIm(2,j) )
         CALL Info('Lorentz Force 2 im', Message, Level=6 )
       END DO

       DO j = 1, NofComponents
         BodyIds => GetComponentBodyIds(j) 

         IF (ASSOCIATED(BodyIds)) THEN
           DO i = 1, 2
             DO k = 1, SIZE(BodyIds)
               bid = BodyIds(k)
               ComponentLorenzForcesRe(i,j) = ComponentLorenzForcesRe(i,j) &
                 + BodyLorentzForcesRe(i,bid)
               ComponentLorenzForcesIm(i,j) = ComponentLorenzForcesIm(i,j) &
                 + BodyLorentzForcesIm(i,bid)
             END DO
           END DO
  
           CALL ListAddConstReal( Model % Simulation,'res: Lorentz Force 1 re & 
                 in Component '//i2s(j), ComponentLorenzForcesRe(1,j) )
                         
           CALL ListAddConstReal( Model % Simulation,'res: Lorentz Force 2 re & 
                 in Component '//i2s(j), ComponentLorenzForcesRe(2,j) )

           CALL ListAddConstReal( Model % Simulation,'res: Lorentz Force 1 im & 
                 in Component '//i2s(j), ComponentLorenzForcesIm(1,j) )
                         
           CALL ListAddConstReal( Model % Simulation,'res: Lorentz Force 2 im & 
                 in Component '//i2s(j), ComponentLorenzForcesIm(2,j) )

         END IF
       END DO

    END IF

    IF (ComplexPowerCompute) THEN
       DO j=1,Model % NumberOfBodies
         DO i = 1, 2
           BodyComplexPower(i,j) = ParallelReduction(BodyComplexPower(i,j)) 
         END DO
         WRITE( Message,'(A,I0,A,ES12.3)') 'Body ',j,' : ',BodyComplexPower(1, j)
         WRITE (bodyNumber, "(I0)") j
         CALL ListAddConstReal( Model % Simulation,'res: Power re in Body '&
              //TRIM(bodyNumber)//':', BodyComplexPower(1,j) )
         CALL Info('Compex Power re', Message, Level=6 )
         WRITE( Message,'(A,I0,A,ES12.3)') 'Body ',j,' : ',BodyComplexPower(2, j)
         WRITE (bodyNumber, "(I0)") j
         CALL ListAddConstReal( Model % Simulation,'res: Power im in Body '&
              //TRIM(bodyNumber)//':', BodyComplexPower(2,j) )
         CALL Info('Compex Power im', Message, Level=6 )
       END DO

       DO j = 1, NofComponents
         BodyIds => GetComponentHomogenizationBodyIds(j) ! this will fall back to GetComponentBodyIds()

         IF (ASSOCIATED(BodyIds)) THEN
           DO i = 1, 2
             DO k = 1, SIZE(BodyIds)
               bid = BodyIds(k)
               CirCompComplexPower(i,j) = CirCompComplexPower(i,j) + BodyComplexPower(i,bid)
             END DO
           END DO
  
           CALL ListAddConstReal( Model % Simulation,'res: Power re & 
                 in Component '//i2s(j), CirCompComplexPower(1,j) )
                         
           CALL ListAddConstReal( Model % Simulation,'res: Power im & 
                 in Component '//i2s(j), CirCompComplexPower(2,j) )
         END IF
       END DO
    END IF

    IF ( BodyVolumesCompute ) THEN
      DO j=1,Model % NumberOfBodies
        BodyVolumes(j) = ParallelReduction(BodyVolumes(j))
      END DO
    END IF

    IF ( CirCompVolumesCompute ) THEN
      DO j=1,NofComponents
         BodyIds => GetComponentHomogenizationBodyIds(j) ! this will fall back to GetComponentBodyIds()
        IF (ASSOCIATED(BodyIds)) THEN
           DO k = 1, SIZE(BodyIds)
             bid = BodyIds(k)
             CirCompVolumes(j) = CirCompVolumes(j) + BodyVolumes(bid)
           END DO
        END IF
      END DO
    END IF
 
    IF (BodyICompute) THEN
      DO j = 1, Model % NumberOfBodies
        BodyCurrent(1, j) = ParallelReduction(BodyCurrent(1, j)) 
        WRITE (bodyNumber, "(I0)") j
        CALL ListAddConstReal( Model % Simulation,'res: Body Current re in Body ' &
                             //TRIM(bodyNumber)//':', BodyCurrent(1,j) )
        WRITE (Message,'(A,I0,A,ES12.3)') 'Body ',j,' : ',BodyCurrent(1,j)
        CALL Info('Body Current re', Message, Level=6 )
 
        IF (FluxDofs==4) THEN
          BodyCurrent(2, j) = ParallelReduction(BodyCurrent(2, j)) 
          CALL ListAddConstReal( Model % Simulation,'res: Body Current im in Body ' &
                               //TRIM(bodyNumber)//':', BodyCurrent(2,j) )
          WRITE (Message,'(A,I0,A,ES12.3)') 'Body ',j,' : ',BodyCurrent(2,j)
          CALL Info('Body Current im', Message, Level=6 )
          END IF
      END DO

      DO j = 1, NofComponents
        BodyIds => GetComponentHomogenizationBodyIds(j) ! this will fall back to GetComponentBodyIds()
        IF (ASSOCIATED(BodyIds)) THEN
          DO i = 1, 2
            DO k = 1, SIZE(BodyIds)
              bid = BodyIds(k)
              CirCompCurrent(i,j) = CirCompCurrent(i,j) + BodyCurrent(i,bid)
            END DO
          END DO
        END IF
      END DO
 
    END IF
 
    IF (AverageBCompute) THEN
      DO j=1,Model % NumberOfBodies 
        IF (.NOT. BodyAverageBCompute(j)) CYCLE
        DO i=1,2
          BodyAvBre(i,j)=ParallelReduction(BodyAvBre(i,j))*ModelDepth/BodyVolumes(j) 
          WRITE (XYNumber, "(I0)") i 
          WRITE (bodyNumber, "(I0)") j
          CALL ListAddConstReal( Model % Simulation,'res: Average Magnetic Flux Density ' &
                               //TRIM(XYNumber)//' in Body ' &
                               //TRIM(bodyNumber)//':', BodyAvBre(i,j) )
          WRITE (Message,'(A,I0,A,ES12.3)') 'Body ',j,' : ',BodyAvBre(i,j)
          CALL Info('Average Magnetic Flux Density '//TRIM(XYNumber), Message, Level=6 )
          IF (Fluxdofs==4) THEN
            BodyAvBim(i,j)=ParallelReduction(BodyAvBim(i,j))*ModelDepth/BodyVolumes(j)
            WRITE (XYNumber, "(I0)") i 
            WRITE (bodyNumber, "(I0)") j
            CALL ListAddConstReal( Model % Simulation,'res: Average Magnetic Flux Density ' &
                                 //TRIM(XYNumber)//' im in Body ' &
                                 //TRIM(bodyNumber)//':', BodyAvBim(i,j) )
            WRITE (Message,'(A,I0,A,ES12.3)') 'Body ',j,' : ',BodyAvBim(i,j)
            CALL Info('Average Magnetic Flux Density '//TRIM(XYNumber)//' im', Message, Level=6 )
          END IF
        END DO
      END DO

      DO j = 1, NofComponents
        BodyIds => GetComponentHomogenizationBodyIds(j) ! this will fall back to GetComponentBodyIds()
        IF (ASSOCIATED(BodyIds)) THEN
          DO i = 1, 2
            DO k = 1, SIZE(BodyIds)
              bid = BodyIds(k)
              CirCompAvBre(i,j) = CirCompAvBre(i,j) & 
                  + BodyVolumes(bid)/CirCompVolumes(j) * BodyAvBre(i,bid)
              CirCompAvBim(i,j) = CirCompAvBim(i,j) &
                  + BodyVolumes(bid)/CirCompVolumes(j) * BodyAvBim(i,bid)
            END DO
          END DO
        END IF
      END DO

    END IF

    IF (HomogenizationParamCompute) THEN
      DO j = 1,Model % NumberOfBodies
        CALL ComputeHomogenizationParams(BodyCurrent(:,j), BodyAvBre(:,j), BodyAvBim(:,j), &
                                         BodyVolumes(j), BodyComplexPower(:,j), Omega, &
                                         BodySkinCond(:,j), BodyProxNu(:,j))
        WRITE (bodyNumber, "(I0)") j
      
        OutputComp = ListGetString(Model % Bodies(j) % Values, 'Homogenization Conductivity Output Component', Found)
        IF (Found) THEN
          CALL ListAddConstReal( Model % Simulation,'res: Homogenization Conductivity '&
                              //TRIM(OutputComp)//' re in Body '//TRIM(bodyNumber)//':', BodySkinCond(1,j) )
          WRITE (Message,'(A,I0,A,ES12.3)') 'Body ',j,' : ',BodySkinCond(1,j)
          CALL Info('Homogenization Conductivity '//TRIM(OutputComp)//' re', Message, Level=6 )

          CALL ListAddConstReal( Model % Simulation,'res: Homogenization Conductivity '&
                             //TRIM(OutputComp)//' im in Body '//TRIM(bodyNumber)//':', BodySkinCond(2,j) )
          WRITE (Message,'(A,I0,A,ES12.3)') 'Body ',j,' : ',BodySkinCond(2,j)
          CALL Info('Homogenization Conductivity '//TRIM(OutputComp)//' im', Message, Level=6 )
       END IF

       OutputComp = ListGetString(Model % Bodies(j) % Values, 'Homogenization Reluctivity Output Component', Found)
       IF (Found) THEN
          CALL ListAddConstReal( Model % Simulation,'res: Homogenization Reluctivity '&
                            //TRIM(OutputComp)//' re in Body '//TRIM(bodyNumber)//':', BodyProxNu(1,j) )
          WRITE (Message,'(A,I0,A,ES12.3)') 'Body ',j,' : ',BodyProxNu(1,j)
          CALL Info('Homogenization Reluctivity '//TRIM(OutputComp)//' re', Message, Level=6 )

          CALL ListAddConstReal( Model % Simulation,'res: Homogenization Reluctivity '&
                           //TRIM(OutputComp)//' im in Body '//TRIM(bodyNumber)//':', BodyProxNu(2,j) )
          WRITE (Message,'(A,I0,A,ES12.3)') 'Body ',j,' : ',BodyProxNu(2,j)
          CALL Info('Homogenization Reluctivity '//TRIM(OutputComp)//' im', Message, Level=6 )
        END IF
      END DO

      DO j = 1, NofComponents

        CALL ComputeHomogenizationParams(CirCompCurrent(:,j), CirCompAvBre(:,j), CirCompAvBim(:,j), &
                                         CirCompVolumes(j), CirCompComplexPower(:,j), Omega, &
                                         CirCompSkinCond(:,j), CirCompProxNu(:,j))

        WRITE (CompNumber, "(I0)") j
  
        OutputComp = ListGetString(Model % Components(j) % Values, 'Homogenization Conductivity Output Component', Found)
        IF (Found) THEN
          CALL ListAddConstReal( Model % Simulation,'res: sigma_'//TRIM(OutputComp)//'_component(' &
                      //TRIM(CompNumber)//') re ', CirCompSkinCond(1,j) )
          CALL ListAddConstReal( Model % Simulation,'res: sigma_'//TRIM(OutputComp)//'_component(' &
                      //TRIM(CompNumber)//') im ', CirCompSkinCond(2,j) )
        END IF
  
        OutputComp = ListGetString(Model % Components(j) % Values, 'Homogenization Reluctivity Output Component', Found)
        IF (Found) THEN
          CALL ListAddConstReal( Model % Simulation,'res: nu_'//TRIM(OutputComp)//'_component(' &
                      //TRIM(CompNumber)//') re ', CirCompProxNu(1,j) )
          CALL ListAddConstReal( Model % Simulation,'res: nu_'//TRIM(OutputComp)//'_component(' &
                      //TRIM(CompNumber)//') im ', CirCompProxNu(2,j) )
        END IF
      END DO
   END IF

   IF (BodyVolumesCompute)         DEALLOCATE(BodyVolumes)
   IF (CirCompVolumesCompute)      DEALLOCATE(CirCompVolumes)
   IF (AverageBCompute)            DEALLOCATE(BodyAvBre, BodyAvBim)
   IF (AverageBCompute)            DEALLOCATE(CirCompAvBre, CirCompAvBim)
   IF (BodyICompute)               DEALLOCATE(BodyCurrent)
   IF (BodyICompute)               DEALLOCATE(CirCompCurrent)
   IF (ComplexPowerCompute)        DEALLOCATE(BodyComplexPower)
   IF (ComplexPowerCompute)        DEALLOCATE(CirCompComplexPower)
   IF (HomogenizationParamCompute) DEALLOCATE(BodySkinCond     ,  &
       BodyProxNu       ,  & 
       CirCompSkinCond,  & 
       CirCompProxNu      )
   IF (LorentzForceCompute)        DEALLOCATE(BodyLorentzForcesRe, &
       BodyLorentzForcesIm, &
       ComponentLorenzForcesRe, &
       ComponentLorenzForcesIm)


    DEALLOCATE( POT, STIFF, FORCE, Basis, dBasisdx, mu, Cond, sigma_33, sigmaim_33, CoreLossUDF)

!------------------------------------------------------------------------------
  END SUBROUTINE BulkAssembly
!------------------------------------------------------------------------------

!-------------------------------------------------------------------
 SUBROUTINE ComputeHomogenizationParams(Current, AvBre, AvBim, Volume, ComplexPower, Omega, &
                                        SkinCond, ProxNu)
!-------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp) :: Current(2), AvBre(2), AvBim(2), Volume
    COMPLEX(KIND=dp) :: imag_value, imag_value2, Bav(2), I
    REAL(KIND=dp) :: SkinCond(2), ProxNu(2), ComplexPower(2)
    REAL(KIND=dp) :: Omega
    COMPLEX(KIND=dp), PARAMETER :: im=(0._dp,1._dp)

    IF (Current(1) > TINY(Omega) .OR. Current(2) > TINY(Omega)) THEN
      imag_value = CMPLX(ComplexPower(1), &
                         ComplexPower(2), &
                         KIND=dp)
      I = CMPLX(Current(1), Current(2))
      imag_value = imag_value*Volume/ABS(I)**2._dp
      imag_value2 = 1._dp/imag_value
      SkinCond(1) = REAL(imag_value2) 
      SkinCond(2) = AIMAG(imag_value2) 
    ELSE
      SkinCond(1) = TINY(Omega)
      SkinCond(2) = TINY(Omega)
    END IF

    IF ( AvBre(1) > TINY(Omega) .OR. AvBre(2) > TINY(Omega) .OR. &
         AvBim(1) > TINY(Omega) .OR. AvBim(2) > TINY(Omega)         ) THEN
      Bav(1) = CMPLX(AvBre(1), AvBim(1), KIND=dp)
      Bav(2) = CMPLX(AvBre(2), AvBim(2), KIND=dp)

      imag_value = CMPLX(ComplexPower(1), &
                         ComplexPower(2), &
                         KIND=dp)
      imag_value = imag_value / im / Volume / Omega / (ABS(Bav(1))**2._dp+ABS(Bav(2))**2._dp)

      ProxNu(1) = REAL(imag_value) 
      ProxNu(2) = AIMAG(-imag_value) 
    ELSE
      ProxNu(1) = HUGE(Omega)
      ProxNu(2) = HUGE(Omega)
    END IF

!-------------------------------------------------------------------
 END SUBROUTINE ComputeHomogenizationParams
!-------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AddLocalFaceTerms(STIFF,FORCE)
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: STIFF(:,:), FORCE(:)

     TYPE(Element_t),POINTER :: P1,P2,Face,Faces(:)
     INTEGER ::t,n,n1,n2,NumberOfFaces,dim

     dim = CoordinateSystemDimension()

     IF (dim==2) THEN
       Faces => Solver % Mesh % Edges
       NumberOfFaces = Solver % Mesh % NumberOfEdges
     ELSE
       Faces => Solver % Mesh % Faces
       NumberOfFaces = Solver % Mesh % NumberOfFaces
     END IF

     DO t=1,NumberOfFaces
       Face => Faces(t)
       IF ( .NOT. ActiveBoundaryElement(Face) ) CYCLE

       P1 => Face % BoundaryInfo % Left
       P2 => Face % BoundaryInfo % Right
       IF ( ASSOCIATED(P2) .AND. ASSOCIATED(P1) ) THEN
          IF(.NOT.ASSOCIATED(GetMaterial(P1),GetMaterial(P2))) CYCLE

          n  = GetElementNOFNodes(Face)
          n1 = GetElementNOFNodes(P1)
          n2 = GetElementNOFNodes(P2)

          CALL LocalJumps( STIFF,Face,n,P1,n1,P2,n2)
          CALL DefaultUpdateEquations( STIFF, FORCE, Face )
       END IF
     END DO
!------------------------------------------------------------------------------
  END SUBROUTINE AddLocalFaceTerms
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE LocalJumps( STIFF,Face,n,P1,n1,P2,n2)
!------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=dp) :: STIFF(:,:)
      INTEGER :: n,n1,n2
      TYPE(Element_t), POINTER :: Face, P1, P2
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: FaceBasis(n), P1Basis(n1), P2Basis(n2)
      REAL(KIND=dp) :: Jump(n1+n2), detJ, U, V, W, S
      LOGICAL :: Stat
      INTEGER :: i, j, p, q, t
      TYPE(GaussIntegrationPoints_t) :: IntegStuff

      TYPE(Nodes_t) :: FaceNodes, P1Nodes, P2Nodes
      SAVE FaceNodes, P1Nodes, P2Nodes
!------------------------------------------------------------------------------
      STIFF = 0._dp

      CALL GetElementNodes(FaceNodes, Face)
      CALL GetElementNodes(P1Nodes, P1)
      CALL GetElementNodes(P2Nodes, P2)
!------------------------------------------------------------------------------
!     Numerical integration over the edge
!------------------------------------------------------------------------------
      IntegStuff = GaussPoints( Face )

      DO t=1,IntegStuff % n
        U = IntegStuff % u(t)
        V = IntegStuff % v(t)
        W = IntegStuff % w(t)
        S = IntegStuff % s(t)

        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo(Face, FaceNodes, U, V, W, detJ, FaceBasis)

        S = S * detJ
        IF( CSymmetry ) THEN
          S = S * SUM( FaceNodes % x(1:n) * FaceBasis(1:n) ) 
        END IF


        ! Find basis functions for the parent elements:
        ! ---------------------------------------------
        CALL GetParentUVW(Face, n, P1, n1, U, V, W, FaceBasis)
        stat = ElementInfo(P1, P1Nodes, U, V, W, detJ, P1Basis)

        CALL GetParentUVW(Face, n, P2, n2, U, V, W, FaceBasis)
        stat = ElementInfo(P2, P2Nodes, U, V, W, detJ, P2Basis)

        ! Integrate jump terms:
        ! ---------------------
        Jump(1:n1) = P1Basis(1:n1)
        Jump(n1+1:n1+n2) = -P2Basis(1:n2)

        DO p=1,n1+n2
          DO q=1,n1+n2
            STIFF(p,q) = STIFF(p,q) + s * Jump(q)*Jump(p)
          END DO
        END DO
      END DO
!------------------------------------------------------------------------------
    END SUBROUTINE LocalJumps
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
END SUBROUTINE BSolver
!------------------------------------------------------------------------------

!> \}

