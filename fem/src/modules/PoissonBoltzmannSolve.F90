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
! *  Authors: Juha Ruokolainen, Peter R�back
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 08 Jun 1997
! *
! *****************************************************************************/
    
!------------------------------------------------------------------------------
!>   Solve the Poisson-Boltzmann equation for the electric potential. Basically this
!> a similar to electrostatic equation except for the free charges on the r.h.s. 
!> that are assumed to obey Boltzmann statistics.
!> \ingroup Solvers
!------------------------------------------------------------------------------
   SUBROUTINE PoissonBoltzmannSolve( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
     USE Types
     USE Lists
     USE Integration
     USE ElementDescription
     USE Differentials
     USE SolverUtils
     USE ElementUtils
     USE Adaptive
     USE DefUtils

     IMPLICIT NONE
!------------------------------------------------------------------------------
 
     TYPE(Model_t) :: Model
     TYPE(Solver_t), TARGET:: Solver
 
     REAL (KIND=DP) :: dt
     LOGICAL :: TransientSimulation
 
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Matrix_t), POINTER  :: StiffMatrix
     TYPE(Element_t), POINTER :: CurrentElement
     TYPE(Variable_t), POINTER :: Var
     TYPE(Nodes_t) :: ElementNodes

     REAL (KIND=DP), POINTER :: ForceVector(:), Potential(:)
     REAL (KIND=DP), POINTER :: ElectricField1(:), ElectricField2(:), ElectricField3(:)
     REAL (KIND=DP), POINTER :: ChargeField(:), EnergyField(:)
     REAL (KIND=DP), POINTER :: Field(:)
     REAL (KIND=DP), ALLOCATABLE ::  Permittivity(:), &
       LocalStiffMatrix(:,:), Load(:), LocalForce(:), LocalPot(:)

#ifdef USE_ISO_C_BINDINGS
     REAL (KIND=DP) :: Norm, RelativeChange, TotEnergy, at0
     REAL (KIND=DP) :: at, st
#else
     REAL (KIND=DP) :: Norm, RelativeChange, TotEnergy, at0, RealTime
     REAL (KIND=DP) :: at, st, CPUTime
#endif

     REAL (KIND=DP) :: Cboltz, Ccharge, Cunit, ReferenceTemperature, Npos, Nneg, &
         NewtonTol, PermittivityOfVacuum, &
         Betapos, Betaneg, Alphapos, Alphaneg
     INTEGER :: Zpos, Zneg, NonlinearIter, NewtonIter, iter, LocalNodes
     LOGICAL :: SymmetricCharges, NewtonLinearization = .FALSE., FluxBC
     TYPE(ValueList_t), POINTER :: SolverParams

     INTEGER, POINTER :: NodeIndexes(:)
     INTEGER, POINTER :: PotentialPerm(:), FieldPerm(:)
     INTEGER, POINTER :: ChargeFieldPerm(:)
     INTEGER :: i, j, k, n, t, istat, bf_id, DIM
 
     LOGICAL :: AllocationsDone = .FALSE., gotIt, Found, NonDimensional
     LOGICAL :: CalculateField, CalculateCharge, CalculateEnergy, ConstantWeights

     CHARACTER(LEN=MAX_NAME_LEN) :: EquationName, Name, NameField, NameEnergy, NameCharge


     SAVE LocalStiffMatrix, Load, LocalForce, LocalPot, &
          ElementNodes, CalculateCharge, CalculateEnergy, &
          AllocationsDone, ElectricField1, ElectricField2, ElectricField3, &
          ChargeField, Permittivity, NameField, NameEnergy, NameCharge, &
          CalculateField, ConstantWeights, NewtonLinearization

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
     Potential     => Solver % Variable % Values
     PotentialPerm => Solver % Variable % Perm
 
     LocalNodes = COUNT( PotentialPerm > 0 )
     IF ( LocalNodes <= 0 ) RETURN

     StiffMatrix => Solver % Matrix
     ForceVector => StiffMatrix % RHS

     Norm = Solver % Variable % Norm
     DIM = CoordinateSystemDimension()
     

     SolverParams => GetSolverParams()

     NonDimensional = .TRUE.

     NonlinearIter = GetInteger(   SolverParams, &
         'Nonlinear System Max Iterations', Found )
     IF ( .NOT.Found ) NonlinearIter = 1

     NewtonTol = GetConstReal( SolverParams, &
         'Nonlinear System Newton After Tolerance',  Found )

     NewtonIter = GetInteger(   SolverParams, &
         'Nonlinear System Newton After Iterations', Found )
     IF ( NewtonIter == 0) NewtonLinearization = .TRUE.

!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone ) THEN
       N = Model % MaxElementNodes
 
       ALLOCATE( ElementNodes % x(N),   &
                 ElementNodes % y(N),   &
                 ElementNodes % z(N),   &
                 Permittivity(N),       &
                 LocalForce(N),         &
                 LocalStiffMatrix(N,N), &
                 Load(N),               &
                 LocalPot(N),           &
                 STAT=istat )
 
       IF ( istat /= 0 ) THEN
         CALL Fatal( 'StatElecSolve', 'Memory allocation error 1' )
       END IF
 
       CalculateField = ListGetLogical( Solver % Values, &
           'Calculate Electric Field', GotIt )
       IF ( .NOT. GotIt )  CalculateField = .TRUE.
       IF ( CalculateField )  THEN
         ALLOCATE( ElectricField1( Model % NumberOfNodes ), &
             ElectricField2( Model % NumberOfNodes ), &
             ElectricField3( Model % NumberOfNodes ), &
             STAT=istat )
         IF ( istat /= 0 ) THEN
           CALL Fatal( 'StatElecSolve', 'Memory allocation error 2' )
         END IF
         WRITE (NameField,'(A,A)') TRIM(ComponentName(Solver % Variable)),' Field'
       END IF

       CalculateCharge = ListGetLogical( Solver % Values, &
           'Calculate Electric Charge', GotIt )
       IF ( .NOT. GotIt )  CalculateCharge = .TRUE.
       IF ( CalculateCharge )  THEN
         ALLOCATE( ChargeField( Model % NumberOfNodes ), STAT=istat )
         WRITE (NameCharge,'(A,A)') TRIM(ComponentName(Solver % Variable)),' Charge'
         IF ( istat /= 0 ) THEN
           CALL Fatal( 'StatElecSolve', 'Memory allocation error 3' )
         END IF
       END IF

       DO i = 1, Model % NumberOfEquations
         CalculateEnergy = ListGetLogical( Model % Equations(i) % Values, &
             'Calculate Electric Energy', GotIt )
         IF ( GotIt ) EXIT
       END DO
       IF ( .NOT. GotIt )  CalculateEnergy = ListGetLogical( Solver % Values, &
           'Calculate Electric Energy', GotIt )

       IF ( CalculateEnergy ) THEN
         ALLOCATE( EnergyField( Model%NumberOfNodes ), STAT=istat )
         WRITE (NameEnergy,'(A,A)') TRIM(ComponentName(Solver % Variable)),' Energy'
         IF ( istat /= 0 ) THEN
           CALL Fatal( 'StatElecSolve', 'Memory allocation error 4' )
         END IF
       END IF

       ConstantWeights = ListGetLogical( Solver % Values, &
           'Constant Weights', GotIt )

!------------------------------------------------------------------------------

       IF ( .NOT.ASSOCIATED( StiffMatrix % MassValues ) ) THEN
         ALLOCATE( StiffMatrix % Massvalues( Model % NumberOfNodes ) )
         StiffMatrix % MassValues = 0.0d0
       END IF

!------------------------------------------------------------------------------
!      Add electric field, charge density and electric energy to the variable list
!------------------------------------------------------------------------------
       IF(CalculateField) THEN         
         Field => ElectricField1
         CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, Solver, &
               TRIM(NameField)//' 1', 1, Field, PotentialPerm)           
         Field => ElectricField2
         CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, Solver, &
             TRIM(NameField)//' 2', 1, Field, PotentialPerm)         
         IF(DIM == 3) THEN
           Field => ElectricField3
           CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, Solver, &
               TRIM(NameField)//' 3', 1, Field, PotentialPerm)
         END IF
       END IF

       IF ( CalculateCharge ) THEN
          Field => ChargeField
          CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, &
               Solver, NameCharge, 1, Field, PotentialPerm)
       END IF
          
       IF ( CalculateEnergy ) THEN
          Field => EnergyField
          CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, &
               Solver, NameEnergy, 1, Field, PotentialPerm )
       END IF
   
       AllocationsDone = .TRUE.
     END IF

!---------------------------------------------------------------------------------
!    Update arrays for derived vars. Rather cumbersome due to possible adaptivity
!---------------------------------------------------------------------------------

     IF ( CalculateField ) THEN
       Var => VariableGet( Model % Variables, TRIM(NameField)//' 1' )
       ElectricField1 => Var % Values
       FieldPerm => Var % Perm
       ElectricField1 = 0.0d0
       
       Var => VariableGet( Model % Variables, TRIM(NameField)//' 2' )
       ElectricField2 => Var % Values
       ElectricField2 = 0.0d0
       
       IF ( DIM == 3 ) THEN
         Var => VariableGet( Model % Variables, TRIM(NameField)//' 3' )
         ElectricField3 => Var % Values
         ElectricField3 = 0.0d0
       END IF
     END IF

     IF ( CalculateCharge ) THEN
       Var => VariableGet( Model % Variables, TRIM(NameCharge) )
       ChargeField => Var % Values
       FieldPerm => Var % Perm
       ChargeField = 0.0d0
     END IF

     IF ( CalculateEnergy ) THEN
       Var => VariableGet( Model % Variables, TRIM(NameEnergy) )
       EnergyField => Var % Values
       FieldPerm => Var % Perm
       EnergyField = 0.0d0
     END IF

     PermittivityOfVacuum = ListGetConstReal( Model % Constants,'Permittivity Of Vacuum',GotIt)
     IF(.NOT. GotIt) PermittivityOfVacuum = 1.0d0
     Cboltz = ListGetConstReal( Model % Constants,'Boltzmann Constant')     
     Ccharge = ListGetConstReal( Model % Constants,'Unit Charge')
     Cunit = Ccharge / Cboltz 

!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------

     EquationName = ListGetString( Solver % Values, 'Equation' )
     

     DO iter=1,NonlinearIter

       at  = CPUTime()
       at0 = RealTime()
       CALL DefaultInitialize()
!------------------------------------------------------------------------------
       CALL Info( 'PoissonBoltzmannSolve', '-------------------------------------',Level=4 )
       WRITE( Message,* ) 'Potential iteration', iter
       CALL Info( 'PoissonBoltzmannSolve', Message, Level=4 )
       CALL Info( 'PoissonBoltzmannSolve', '-------------------------------------',Level=4 )
       CALL Info( 'PoissonBoltzmannSolve', 'Starting Assembly...', Level=4 )
       
       IF(iter > NewtonIter) NewtonLinearization = .TRUE.

!------------------------------------------------------------------------------
!    Do the assembly
!------------------------------------------------------------------------------
       DO t = 1, Solver % NumberOfActiveElements
         
         IF ( RealTime() - at0 > 1.0 ) THEN
           WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
               (Solver % Mesh % NumberOfBulkElements-t) / &
               (1.0*Solver % Mesh % NumberOfBulkElements)), ' % done'
           
           CALL Info( 'PoissonBoltzmannSolve', Message, Level=5 )
           
           at0 = RealTime()
         END IF
         
!------------------------------------------------------------------------------
!        Check if this element belongs to a body where potential
!        should be calculated
!------------------------------------------------------------------------------
         CurrentElement => Solver % Mesh % Elements(Solver % ActiveElements(t))
         Model % CurrentElement => CurrentElement
         
         NodeIndexes => CurrentElement % NodeIndexes         
         n = CurrentElement % TYPE % NumberOfNodes
         
         ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
         ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
         ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)

         LocalPot(1:n) = Potential(PotentialPerm(NodeIndexes(1:n)))
!------------------------------------------------------------------------------

         bf_id = ListGetInteger( Model % Bodies(CurrentElement % BodyId) % &
             Values, 'Body Force',gotIt, minv=1, maxv=Model % NumberOfBodyForces )
         
         IF ( gotIt ) THEN
           Load(1:n) = ListGetReal( Model % BodyForces(bf_id) % Values, &
               'Charge Density', n, NodeIndexes, GotIt )
         ELSE
           Load(1:n) = 0.0d0
         END IF
         
         CALL ElectrolyteMaterialParameters() 

!------------------------------------------------------------------------------
!      Get element local matrix, and rhs vector
!------------------------------------------------------------------------------
         CALL PoissonBoltzmannCompose( LocalStiffMatrix,LocalForce, &
             Permittivity,Load,CurrentElement,n,ElementNodes )
!------------------------------------------------------------------------------
!      Update global matrix and rhs vector from local matrix & vector
!------------------------------------------------------------------------------
         CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce )
!------------------------------------------------------------------------------
       END DO
       CALL DefaultFinishBulkAssembly()

!------------------------------------------------------------------------------
!     Neumann boundary conditions
!------------------------------------------------------------------------------
       DO t=Solver % Mesh % NumberOfBulkElements + 1, &
           Solver % Mesh % NumberOfBulkElements + &
           Solver % Mesh % NumberOfBoundaryElements
         
         CurrentElement => Solver % Mesh % Elements(t)
         Model % CurrentElement => CurrentElement

!------------------------------------------------------------------------------
         DO i=1,Model % NumberOfBCs
           IF ( CurrentElement % BoundaryInfo % Constraint == &
               Model % BCs(i) % Tag ) THEN

              FluxBC = ListGetLogical(Model % BCs(i) % Values,'Electric Flux BC',gotIt) 
              IF( GotIt .AND. .NOT. FluxBC) CYCLE

!------------------------------------------------------------------------------
!             Set the current element pointer in the model structure to
!             reflect the element being processed
!------------------------------------------------------------------------------
             Model % CurrentElement => Solver % Mesh % Elements(t)            
!------------------------------------------------------------------------------
             n = CurrentElement % TYPE % NumberOfNodes
             NodeIndexes => CurrentElement % NodeIndexes
             IF ( ANY( PotentialPerm(NodeIndexes) <= 0 ) ) CYCLE

!------------------------------------------------------------------------------
!             BC: epsilon@Phi/@n = sigma
!------------------------------------------------------------------------------
             Load(1:n) = ListGetReal( Model % BCs(i) % Values,'Surface Charge', &
                 n,NodeIndexes,gotIt )
             IF(.NOT. GotIt) CYCLE
            
             ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
             ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
             ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)

!------------------------------------------------------------------------------
!             Get element matrix and rhs due to boundary conditions ...
!------------------------------------------------------------------------------
             CALL StatElecBoundary( LocalStiffMatrix, LocalForce,  &
                 Load, CurrentElement, n, ElementNodes )
!------------------------------------------------------------------------------
!             Update global matrices from local matrices
!------------------------------------------------------------------------------
             CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce )
!------------------------------------------------------------------------------
           END IF ! of currentelement bc == bcs(i)
         END DO ! of i=1,model bcs
       END DO   ! Neumann BCs
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!    FinishAssembly must be called after all other assembly steps, but before
!    Dirichlet boundary settings. Actually no need to call it except for
!    transient simulations.
!------------------------------------------------------------------------------
       CALL DefaultFinishAssembly( )
!------------------------------------------------------------------------------
!    Dirichlet boundary conditions
!------------------------------------------------------------------------------

       CALL DefaultDirichletBCs( )

       at = CPUTime() - at
       WRITE( Message, * ) 'Assembly (s)          :',at
       CALL Info( 'PoissonBoltzmannSolve', Message, Level=4 )
!------------------------------------------------------------------------------
!    Solve the system and we are done.
!------------------------------------------------------------------------------
       st = CPUTime()
       
       Norm = DefaultSolve()
       RelativeChange = Solver % Variable % NonlinChange

       st = CPUTime() - st
       WRITE( Message, * ) 'Solve (s)             :',st
       CALL Info( 'PoissonBoltzmannSolve', Message, Level=4 )

       WRITE( Message, * ) 'Result Norm   : ',Norm
       CALL Info( 'PoissonBoltzmannSolve', Message, Level=4 )
       Solver % Variable % Norm = Norm

       WRITE( Message, * ) 'Relative Change : ',RelativeChange
       CALL Info( 'PoissonBoltzmannSolve', Message, Level=4 )
       
       IF ( RelativeChange < NewtonTol ) NewtonLinearization = .TRUE.

       IF ( Solver % Variable % NonlinConverged == 1 ) EXIT
     END DO


!------------------------------------------------------------------------------
!    Compute the electric field from the potential: E = -grad Phi
!    Compute the electric flux: D = epsilon (-grad Phi)
!    Compute the total electric energy: W_e,tot = Integral (E . D)dV
!------------------------------------------------------------------------------

     IF ( CalculateField .OR. CalculateCharge .OR. CalculateEnergy) THEN 
       CALL GeneralElectricFlux( Model, Potential, PotentialPerm )
     END IF

     IF ( CalculateEnergy ) THEN
       WRITE( Message, * ) 'Tot. Electric Energy  :', TotEnergy
       CALL Info( 'PoissonBoltzmannSolve', Message, Level=4 )
       CALL ListAddConstReal( Model % Simulation, &
           'RES: Electric Energy', TotEnergy )
     END IF


!------------------------------------------------------------------------------
 
CONTAINS

  SUBROUTINE ElectrolyteMaterialParameters() 
    
    k = ListGetInteger( Model % Bodies(CurrentElement % BodyId) % &
        Values, 'Material', minv=1, maxv=Model % NumberOfMaterials )
    

    Betapos = ListGetConstReal(Model % Materials(k) % Values, &
         'Poisson Boltzmann Beta',GotIt)         
    IF(.NOT. GotIt) NonDimensional = .FALSE.
 
    Alphapos = ListGetConstReal(Model % Materials(k) % Values, &
        'Poisson Boltzmann Alpha',GotIt)         
    IF(.NOT. GotIt) NonDimensional = .FALSE.

    SymmetricCharges = .TRUE.
    IF(.NOT. NonDimensional) THEN
      Permittivity(1:n) = ListGetReal( Model % Materials(k) % Values, &
          'Relative Permittivity',n, NodeIndexes )
    
      Zpos = ListGetInteger(Model % Materials(k) % Values, &
          'Charge Number',GotIt)
      IF(.NOT. GotIt) THEN
        Zpos = ListGetInteger(Model % Materials(k) % Values, &
            'Positive Charge Number')
        Zneg = ListGetInteger(Model % Materials(k) % Values, &
            'Negative Charge Number')
        Zneg = -ABS(Zneg)
        SymmetricCharges = .FALSE.
      ELSE
        Zneg = -Zpos
      END IF
      ReferenceTemperature = ListGetConstReal(Model % Materials(k) % Values, &
         'Reference Temperature')
      Betapos = Cunit * Zpos / ReferenceTemperature

      Npos = ListGetConstReal(Model % Materials(k) % Values, &
          'Ion Density',GotIt)
      IF(.NOT. GotIt) THEN
        Npos = ListGetConstReal(Model % Materials(k) % Values, &
            'Positive Ion Density')
        Nneg = ListGetConstReal(Model % Materials(k) % Values, &
            'Negative Ion Density')
        SymmetricCharges = .FALSE.
      ELSE
        Nneg = Npos
      END IF
      Alphapos = 2.0 * Npos * Zpos * Ccharge 
    
      IF(.NOT. SymmetricCharges) THEN
        Betapos = Cunit * Zpos / ReferenceTemperature
        Betaneg = Cunit * Zneg / ReferenceTemperature
        Alphapos = Npos * Zpos * Ccharge 
        Alphaneg = Nneg * Zneg * Ccharge 
      END IF
    END IF

  END SUBROUTINE ElectrolyteMaterialParameters



!------------------------------------------------------------------------------
  SUBROUTINE PoissonBoltzmannCompose( StiffMatrix,Force,Permittivity, &
      Load,Element,n,Nodes )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: StiffMatrix(:,:),Force(:),Load(:), Permittivity(:)
    INTEGER :: n
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
 
    REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3)
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3)
    REAL(KIND=dp) :: SqrtElementMetric,U,V,W,S,A,B,L,C,x,y,z
    REAL(KIND=dp) :: Csinh, Ccosh, Pot
    LOGICAL :: Stat    
    INTEGER :: i,p,q,t,DIM    
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
 
!------------------------------------------------------------------------------
    DIM = CoordinateSystemDimension()
    
    Force = 0.0d0
    StiffMatrix = 0.0d0

!------------------------------------------------------------------------------
!      Numerical integration
!------------------------------------------------------------------------------
    IntegStuff = GaussPoints( Element )
    
    DO t=1,IntegStuff % n
      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)
      S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element,Nodes,U,V,W,SqrtElementMetric, &
          Basis,dBasisdx )
!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
      IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
        x = SUM( ElementNodes % x(1:n)*Basis(1:n) )
        y = SUM( ElementNodes % y(1:n)*Basis(1:n) )
        z = SUM( ElementNodes % z(1:n)*Basis(1:n) )
      END IF
      
      CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )
      
      S = S * SqrtElementMetric * SqrtMetric
      
      IF(NonDimensional) THEN
        C = 1.0_dp
      ELSE
        C = PermittivityOfVacuum * SUM( Permittivity(1:n) * Basis(1:n) )
      END IF
      Pot = SUM( LocalPot(1:n) * Basis(1:n) )

!------------------------------------------------------------------------------
!        The rhs at integration point 
!------------------------------------------------------------------------------

      L = SUM(Load(1:n) * Basis(1:n))

      IF(SymmetricCharges) THEN
        Csinh = SINH(Betapos * Pot)
        Ccosh = COSH(Betapos * Pot) 
      
        L = L - Alphapos * Csinh       
        IF(NewtonLinearization) THEN
          L = L + Alphapos * Betapos * Ccosh * Pot 
        END IF
      ELSE
        L = L + Alphapos * EXP(-Betapos * Pot) + Alphaneg * EXP(-Betaneg * Pot)
        IF(NewtonLinearization) THEN
          L = L + Alphapos * Betapos *  EXP(-Betapos * Pot) * Pot + &
              Alphaneg * Betaneg *  EXP(-Betaneg * Pot) * Pot
        END IF       
      END IF
      
!------------------------------------------------------------------------------
!        The Poisson equation
!------------------------------------------------------------------------------

      DO p=1,n
        DO q=1,n
          
          A = C * SUM( dBasisdx(p,1:DIM) * dBasisdx(q,1:DIM))

          IF(NewtonLinearization) THEN
            IF(SymmetricCharges) THEN
              A = A + Alphapos * Betapos * Ccosh * Basis(p) * Basis(q)
            ELSE
              A = A + Basis(p) * Basis(q) * &
                  ( Alphapos * Betapos * EXP(-Betapos * Pot) + &
                  Alphaneg * Betaneg * EXP(-Betaneg * Pot) )
            END IF
          END IF
          
          StiffMatrix(p,q) = StiffMatrix(p,q) + S * A
        END DO
        Force(p) = Force(p) + S * L * Basis(p)
      END DO
!------------------------------------------------------------------------------
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE PoissonBoltzmannCompose
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE StatElecBoundary( BoundaryMatrix, BoundaryVector, &
      LoadVector, Element, n, Nodes )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: BoundaryMatrix(:,:), BoundaryVector(:), LoadVector(:)
    TYPE(Nodes_t)   :: Nodes
    TYPE(Element_t) :: Element
    INTEGER :: n
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n)
    REAL(KIND=dp) :: dBasisdx(n,3),SqrtElementMetric
    REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3)       
    REAL(KIND=dp) :: u,v,w,s,x,y,z
    REAL(KIND=dp) :: Force
    REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)       
    INTEGER :: t,q,N_Integ       
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff       
    LOGICAL :: stat
!------------------------------------------------------------------------------

    BoundaryVector = 0.0d0
    BoundaryMatrix = 0.0d0
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
    IntegStuff = GaussPoints( Element )
    U_Integ => IntegStuff % u
    V_Integ => IntegStuff % v
    W_Integ => IntegStuff % w
    S_Integ => IntegStuff % s
    N_Integ =  IntegStuff % n
    
!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------
    DO t=1,N_Integ
      u = U_Integ(t)
      v = V_Integ(t)
      w = W_Integ(t)
!------------------------------------------------------------------------------
!     Basis function values & derivates at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
          Basis,dBasisdx )

!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
      IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
        x = SUM( ElementNodes % x(1:n)*Basis(1:n) )
        y = SUM( ElementNodes % y(1:n)*Basis(1:n) )
        z = SUM( ElementNodes % z(1:n)*Basis(1:n) )
      END IF
      
      CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )
      
      s = S_Integ(t) * SqrtElementMetric * SqrtMetric
      
!------------------------------------------------------------------------------
      Force = SUM( LoadVector(1:n)*Basis ) 
      
      DO q=1,N
        BoundaryVector(q) = BoundaryVector(q) + s * Basis(q) * Force
      END DO
    END DO
  END SUBROUTINE StatElecBoundary
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
! Compute the Electric Flux, Electric Field and Electric Energy at model nodes
!------------------------------------------------------------------------------
     SUBROUTINE GeneralElectricFlux( Model, Potential, Reorder )
!------------------------------------------------------------------------------
       TYPE(Model_t) :: Model
       REAL(KIND=dp) :: Potential(:)
       INTEGER :: Reorder(:)
!------------------------------------------------------------------------------
       TYPE(Element_t), POINTER :: Element
       TYPE(Nodes_t) :: Nodes 
       TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
       
       REAL(KIND=dp), POINTER :: U_Integ(:), V_Integ(:), W_Integ(:), S_Integ(:)
       REAL(KIND=dp), ALLOCATABLE :: SumOfWeights(:)
       REAL(KIND=dp) :: PermittivityOfVacuum
       REAL(KIND=dp) :: Basis(Model % MaxElementNodes)
       REAL(KIND=dp) :: dBasisdx(Model % MaxElementNodes,3)
       REAL(KIND=DP) :: SqrtElementMetric
       REAL(KIND=dp) :: ElementPot(Model % MaxElementNodes)
       REAL(KIND=dp) :: EnergyDensity
       REAL(KIND=dp) :: Flux(3), Field(3), ElemVol
       REAL(KIND=dp) :: s, ug, vg, wg, Grad(3), EpsGrad(3)
       REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)
       REAL(KIND=dp) :: x, y, z, L
       INTEGER, POINTER :: NodeIndexes(:)
       INTEGER :: n, N_Integ, t, tg, i, j, k, DIM
       LOGICAL :: Stat
       
!------------------------------------------------------------------------------

       ALLOCATE( Nodes % x( Model % MaxElementNodes ) )
       ALLOCATE( Nodes % y( Model % MaxElementNodes ) )
       ALLOCATE( Nodes % z( Model % MaxElementNodes ) )
       ALLOCATE( SumOfWeights( Model % NumberOfNodes ) )
       
       SumOfWeights = 0.0d0
       
       TotEnergy = 0.0d0

       PermittivityOfVacuum = ListGetConstReal( Model % Constants, &
           'Permittivity Of Vacuum',gotIt )
       IF ( .NOT.gotIt ) PermittivityOfVacuum = 1
       
       DIM = CoordinateSystemDimension()

!------------------------------------------------------------------------------
!   Go through model elements, we will compute on average of elementwise
!   fluxes to nodes of the model
!------------------------------------------------------------------------------
       DO t=1,Model % NumberOfBulkElements
!------------------------------------------------------------------------------
!        Check if this element belongs to a body where electrostatics
!        should be calculated
!------------------------------------------------------------------------------
         Element => Model % Elements(t)
         CurrentElement => Element

         NodeIndexes => Element % NodeIndexes
         IF ( .NOT. CheckElementEquation( Model, Element, EquationName ) ) CYCLE

         n = Element % TYPE % NumberOfNodes
         
         IF ( ANY(Reorder(NodeIndexes) == 0) ) CYCLE
         
         ElementPot(1:n) = Potential( Reorder( NodeIndexes(1:n) ) )

         bf_id = ListGetInteger( Model % Bodies(CurrentElement % BodyId) % &
             Values, 'Body Force',gotIt, minv=1, maxv=Model % NumberOfBodyForces )

         IF ( gotIt ) THEN
           Load(1:n) = ListGetReal( Model % BodyForces(bf_id) % Values, &
               'Charge Density', n, NodeIndexes, GotIt )
         ELSE
           Load(1:n) = 0.0d0
         END IF

         CALL ElectrolyteMaterialParameters() 
        
         Nodes % x(1:n) = Model % Nodes % x( NodeIndexes )
         Nodes % y(1:n) = Model % Nodes % y( NodeIndexes )
         Nodes % z(1:n) = Model % Nodes % z( NodeIndexes )
         
!------------------------------------------------------------------------------
!    Gauss integration stuff
!------------------------------------------------------------------------------
         IntegStuff = GaussPoints( Element )
         U_Integ => IntegStuff % u
         V_Integ => IntegStuff % v
         W_Integ => IntegStuff % w
         S_Integ => IntegStuff % s
         N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------

         k = ListGetInteger( Model % Bodies( Element % BodyId ) % &
             Values, 'Material', minv=1, maxv=Model % NumberOfMaterials )
         
         IF(.NOT. NonDimensional) THEN
           Permittivity(1:n) = ListGetReal( Model % Materials(k) % Values, &
               'Relative Permittivity',n, NodeIndexes )
         END IF
         
         EnergyDensity = 0.0d0
         Flux = 0.0d0
         Field = 0.0d0
         ElemVol = 0.0d0

!------------------------------------------------------------------------------
! Loop over Gauss integration points
!------------------------------------------------------------------------------
         DO tg=1,N_Integ
           
           ug = U_Integ(tg)
           vg = V_Integ(tg)
           wg = W_Integ(tg)
        
!------------------------------------------------------------------------------
! Need SqrtElementMetric and Basis at the integration point
!------------------------------------------------------------------------------
           stat = ElementInfo( Element, Nodes,ug,vg,wg, &
               SqrtElementMetric,Basis,dBasisdx )
        
!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
           s = 1
           IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
             x = SUM( Nodes % x(1:n)*Basis(1:n) )
             y = SUM( Nodes % y(1:n)*Basis(1:n) )
             z = SUM( Nodes % z(1:n)*Basis(1:n) )
             s = 2*PI
           END IF
           
           CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )
           
           s = s * SqrtMetric * SqrtElementMetric * S_Integ(tg)

!------------------------------------------------------------------------------

           EpsGrad = 0.0d0
           DO j = 1, DIM
             Grad(j) = SUM( dBasisdx(1:n,j) * ElementPot(1:n) )
           END DO
           EpsGrad(1:DIM) = SUM( Permittivity(1:n) * Basis(1:n) ) * Grad(1:DIM)
           
           EnergyDensity = EnergyDensity + s * SUM(Grad(1:DIM) * EpsGrad(1:DIM))
           DO j = 1,DIM
             Field(j) = Field(j) - Grad(j) * s
           END DO
           
           ElemVol = ElemVol + s
         END DO

         IF(CalculateEnergy) THEN
           TotEnergy = TotEnergy + EnergyDensity 
         END IF
          
!------------------------------------------------------------------------------
!   Weight with element area if required
!------------------------------------------------------------------------------

         IF ( ConstantWeights ) THEN
           EnergyDensity = EnergyDensity / ElemVol
           Field(1:DIM) = Field(1:DIM) / ElemVol
           SumOfWeights( Reorder( NodeIndexes(1:n) ) ) = &
               SumOfWeights( Reorder( NodeIndexes(1:n) ) ) + 1
         ELSE
           SumOfWeights( Reorder( NodeIndexes(1:n) ) ) = &
               SumOfWeights( Reorder( NodeIndexes(1:n) ) ) + ElemVol
         END IF

!------------------------------------------------------------------------------

         IF(CalculateEnergy) THEN
           EnergyField( FieldPerm(NodeIndexes(1:n)) ) = &
               EnergyField( FieldPerm(NodeIndexes(1:n)) ) + EnergyDensity
         END IF
         
         IF(CalculateField) THEN
           ElectricField1( FieldPerm( NodeIndexes(1:n) ) ) = &
               ElectricField1( FieldPerm( NodeIndexes(1:n) ) ) + Field(1)
           ElectricField2( FieldPerm( NodeIndexes(1:n) ) ) = &
               ElectricField2( FieldPerm( NodeIndexes(1:n) ) ) + Field(2)
           IF ( DIM == 3 )  THEN
             ElectricField3( FieldPerm( NodeIndexes(1:n) ) ) = &
                 ElectricField3( FieldPerm( NodeIndexes(1:n) ) ) + Field(3)
           END IF
         END IF
         
         IF(CalculateCharge) THEN
           ChargeField( FieldPerm (NodeIndexes(1:n)) ) = Load(1:n) 
           
           DO j=1,n
             IF(SymmetricCharges) THEN
               L = -Alphapos * SINH(Betapos * ElementPot(j))
             ELSE
               L = Alphapos * EXP(-Betapos * ElementPot(j)) + &
                   Alphaneg * EXP(-Betaneg * ElementPot(j))
             END IF
             ChargeField( FieldPerm(NodeIndexes(j))) = Load(j) + L
           END DO
         END IF

       END DO

!------------------------------------------------------------------------------
!   Finally, compute average of the fluxes at nodes
!------------------------------------------------------------------------------

       DO i = 1, Model % NumberOfNodes
         IF ( Reorder(i) == 0 )  CYCLE
         IF ( ABS( SumOfWeights(Reorder(i)) ) > AEPS ) THEN
           IF ( CalculateEnergy )  EnergyField(FieldPerm(i)) = &
               EnergyField(FieldPerm(i)) / SumOfWeights(Reorder(i))
           
           IF ( CalculateField ) THEN
             ElectricField1( FieldPerm(i) ) = ElectricField1( FieldPerm(i) ) / &
                 SumOfWeights( Reorder(i) )
             ElectricField2( FieldPerm(i) ) = ElectricField2( FieldPerm(i) ) / &
                 SumOfWeights( Reorder(i) )
             IF ( DIM == 3 ) THEN
               ElectricField3( FieldPerm(i) ) = ElectricField3( FieldPerm(i) ) / &
                   SumOfWeights( Reorder(i) )
             END IF
           END IF
           
         END IF
       END DO
       
       TotEnergy = PermittivityOfVacuum * TotEnergy / 2.0d0
       
       DEALLOCATE( Nodes % x, Nodes % y, Nodes % z, SumOfWeights)

!------------------------------------------------------------------------------
     END SUBROUTINE GeneralElectricFlux
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 END SUBROUTINE PoissonBoltzmannSolve
!------------------------------------------------------------------------------

