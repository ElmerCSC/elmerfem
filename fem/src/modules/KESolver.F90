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
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 08 Jun 1997
! *
! ****************************************************************************/

!------------------------------------------------------------------------------
!>  Subroutine for the k-e turbulence model.
!> \ingroup Solvers
!------------------------------------------------------------------------------
   SUBROUTINE KESolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
     USE DefUtils

     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Model_t)  :: Model
     TYPE(Solver_t) :: Solver
     REAL(KIND=dp) :: dt
     LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Matrix_t),POINTER  :: StiffMatrix
     INTEGER :: i,j,k,l,n,nd,nb,t,iter,k1,k2,body_id,eq_id,istat,LocalNodes,bf_id,DOFs

     TYPE(Nodes_t)   :: ElementNodes
     TYPE(Element_t),POINTER :: Element

     REAL(KIND=dp) :: RelativeChange,Norm,PrevNorm,S,C

     INTEGER, POINTER :: NodeIndexes(:)
     LOGICAL :: NewtonLinearization = .FALSE.,gotIt
!
     LOGICAL :: AllocationsDone = .FALSE., Bubbles

     CHARACTER(LEN=MAX_NAME_LEN) :: KEModel, V2FModel

     TYPE(Variable_t), POINTER :: FlowSol, KE

     INTEGER, POINTER :: FlowPerm(:),KinPerm(:)

     INTEGER :: NSDOFs,NewtonIter,NonlinearIter,NoActive
     REAL(KIND=dp) :: NewtonTol, Clip, V2FCp

     REAL(KIND=dp), POINTER :: KEpsilon(:),KineticDissipation(:), &
         FlowSolution(:), ElectricCurrent(:), ForceVector(:)

     REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), &
         STIFF(:,:),LayerThickness(:), &
         LOAD(:,:),FORCE(:),U(:),V(:),W(:), &
         Density(:),Viscosity(:),EffectiveVisc(:,:),Work(:),  &
         TurbulentViscosity(:),LocalDissipation(:), &
         LocalKinEnergy(:),KESigmaK(:),KESigmaE(:),KECmu(:),KEC1(:),&
         KEC2(:),C0(:,:), SurfaceRoughness(:), TimeForce(:),LocalV2(:),V2FCT(:)

     TYPE(ValueList_t), POINTER :: BC, Equation, Material

     SAVE U,V,W,MASS,STIFF,LOAD,FORCE, &
         ElementNodes,LayerThickness,Density,&
         AllocationsDone,Viscosity,LocalNodes,Work,TurbulentViscosity, &
         LocalDissipation,LocalKinEnergy,KESigmaK,KESigmaE,KECmu,C0, &
         SurfaceRoughness, TimeForce, KEC1, KEC2, EffectiveVisc, LocalV2, V2FCT

     REAL(KIND=dp), POINTER :: SecInv(:)
     SAVE SecInv

#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: at,at0,KMax, EMax, KVal, EVal
#else
     REAL(KIND=dp) :: at,at0,CPUTime,RealTime, KMax, EMax, KVal, EVal
#endif
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
     IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

     KE => Solver % Variable
     IF ( ASSOCIATED( KE ) ) THEN
       DOFs     =  KE % DOFs
       KinPerm  => KE % Perm
       KEpsilon => KE % Values
     END IF

     LocalNodes = COUNT( KinPerm > 0 )
     IF ( LocalNodes <= 0 ) RETURN

     FlowSol => VariableGet( Model % Variables, 'Flow Solution' )
     IF ( ASSOCIATED( FlowSol ) ) THEN
       FlowPerm     => FlowSol % Perm
       NSDOFs       =  FlowSol % DOFs
       FlowSolution => FlowSol % Values
     END IF

     StiffMatrix => Solver % Matrix
     ForceVector => StiffMatrix % RHS
     Norm = KE % Norm
!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone ) THEN
       N = Model % MaxElementNodes

       ALLOCATE( U( N ), V( N ), W( N ),  &
                 Density( N ),Work( N ),  &
                 Viscosity(N), &
                 EffectiveVisc(2,N), &
                 TurbulentViscosity(N), C0(DOFs,N), &
                 LayerThickness(N), &
                 SurfaceRoughness(N), &
                 KEC1(N), KEC2(N),      &
                 KESigmaK(N), KESigmaE(N),KECmu(N), &
                 LocalKinEnergy( N ),     &
                 LocalDissipation( N ),&
                 LocalV2(N), V2FCT(N), &
                 MASS( 2*DOFs*N,2*DOFs*N ), &
                 STIFF( 2*DOFs*N,2*DOFs*N ),LOAD( DOFs,N ), &
                 FORCE( 2*DOFs*N ), TimeForce( 2*DOFs*N ), STAT=istat )

       IF ( istat /= 0 ) THEN
         CALL Fatal( 'KESolver', 'Memory allocation error.' )
       END IF

       NULLIFY(SecInv)
       AllocationsDone = .TRUE.
     END IF

!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------
     NewtonTol = ListGetConstReal( Solver % Values, &
        'Nonlinear System Newton After Tolerance',gotIt )

     NewtonIter = ListGetInteger( Solver % Values, &
        'Nonlinear System Newton After Iterations',gotIt )

     NonlinearIter = ListGetInteger( Solver % Values, &
         'Nonlinear System Max Iterations',GotIt )

     IF ( .NOT.GotIt ) NonlinearIter = 1

     Bubbles = GetString(GetSolverParams(), &
               'Stabilization method', GotIt ) == 'bubbles'
     IF ( .NOT. GotIt ) Bubbles = .TRUE.

!------------------------------------------------------------------------------
      DO i=1,Model % NumberOFBCs
        BC => Model % BCs(i) % Values
        IF ( GetLogical(  BC, 'Noslip wall BC', gotit ) ) THEN
          CALL ListAddConstReal( BC, 'Kinetic Energy', 0.0_dp )
        END IF
      END DO
!------------------------------------------------------------------------------

     DO iter=1,NonlinearIter

       at  = CPUTime()
       at0 = RealTime()

       CALL Info( 'KESolver', ' ', Level=4 )
       CALL Info( 'KESolver', ' ', Level=4 )
       CALL Info( 'KESolver', &
          '-------------------------------------', Level=4 )
       WRITE( Message, * ) 'KEpsilon iteration: ', iter
       CALL Info( 'KESolver', Message, Level=4 )
       CALL Info( 'KESolver', &
          '-------------------------------------', Level=4 )
       CALL Info( 'KESolver', ' ', Level=4 )
       CALL Info( 'KESolver', 'Starting Assembly...', Level=4 )

       CALL DefaultInitialize()

!------------------------------------------------------------------------------
!      Bulk elements
!------------------------------------------------------------------------------
       body_id = -1
       CALL StartAdvanceOutput('KESolver','Assembly' )
       
       NoActive = GetNOFActive()
       
       DO t=1,NoActive

!------------------------------------------------------------------------------

         CALL AdvanceOutput(t,NoActive)

!------------------------------------------------------------------------------
!        Check if this element belongs to a body where kinetic energy
!        should be calculated
!------------------------------------------------------------------------------
         Element => GetActiveElement(t)
         IF ( Element % BodyId /= body_id ) THEN
            Material => GetMaterial()
            Equation => GetEquation()

            Clip = GetConstReal( Material, 'KE Clip', GotIt )
            IF ( .NOT.GotIt ) Clip = 1.0d-6

            KEModel = GetString( Material, 'KE Model', GotIt )
            IF ( .NOT. GotIt ) KEModel = 'standard'

            V2FModel = GetString( Material, 'V2-F Model', GotIt )
         END IF
!------------------------------------------------------------------------------
         CALL GetElementNodes( ElementNodes )
         n  = GetElementNOFNodes()
         nd = GetElementNOFDOFs()
         IF ( Bubbles ) nd=2*n
         nb = GetElementNOFBDOFs()
         NodeIndexes => Element % NodeIndexes
!------------------------------------------------------------------------------
         CALL GetScalarLocalSolution( LocalV2, 'V2' )
         CALL GetScalarLocalSolution( LocalKinEnergy, 'Kinetic Energy' )
         CALL GetScalarLocalSolution( LocalDissipation, 'Kinetic Dissipation' )

         CALL GetScalarLocalSolution( U, 'Velocity 1' )
         CALL GetScalarLocalSolution( V, 'Velocity 2' )
         CALL GetScalarLocalSolution( W, 'Velocity 3' )


         KESigmaK(1:n) = GetReal( Material, 'KE SigmaK', GotIt )
         IF ( .NOT. GotIt ) THEN
            KESigmaK = 1.0d0
            CALL ListAddConstReal( Material, 'KE SigmaK', KESigmaK(1) )
          END IF

         KESigmaE(1:n) = GetReal( Material, 'KE SigmaE', GotIt )
         IF ( .NOT. GotIt ) THEN
            SELECT CASE( KEModel )
            CASE( 'standard','v2-f' )
              KESigmaE = 1.3_dp
            CASE( 'rng' )
              KESigmaE = 1.0_dp
            CASE DEFAULT
               CALL Fatal( 'KESolver', 'Unknown K-Epsilon model' )
            END SELECT
            CALL ListAddConstReal( Material, 'KE SigmaE', KESigmaE(1) )
         END IF

         KECmu(1:n) = ListGetConstReal( Material, 'KE Cmu', GotIt )
         IF ( .NOT. GotIt ) THEN
            SELECT CASE( KEModel )
            CASE( 'standard' )
              KECmu = 0.09_dp
            CASE( 'v2-f')
              KECmu = 0.22_dp
            CASE( 'rng' )
              KECmu = 0.0845_dp
            CASE DEFAULT
               CALL Fatal( 'KESolver', 'Unknown K-Epsilon model' )
            END SELECT
            CALL ListAddConstReal( Material, 'KE Cmu', KECmu(1) )
         END IF

         KEC1(1:n) = GetReal( Material, 'KE C1', GotIt )
         IF ( .NOT. GotIt ) THEN
            SELECT CASE( KEModel )
            CASE( 'standard' )
              KEC1 = 1.44_dp
            CASE( 'v2-f' )
              KEC1(1:n) = 1.4_dp
            CASE( 'rng' )
              KEC1 = 1.42_dp
            CASE DEFAULT
               CALL Fatal( 'KESolver', 'Unknown K-Epsilon model' )
            END SELECT
            CALL ListAddConstReal( Material, 'KE C1', KEC1(1) )
         END IF

         KEC2(1:n) = GetReal( Material, 'KE C2', GotIt )
         IF ( .NOT. GotIt ) THEN
            SELECT CASE( KEModel )
            CASE( 'standard' )
              KEC2 = 1.92_dp
            CASE( 'v2-f' )
              KEC2 = 1.9_dp
            CASE( 'rng' )
              KEC2 = 1.68_dp
            CASE DEFAULT
               CALL Fatal( 'KESolver', 'Unknown K-Epsilon model' )
            END SELECT
            CALL ListAddConstReal( Material, 'KE C2', KEC2(1) )
         END IF

         IF ( KEModel == 'v2-f' ) THEN
           V2FCT(1:n) = GetReal( Material, 'V2-F CT', GotIt )
           IF ( .NOT. GotIt ) THEN
             V2FCT(1:n) = 6.0_dp
             CALL ListAddConstReal( Material, 'V2-F CT', V2FCT(1) )
           END IF

           V2FCp = GetCReal( Material, 'V2-F Cp', GotIt )
           IF ( .NOT. GotIt ) V2FCp = 0.05_dp
         END IF
!------------------------------------------------------------------------------
         Density(1:n)   = GetReal( Material,'Density' )
         Viscosity(1:n) = GetReal( Material,'Viscosity' )

!------------------------------------------------------------------------------
!        Get element local matrices, and RHS vectors
!------------------------------------------------------------------------------
         CALL LocalMatrix( MASS,STIFF,FORCE,LOAD, &
           U,V,W, Element,n,nd+nb,ElementNodes )
!------------------------------------------------------------------------------
         TimeForce = 0.0_dp
         IF ( TransientSimulation ) THEN
            CALL Default1stOrderTime( MASS, STIFF, FORCE )
         END IF
!------------------------------------------------------------------------------
!        Update global matrices from local matrices
!------------------------------------------------------------------------------
         IF ( Bubbles ) THEN
           CALL Condensate( DOFs*N, STIFF, FORCE, TimeForce )
         ELSE IF ( nb > 0 ) THEN
           CALL CondensateP( DOFs*nd, DOFs*nb, STIFF, FORCE, TimeForce )
         END IF
         CALL DefaultUpdateEquations( STIFF, FORCE )
!------------------------------------------------------------------------------
      END DO     !  Bulk elements
      CALL DefaultFinishBulkAssembly()      
      CALL Info( 'KESolver', 'Assembly done', Level=4 )

    
      IF(ListGetLogicalAnyBC(Model,'Epsilon Wall BC') .OR. &
          ListGetLogicalAnyBC(Model,'Noslip Wall BC') ) THEN
        DO t = 1, Solver % Mesh % NumberOfBoundaryElements
          Element => GetBoundaryElement(t)
          n = GetElementNOFNodes()

          BC => GetBC()
          IF ( ASSOCIATED( BC ) ) THEN
            IF ( GetLogical( BC, 'Epsilon Wall BC', gotIt ) .OR. &
                 GetLogical( BC, 'Noslip Wall BC',  gotIt ) ) THEN
              DO i=1,n
                j = KinPerm(Element % NodeIndexes(i))
                CALL ZeroRow( Solver % Matrix, 2*j )
                Solver % Matrix % RHS(2*j) = 0.0_dp
              END DO
             END IF
          END IF
        END DO
     END IF

!------------------------------------------------------------------------------
      DO t = 1, Solver % Mesh % NumberOfBoundaryElements
        Element => GetBoundaryElement(t)
        IF ( .NOT. ActiveBoundaryElement() ) CYCLE

!------------------------------------------------------------------------------
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes

        BC => GetBC()

        IF ( ASSOCIATED( BC ) ) THEN
          IF ( GetLogical( BC, 'Wall Law',gotIt ) ) THEN

            Density(1:n)   = GetParentMatProp( 'Density' )
            Viscosity(1:n) = GetParentMatProp( 'Viscosity' )

            SurfaceRoughness(1:n) = GetReal( BC, 'Surface Roughness', gotIt )
            LayerThickness(1:n)   = GetReal( BC, 'Boundary Layer Thickness' )

            DO j=1,n
              k = FlowPerm(NodeIndexes(j))
              IF ( k > 0 ) THEN
                SELECT CASE( NSDOFs )
                  CASE(3)
                    U(j) = FlowSolution( NSDOFs*k-2 )
                    V(j) = FlowSolution( NSDOFs*k-1 )
                    W(j) = 0.0D0

                  CASE(4)
                    U(j) = FlowSolution( NSDOFs*k-3 )
                    V(j) = FlowSolution( NSDOFs*k-2 )
                    W(j) = FlowSolution( NSDOFs*k-1 )
                END SELECT
              ELSE
                U(j) = 0.0d0
                V(j) = 0.0d0
                W(j) = 0.0d0
              END IF
            END DO

            DO j=1,n
              CALL KEWall( Work(1), Work(2), Work(3), SQRT(U(j)**2+V(j)**2+W(j)**2), &
               LayerThickness(j), SurfaceRoughness(j), Viscosity(j), &
                 Density(j) )

              k = DOFs*(KinPerm(NodeIndexes(j))-1)

              !ForceVector(k+1) = Work(1)
              !CALL ZeroRow( StiffMatrix,k+1 )
              !CALL SetMatrixElement( StiffMatrix,k+1,k+1,1.0d0 )

               CALL UpdateDirichletDof( StiffMatrix, k+1, Work(1) )
               CALL UpdateDirichletDof( StiffMatrix, k+2, Work(2) )

              !ForceVector(k+2) = Work(2)
              !CALL ZeroRow( StiffMatrix,k+2 )
              !CALL SetMatrixElement( StiffMatrix,k+2,k+2,1.0d0 )
            END DO
          END IF

          IF ( GetLogical( BC, 'Epsilon Wall BC', gotIt ) .OR. &
               GetLogical( BC, 'Noslip Wall BC',  gotIt ) ) THEN
            CALL EpsilonWall( Element, n, STIFF, FORCE )
          END IF
        END IF
      END DO
!------------------------------------------------------------------------------

      CALL DefaultFinishAssembly()
!------------------------------------------------------------------------------
!     Dirichlet boundary conditions
!------------------------------------------------------------------------------
      CALL DefaultDirichletBCs()
!------------------------------------------------------------------------------
      CALL Info( 'KESolver', 'Set boundaries done', Level=4 )
!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
      PrevNorm = Norm

      Norm = DefaultSolve()
!------------------------------------------------------------------------------
!      Kinetic Energy Solution should be positive
!------------------------------------------------------------------------------
      n = Solver % Mesh % NumberOfNodes
      Kmax = MAXVAL( Solver % Variable % Values(1:n:2) )
      Emax = MAXVAL( Solver % Variable % Values(2:n:2) )
      DO i=1,n
         k = Solver % Variable % Perm(i)
         IF ( k <= 0 ) CYCLE

         Kval = Solver % Variable % Values(2*k-1)
         Eval = Solver % Variable % Values(2*k-0)

         IF ( KVal < Clip*Kmax ) Kval = Clip*KMax

         IF ( Eval < Clip*EMax ) THEN
            KVal = Clip*EMax
            Eval = MAX(Density(1)*KECmu(1)*KVal**2/Viscosity(1),Clip*EMax)
         END IF

         Solver % Variable % Values(2*k-1) = MAX( KVal, 1.0d-10 )
         Solver % Variable % Values(2*k-0) = MAX( EVal, 1.0d-10 )
      END DO
!------------------------------------------------------------------------------
      RelativeChange = Solver % Variable % NonlinChange

      WRITE( Message,* ) 'Result Norm   : ',Norm
      CALL Info( 'KESolver', Message, Level = 4 )
      WRITE( Message,* ) 'Relative Change : ',RelativeChange
      CALL Info( 'KESolver', Message, Level = 4 )

      IF ( RelativeChange < NewtonTol .OR. &
              iter > NewtonIter ) NewtonLinearization = .TRUE.

      IF ( Solver % Variable % NonlinConverged == 1 ) EXIT
!------------------------------------------------------------------------------
    END DO
!------------------------------------------------------------------------------

    n = SIZE( Solver % Variable % Values )
    KE => VariableGet( Solver % Mesh % Variables, 'Kinetic Energy' )
    IF (ASSOCIATED(KE)) KE % Values = Solver % Variable % Values(1:n:2)

    KE => VariableGet( Solver % Mesh % Variables, 'Kinetic Dissipation' )
    IF (ASSOCIATED(KE)) KE % Values = Solver % Variable % Values(2:n:2)

CONTAINS

!------------------------------------------------------------------------------
   SUBROUTINE LocalMatrix( MASS,STIFF,FORCE, &
         LOAD, UX,UY,UZ,Element,n,nd,Nodes )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Return element local matrices and RSH vector for diffusion-convection
!  equation: 
!
!  ARGUMENTS:
!
!  REAL(KIND=dp) :: MASS(:,:)
!     OUTPUT: time derivative coefficient matrix
!
!  REAL(KIND=dp) :: STIFF(:,:)
!     OUTPUT: rest of the equation coefficients
!
!  REAL(KIND=dp) :: FORCE(:)
!     OUTPUT: RHS vector
!
!  REAL(KIND=dp) :: LOAD(:)
!     INPUT:
!
!  REAL(KIND=dp) :: UX(:),UY(:),UZ(:)
!     INPUT: Nodal values of velocity components from previous iteration
!           used only if coefficient of the convection term (C1) is nonzero
!
!  TYPE(Element_t) :: Element
!       INPUT: Structure describing the element (dimension,nof nodes,
!               interpolation degree, etc...)
!
!  INTEGER :: n
!       INPUT: Number of element nodes
!
!  TYPE(Nodes_t) :: Nodes
!       INPUT: Element node coordinates
!
!******************************************************************************
     USE MaterialModels

     IMPLICIT NONE

     REAL(KIND=dp), DIMENSION(:)   :: FORCE,UX,UY,UZ
     REAL(KIND=dp), DIMENSION(:,:) :: MASS,STIFF,LOAD

     INTEGER :: n,nd

     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t) :: Element

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
!
     REAL(KIND=dp) :: ddBasisddx(nd,3,3)
     REAL(KIND=dp) :: Basis(nd)
     REAL(KIND=dp) :: dBasisdx(nd,3),detJ

     REAL(KIND=dp) :: Velo(3),dVelodx(3,3)

     REAL(KIND=dp) :: A(2,2),M(2,2),ProdK,ProdE,div,ProdTensor(3,3)
     INTEGER :: i,j,c,p,q,t,dim,NBasis
     REAL(KIND=dp) :: LoadatIp(2),Cmu,Rho,mu,Tmu,Effmu(2),TimeScale,Re_T,LC1,LC2,LC3

     REAL(KIND=dp) :: s,u,v,w, K,E,Eta,Strain(3,3), GradAsymm(3,3), nu, Cv,&
             alpha,oldalpha,dalpha,err,ww,olderr,derr

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     REAL(KIND=dp) :: SpecificHeatRatio, Pressure(n), ReferencePressure, &
          Sound_speed_sq, Mach_number_sq
     REAL(KIND=dp) :: Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3), &
              SqrtMetric,Gravity(3),Pr_rho(n),Pr,rho_g, c3(n)
     REAL(KIND=dp) :: C0(2),C1,CT,C2(2),dC2dx(3),SecInv,X,Y,Z,LV2,LCT,SigmaK,SigmaE

     REAL(KIND=dp), POINTER :: gWork(:,:)

     LOGICAL :: stat,Convection, UseRNGModel

!------------------------------------------------------------------------------

     dim = CoordinateSystemDimension()

     FORCE = 0.0D0
     STIFF = 0.0D0
     MASS  = 0.0D0

     NBasis = nd
     IF ( Bubbles ) NBasis = 2*n

     UseRNGModel = KEModel == 'rng'

     gWork => ListGetConstRealArray( Model % Constants,'Gravity',GotIt)
     IF ( GotIt ) THEN
       Gravity = gWork(1:3,1)*gWork(4,1)
     ELSE
       Gravity    =  0.00_dp
       Gravity(2) = -9.81_dp
     END IF

     Pr_rho(1:n) = GetReal( Material, 'Turbulent Prandtl Number', stat )
     IF ( .NOT. stat ) Pr_rho(1:n) = 0.85_dp

     c3(1:n) = GetReal( Material, 'Dissipation buoyancy coefficient', stat )
     IF ( .NOT. stat ) c3(1:n) = 0.0_dp


     CALL ElementDensity( Density, n )
     SpecificheatRatio = GetCReal( Material, 'Specific Heat Ratio', stat )
     CALL getScalarLocalSolution( Pressure, 'Pressure' )
     ReferencePressure = GetCReal( Material, 'Reference Pressure', stat )

!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IF ( Bubbles ) THEN
        IntegStuff = GaussPoints( element, element % TYPE % GaussPoints2 )
     ELSE
        IntegStuff = GaussPoints( element )
     END IF

!------------------------------------------------------------------------------
!    Now we start integrating
!------------------------------------------------------------------------------
     DO t=1,IntegStuff % n
       u = IntegStuff % u(t)
       v = IntegStuff % v(t)
       w = IntegStuff % w(t)
!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,detJ, &
             Basis,dBasisdx,Bubbles=Bubbles )
!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
       s = detJ * IntegStuff % s(t)
       IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
         X = SUM( Nodes % x(1:n)*Basis(1:n) )
         Y = SUM( Nodes % y(1:n)*Basis(1:n) )
         Z = SUM( nodes % z(1:n)*Basis(1:n) )
         CALL CoordinateSystemInfo(Metric,SqrtMetric,Symb,dSymb,X,Y,Z)

         s = s * SqrtMetric
       END IF

!      Velocity from previous iteration at the integration point
!------------------------------------------------------------------------------
       Velo = 0.0D0
       Velo(1) = SUM( UX(1:n)*Basis(1:n) )
       Velo(2) = SUM( UY(1:n)*Basis(1:n) )
       Velo(3) = SUM( UZ(1:n)*Basis(1:n) )

       dVelodx = 0.0d0
       DO i=1,dim
         dVelodx(1,i) = SUM( UX(1:n)*dBasisdx(1:n,i) )
         dVelodx(2,i) = SUM( UY(1:n)*dBasisdx(1:n,i) )
         dVelodx(3,i) = SUM( UZ(1:n)*dBasisdx(1:n,i) )
       END DO

       IF ( CurrentCoordinateSystem() == Cartesian ) THEN
         Strain = 0.5d0 * ( dVelodx + TRANSPOSE(dVelodx) )
         Secinv = 2 * SUM( Strain * Strain )
       ELSE
         SecInv = SecondInvariant( Velo,dVelodx,Metric,Symb ) / 2
       END IF

!------------------------------------------------------------------------------
!      Coefficient of the convection and time derivative terms
!      at the integration point
!------------------------------------------------------------------------------
       K = SUM( LocalKinEnergy(1:n) * Basis(1:n) )
       E = SUM( LocalDissipation(1:n) * Basis(1:n) )
       Eta =  SQRT(SecInv) * K / E

       Pr = SUM(Pr_rho(1:n) * Basis(1:n) )
       rho_g = 0._dp
       DO i=1,dim
         rho_g = rho_g + SUM(Density(1:n) * dBasisdx(1:n,i)) * Gravity(i)
       END DO
       rho_g = 0._dp

       mu  = SUM( Viscosity(1:n) * Basis(1:n) )
       rho = SUM( Density(1:n) * Basis(1:n) )

       SigmaK = SUM( KESigmaK(1:n) * Basis(1:n) )
       SigmaE = SUM( KESigmaE(1:n) * Basis(1:n) )

       Cmu = SUM( KECMu(1:n) * Basis(1:n) )

       LC1 = SUM( KEC1(1:n) * Basis(1:n) )
       LC2 = SUM( KEC2(1:n) * Basis(1:n) )
       LC3 = SUM( c3(1:n) * Basis(1:n) )

       Sound_speed_sq = (SUM(Basis(1:n)*Pressure(1:n))+ReferencePressure) * &
                     SpecificHeatRatio / rho
       Mach_number_sq = 0._dp
       IF ( Sound_speed_sq > 0._dp ) Mach_number_sq = K/Sound_speed_sq

       IF ( KEModel=='v2-f' ) THEN
         LV2 = SUM( LocalV2(1:n) * Basis(1:n) )
         LCT = SUM( V2FCT(1:n) * Basis(1:n) )
         LC1 = 1.4_dp * (1+V2FCp*SQRT(K/LV2))

         Timescale = MAX(K/E,LCT*SQRT(mu/rho/E))
         Tmu = Rho * Cmu * LV2 * TimeScale

!        div = Strain(1,1) + Strain(2,2) + Strain(3,3)
!        ProdTensor = 2*Tmu*Strain
!        DO i=1,3
!          ProdTensor(i,i) = ProdTensor(i,i) - 2.0_dp*(Tmu*div+K)/3.0_dp
!        END DO
!        Prod = SUM(ProdTensor * dVelodx) / Rho
       ELSE
         Tmu  = Rho * Cmu*K**2  / E
       END IF
       ProdK = Tmu * (SecInv-rho_g/(rho*Pr)) / rho
       ProdE = Tmu * (SecInv-LC3*rho_g/(rho*Pr)) / rho

       Effmu(1) = mu + Tmu / SigmaK
       Effmu(2) = mu + Tmu / SigmaE

       C0(1) = Rho
       IF ( KEModel == 'v2-f' ) THEN
         C0(2) = Rho * LC2 / TimeScale
       ELSE
         C0(2) = Rho * LC2 * E / K
       END IF

! moved to ---> rhs
!      C0(2) = C0(2) + Cmu*Rho*Eta**3*(1-Eta/4.38d0) / &
!                 (1.0d0 + 0.012d0*Eta**3) * E / K

       C1 = Rho
       CT = Rho
!------------------------------------------------------------------------------
!      Coefficient of the diffusion term &
!------------------------------------------------------------------------------
       Alpha = 1.0d0

       IF ( UseRNGModel ) THEN
          ww = mu / Effmu(1)
          alpha = 1.3929d0
          oldalpha = 1

          olderr = ABS((oldalpha-1.3929d0)/(1.0d0-1.3929d0))**0.6321d0
          olderr = olderr * ABS((oldalpha+2.3929d0)/(1.0d0+2.3929d0))**0.3679d0
          olderr = olderr - ww

          DO i=1,100
             err = ABS((alpha-1.3929d0)/(1.0d0-1.3929d0))**0.6321d0
             err = err * ABS((alpha+2.3929d0)/(1.0d0+2.3929d0))**0.3679d0
             err = err - ww
             derr = olderr - err
             olderr = err
             dalpha = oldalpha - alpha
             oldalpha = alpha
             alpha = alpha - 0.5 * err * dalpha / derr
             IF ( ABS(err) < 1.0d-8 ) EXIT
          END DO

          IF ( ABS(err) > 1.0d-8 ) THEN
             PRINT*,'huh: ', alpha
             alpha = 1.3929d0
          END IF
       END IF

       C2(1) = Alpha * Effmu(1)
       C2(2) = Alpha * Effmu(2)
!------------------------------------------------------------------------------
!       Loop over basis functions of both unknowns and weights
!------------------------------------------------------------------------------
       DO p=1,NBasis
       DO q=1,NBasis
!------------------------------------------------------------------------------
!         The diffusive-convective equation without stabilization
!------------------------------------------------------------------------------
          M = 0.0d0
          A = 0.0d0

          M(1,1) = CT * Basis(q) * Basis(p)
          M(2,2) = CT * Basis(q) * Basis(p)

          A(1,2) = C0(1) * Basis(q) * Basis(p)
          A(2,2) = C0(2) * Basis(q) * Basis(p)
!------------------------------------------------------------------------------
!         The diffusion term
!------------------------------------------------------------------------------
          IF ( CurrentCoordinateSystem() == Cartesian ) THEN
             DO i=1,dim
               A(1,1) = A(1,1) + C2(1) * dBasisdx(q,i) * dBasisdx(p,i)
               A(2,2) = A(2,2) + C2(2) * dBasisdx(q,i) * dBasisdx(p,i)
             END DO
          ELSE
             DO i=1,dim
               DO j=1,dim
                  A(1,1) = A(1,1) + Metric(i,j) * C2(1) * &
                       dBasisdx(q,i) * dBasisdx(p,i)

                  A(2,2) = A(2,2) + Metric(i,j) * C2(2) * &
                       dBasisdx(q,i) * dBasisdx(p,i)
               END DO
             END DO
          END IF

!------------------------------------------------------------------------------
!           The convection term
!------------------------------------------------------------------------------
          DO i=1,dim
            A(1,1) = A(1,1) + C1 * Velo(i) * dBasisdx(q,i) * Basis(p)
            A(2,2) = A(2,2) + C1 * Velo(i) * dBasisdx(q,i) * Basis(p)
          END DO

          DO i=1,2
             DO j=1,2
               STIFF(2*(p-1)+i,2*(q-1)+j) = STIFF(2*(p-1)+i,2*(q-1)+j)+s*A(i,j)
               MASS(2*(p-1)+i,2*(q-1)+j)  = MASS(2*(p-1)+i,2*(q-1)+j) +s*M(i,j)
             END DO
          END DO
        END DO
        END DO

        ! Load at the integration point:
        !-------------------------------
        LoadAtIP(1) = rho*(ProdK-2*E*Mach_number_sq)
        IF ( KEModel=='v2-f' ) THEN
          LoadAtIP(2) = Rho*LC1*ProdE/TimeScale
        ELSE
          LoadAtIP(2) = Rho*LC1*ProdE*E/K
        END IF

        IF ( UseRNGModel ) &
            LoadatIP(2) = LoadatIP(2) - Cmu*Rho*Eta**3*(1-Eta/4.38d0) / &
                       (1.0d0 + 0.012d0*Eta**3) * E**2 / K
!------------------------------------------------------------------------------
        DO p=1,NBasis
          FORCE(2*(p-1)+1) = FORCE(2*(p-1)+1)+s*LoadAtIp(1)*Basis(p)
          FORCE(2*(p-1)+2) = FORCE(2*(p-1)+2)+s*LoadAtIp(2)*Basis(p)
        END DO
      END DO
!------------------------------------------------------------------------------
   END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Wall law for the k-e turbulence model.
!------------------------------------------------------------------------------
   SUBROUTINE EpsilonWall( Element, n, STIFF, FORCE )
!------------------------------------------------------------------------------
     TYPE(Element_t), POINTER :: Element
     INTEGER :: n
     REAL(KIND=dp) :: STIFF(:,:), FORCE(:)
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: BasisB(32), Basis(32), BasisK(32), dBasisdx(32,3), &
                 EVals(32), KVals(32), detJ

     INTEGER :: i,j,c,p,q,t,np, dim,N_Integ,NBasis

     REAL(KIND=dp) :: s,u,v,w,E,K,Kder,X,Y,Z,Normal(3),mu,rho,Relax

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     LOGICAL :: stat
     TYPE(Element_t), POINTER :: Parent

     TYPE(Nodes_t) :: Nodes, ParentNodes
     SAVE Nodes, ParentNodes

!------------------------------------------------------------------------------
     dim = CoordinateSystemDimension()

     FORCE = 0.0_dp
     STIFF = 0.0_dp

     ! Integration stuff:
     !-------------------
     IntegStuff = GaussPoints( Element )

     Relax = GetCReal( BC, 'Epsilon Relax', stat )
     IF (.NOT. stat ) Relax = 1

     Parent => Element % BoundaryInfo % Left
     IF ( .NOT. ASSOCIATED(Parent) ) &
       Parent => Element % BoundaryInfo % Right
     IF(.NOT.ASSOCIATED(Parent))RETURN

     np = GetElementNOFDOFs(Parent)
     CALL GetElementNodes( Nodes )
     CALL GetElementNodes( ParentNodes, Parent )

     Density(1:np)   = GetReal( GetMaterial(Parent), 'Density',UElement=Parent )
     Viscosity(1:np) = GetReal( GetMaterial(Parent), 'Viscosity',UElement=Parent )

     CALL GetScalarLocalSolution( KVals, 'Kinetic Energy', Parent )
     CALL GetScalarLocalSolution( EVals, 'Kinetic dissipation', Parent )

     DO t=1,IntegStuff % n
       u = IntegStuff % u(t)
       v = IntegStuff % v(t)
       w = IntegStuff % w(t)

       ! Basis function values & derivatives at the integration point:
       !--------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,detJ, BasisB )

       s = detJ * IntegStuff % s(t)
       IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
         x = SUM( Nodes % x(1:n)*Basis(1:n) )
         y = SUM( Nodes % y(1:n)*Basis(1:n) )
         z = SUM( Nodes % z(1:n)*Basis(1:n) )
         s = s *  CoordinateSqrtMetric(x,y,z)
       END IF

       Normal = NormalVector( Element, Nodes, U, V, .TRUE. )

       CALL GetParentUVW( Element,n,Parent,np,U,V,W,BasisB )
       stat = ElementInfo( Parent,ParentNodes,U,V,W,detJ, &
             Basis,dBasisdx )

       IF ( ABS(u)>0.999_dp ) THEN
         u = -u
       ELSE IF ( ABS(v)>0.999_dp ) THEN
         v = -v
       ELSE IF ( ABS(w)>0.999_dp ) THEN
         w = -w
       END IF
       stat = ElementInfo( Parent,ParentNodes,U,V,W,detJ,BasisK )

       rho = SUM( Basis(1:np) * Density(1:np) )
       mu  = SUM( Basis(1:np) * Viscosity(1:np) )

       E = SUM( Basis(1:np) *Evals(1:np) ) 
       K = SUM( BasisK(1:np)*Kvals(1:np) )
       KVals(1:np) = SQRT(KVals(1:np))
       Kder = 0.0_dp
       DO i=1,3
         Kder = Kder+SUM(dBasisdx(1:np,i)*Kvals(1:np))*Normal(i)
       END DO

       DO p=1,np
         DO q=1,np
           STIFF(2*p,2*q)   = STIFF(2*p,2*q) + s*Basis(q)*Basis(p)
           STIFF(2*p,2*q-1) = STIFF(2*p,2*q-1) - s*Relax*2*mu/rho*Kder**2*BasisK(q)/K*Basis(p)
         END DO
         FORCE(2*p) = FORCE(2*p) + s*(1-Relax)*E*Basis(p)
       END DO
     END DO

     IF ( TransientSimulation ) THEN
       MASS = 0.0_dp
       CALL Default1stOrderTime(MASS, STIFF, FORCE, UElement=Parent)
     END IF
     CALL DefaultUpdateEquations(STIFF, FORCE, UElement=Parent)
!------------------------------------------------------------------------------
   END SUBROUTINE EpsilonWall
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  END SUBROUTINE KESolver
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Initialization for the primary solver: KESolver
!> \ingroup Solvers
!------------------------------------------------------------------------------
   SUBROUTINE KESolver_Init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
     USE DefUtils

     IMPLICIT NONE
!------------------------------------------------------------------------------

     TYPE(Model_t)  :: Model
     TYPE(Solver_t) :: Solver

     REAL(KIND=dp) :: dt
     LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: SolverParams
     LOGICAL :: Found
     CHARACTER(LEN=MAX_NAME_LEN) :: str
!------------------------------------------------------------------------------
     SolverParams => GetSolverParams()

     str = GetString( SolverParams,'Variable', Found )
     IF ( .NOT. Found ) str = 'K-eps'
     IF ( INDEX( str, '[' ) <= 0 ) THEN
       CALL ListAddString( SolverParams, 'Variable', &
             TRIM(str) // '[Kinetic Energy:1 Kinetic Dissipation:1]' )
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE KESolver_Init
!------------------------------------------------------------------------------

