!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 10 Nov 1997
! *
! ****************************************************************************/

!------------------------------------------------------------------------------
!>  Solver for the K-omega turbulence model.
!> \ingroup Solvers
!------------------------------------------------------------------------------
   SUBROUTINE KOmega( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
     USE DefUtils
     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Model_t)  :: Model
     TYPE(Solver_t) :: Solver
     REAL(KIND=dp) :: dt
     LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Matrix_t),POINTER  :: StiffMatrix
     INTEGER :: i,j,k,n,nd,nb,iter,t,body_id,eq_id,istat,LocalNodes,bf_id,DOFs

     TYPE(Nodes_t)   :: ElementNodes
     TYPE(Element_t),POINTER :: Element

     REAL(KIND=dp) :: RelativeChange,Norm
     LOGICAL :: Stabilize = .TRUE.,NewtonLinearization = .FALSE.,gotIt

     LOGICAL :: AllocationsDone = .FALSE.
     LOGICAL :: Bubbles, BubblesDefault

     TYPE(Variable_t), POINTER :: FlowSol, KE

     INTEGER, POINTER :: KinPerm(:)

     INTEGER :: NewtonIter,NonlinearIter,NoActive
     REAL(KIND=dp) :: NewtonTol

     REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), &
       STIFF(:,:), LOAD(:,:),FORCE(:), LocalKinEnergy(:), TimeForce(:), &
       LocalDissipation(:), PrevKinEnergy(:), PrevDissipation(:), xl(:), xlprev(:)

     TYPE(ValueList_t), POINTER :: BC, Equation, Material
     REAL(KIND=dp) :: at,at0,KMax, EMax, KVal, EVal

     SAVE MASS,STIFF,LOAD,FORCE, ElementNodes,AllocationsDone,TimeForce, &
       LocalKinEnergy, LocalDissipation, PrevKinEnergy, PrevDissipation, xl, xlprev

     ! Per-element bubble history (current and previous timestep), needed to
     ! form a consistent BDF(1) time derivative for a condensed bubble: see
     ! CondensatePTransient in MatrixAssembly.F90 and IncompressibleNSVec's
     ! LCondensate, which this mirrors (also mirrored in KESolver.F90,
     ! V2FSolver.F90 and SSTKomega.F90). Indexed by Element % ElementIndex
     ! with stride bxStride = DOFs*Mesh % MaxBDOFs (not DOFs*nb of any one
     ! element, so blocks stay aligned on a mesh with mixed bubble counts).
     REAL(KIND=dp), ALLOCATABLE, SAVE :: bx(:), bxprev(:)
     INTEGER, SAVE :: bxStride = 0, BubbleTimestep = -1
     INTEGER :: boff

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
     IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

     KE => Solver % Variable
     IF ( ASSOCIATED( KE ) ) THEN
       DOFs     =  KE % DOFs
       KinPerm  => KE % Perm
     END IF

     LocalNodes = COUNT( KinPerm > 0 )
     IF ( LocalNodes <= 0 ) RETURN

     Norm = KE % Norm
!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone ) THEN
       N = Solver % Mesh % MaxElementDOFs

       ALLOCATE( MASS( 2*DOFs*N,2*DOFs*N ), &
                 STIFF( 2*DOFs*N,2*DOFs*N ),LOAD( DOFs,N ), &
                 FORCE( 2*DOFs*N ), TimeForce( 2*DOFs*N ), &
                 LocalKinEnergy(N), LocalDissipation(N), &
                 PrevKinEnergy(N), PrevDissipation(N), &
                 xl(DOFs*N), xlprev(DOFs*N), STAT=istat )

       IF ( istat /= 0 ) THEN
         CALL Fatal( 'KOmega', 'Memory allocation error.' )
       END IF

       AllocationsDone = .TRUE.
     END IF

     ! Per-element bubble history for the transient condensed-bubble case
     ! (TransientSimulation, bubbles present, "Bubbles in Global System =
     ! False"): allocate once, sized by the mesh's own worst-case bubble
     ! count (not this solver's nb, which can vary element to element) times
     ! the number of BULK elements. See the matching block and its rationale
     ! in KESolver.F90.
     IF ( TransientSimulation .AND. .NOT. Solver % GlobalBubbles .AND. &
          Solver % Mesh % MaxBDOFs > 0 .AND. .NOT. ALLOCATED(bx) ) THEN
       bxStride = DOFs * Solver % Mesh % MaxBDOFs
       ALLOCATE( bx( bxStride * Solver % Mesh % NumberOfBulkElements ), &
                 bxprev( bxStride * Solver % Mesh % NumberOfBulkElements ) )
       bx = 0.0_dp
       bxprev = 0.0_dp
     END IF

     ! A new timestep started: the bubble part left over from the last solve
     ! of the previous timestep becomes "previous" for this one. Must happen
     ! only once per timestep, not once per call -- this solver may be called
     ! several times per timestep by the outer (Steady State) coupled
     ! iteration, and only the first such call should shift the history.
     IF ( TransientSimulation .AND. ALLOCATED(bx) .AND. &
          GetTimestep() /= BubbleTimestep ) THEN
       bxprev = bx
       BubbleTimestep = GetTimestep()
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

     BubblesDefault = ListGetLogical( Solver % Values, 'Bubbles', GotIt )
     IF ( .NOT.GotIt ) BubblesDefault = .TRUE.

!------------------------------------------------------------------------------
      DO i=1,Model % NumberOFBCs
        BC => Model % BCs(i) % Values
        IF ( GetLogical( BC, 'Noslip wall BC', gotit ) ) THEN
          CALL ListAddConstReal( BC, 'Kinetic Energy', 0.0_dp )
        END IF
      END DO
!------------------------------------------------------------------------------

     DO iter=1,NonlinearIter

       at  = CPUTime()
       at0 = RealTime()

       CALL Info( 'KOmega', ' ', Level=4 )
       CALL Info( 'KOmega', ' ', Level=4 )
       CALL Info( 'KOmega', &
          '-------------------------------------', Level=4 )
       WRITE( Message, * ) 'Komega iteration: ', iter
       CALL Info( 'KOmega', Message, Level=4 )
       CALL Info( 'KOmega', &
          '-------------------------------------', Level=4 )
       CALL Info( 'KOmega', ' ', Level=4 )
       CALL Info( 'KOmega', 'Starting Assembly...', Level=4 )

       CALL DefaultInitialize()

!------------------------------------------------------------------------------
!      Bulk elements
!------------------------------------------------------------------------------
       body_id = -1
       CALL StartAdvanceOutput( 'Komega', 'Assembly:')
       NoActive = GetNOFActive()

       DO t=1,NoActive

         CALL AdvanceOutput(t,NoActive)
!------------------------------------------------------------------------------
!        Check if this element belongs to a body where kinetic energy
!        should be calculated
!------------------------------------------------------------------------------
         Element => GetActiveElement(t)
         Bubbles = BubblesDefault .AND. .NOT. ASSOCIATED( Element % PDefs )
         Material => GetMaterial()

         n = GetElementNOFNodes()
         nd = GetElementNOFDOFs()
         IF ( Bubbles ) nd = 2*n
         nb = GetElementNOFBDOFs()
         CALL GetElementNodes( ElementNodes )

         IF ( TransientSimulation .AND. .NOT. Solver % GlobalBubbles .AND. nb > 0 ) THEN
           CALL GetScalarLocalSolution( LocalKinEnergy, 'Kinetic energy' )
           CALL GetScalarLocalSolution( LocalDissipation, 'Kinetic Dissipation' )
           CALL GetScalarLocalSolution( PrevKinEnergy, 'Kinetic energy', tStep=-1 )
           CALL GetScalarLocalSolution( PrevDissipation, 'Kinetic Dissipation', tStep=-1 )
         END IF
!------------------------------------------------------------------------------
!        Get element local matrices, and RHS vectors
!------------------------------------------------------------------------------
         CALL LocalMatrix( MASS,STIFF,FORCE,LOAD,Element,n,nd+nb,ElementNodes )
         TimeForce = 0.0_dp
         IF ( Bubbles ) THEN
           IF ( TransientSimulation ) CALL Default1stOrderTime( MASS, STIFF, FORCE )
           CALL Condensate( DOFs*N, STIFF, FORCE, TimeForce )
         ELSE IF ( nb > 0 ) THEN
           IF ( TransientSimulation .AND. .NOT. Solver % GlobalBubbles ) THEN
             ! A condensed bubble's own value from the previous timestep is not
             ! in the global solution vector (it was eliminated from it), so
             ! Default1stOrderTime cannot form its time derivative -- it would
             ! silently treat that history as zero. CondensatePTransient forms
             ! M/dt and M*xprev/dt over the FULL bubble-augmented block instead,
             ! using this element's own recorded bubble history, before
             ! eliminating the bubble rows/columns. Calling Default1stOrderTime
             ! as well would add M/dt to the retained block a second time. See
             ! the matching comment in KESolver.F90.
             xl(1:2*nd-1:2)     = LocalKinEnergy(1:nd)
             xl(2:2*nd:2)       = LocalDissipation(1:nd)
             xlprev(1:2*nd-1:2) = PrevKinEnergy(1:nd)
             xlprev(2:2*nd:2)   = PrevDissipation(1:nd)

             boff = (Element % ElementIndex - 1) * bxStride
             CALL CondensatePTransient( nd, nb, DOFs, dt, MASS, STIFF, FORCE, &
                 xlprev(1:DOFs*nd), xl(1:DOFs*nd), &
                 bxprev(boff+1:boff+DOFs*nb), bx(boff+1:boff+DOFs*nb) )
           ELSE
             IF ( TransientSimulation ) CALL Default1stOrderTime( MASS, STIFF, FORCE )
             CALL CondensateP( DOFs*nd, DOFs*nb, STIFF, FORCE, TimeForce )
           END IF
         ELSE
           IF ( TransientSimulation ) CALL Default1stOrderTime( MASS, STIFF, FORCE )
         END IF
!------------------------------------------------------------------------------
!        Update global matrices from local matrices
!------------------------------------------------------------------------------
         CALL DefaultUpdateEquations( STIFF, FORCE )

!------------------------------------------------------------------------------
      END DO     !  Bulk elements
      CALL DefaultFinishBulkAssembly()
      CALL Info( 'KOmega', 'Assembly done', Level=4 )

!------------------------------------------------------------------------------
      CALL DefaultFinishAssembly()

!------------------------------------------------------------------------------
!     Dirichlet boundary conditions
!------------------------------------------------------------------------------
      DO t=1,Solver % Mesh % NumberOfBoundaryElements
        Element => GetBoundaryElement(t) 
        IF ( .NOT. ActiveBoundaryElement() ) CYCLE
        n = GetElementNOFNodes()
        BC => GetBC()
        IF ( .NOT. ASSOCIATED(BC) ) CYCLE
        IF (GetLogical(BC, 'Omega Wall BC', gotIt ) .OR. &
            GetLogical(BC, 'Noslip Wall BC',  gotIt)) CALL OmegaWall(Element,n)
      END DO

      CALL DefaultDirichletBCs()
!------------------------------------------------------------------------------
      CALL Info( 'KOmega', 'Set boundaries done', Level=4 )
!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------

      Norm = DefaultSolve()
!------------------------------------------------------------------------------
!      Kinetic Energy Solution should be positive
!------------------------------------------------------------------------------
      n = SIZE( Solver % Variable % Values)
      Kmax = MAXVAL( Solver % Variable % Values(1:n:2) )
      Emax = MAXVAL( Solver % Variable % Values(2:n:2) )
      DO i=1,SIZE(Solver % Variable % Perm)
         k = Solver % Variable % Perm(i)
         IF ( k <= 0 ) CYCLE
         Kval = Solver % Variable % Values(2*k-1)
         Eval = Solver % Variable % Values(2*k-0)
         Solver % Variable % Values(2*k-1) = MAX( KVal, 1.0d-12 )
         Solver % Variable % Values(2*k-0) = MAX( EVal, 1.0d-12 )
      END DO

!------------------------------------------------------------------------------
      WRITE( Message,* ) 'Result Norm   : ',Norm
      CALL Info( 'KOmega', Message, Level = 4 )

      RelativeChange = Solver % Variable % NonlinChange
      WRITE( Message,* ) 'Relative Change : ',RelativeChange
      CALL Info( 'KOmega', Message, Level = 4 )

      IF ( Solver % Variable % NonlinConverged == 1 ) EXIT
!------------------------------------------------------------------------------
    END DO
!------------------------------------------------------------------------------

CONTAINS

!------------------------------------------------------------------------------
   SUBROUTINE LocalMatrix( MASS,STIFF,FORCE, LOAD, Element,n,nd,Nodes )
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

     REAL(KIND=dp), DIMENSION(:)   :: FORCE
     REAL(KIND=dp), DIMENSION(:,:) :: MASS,STIFF,LOAD

     INTEGER :: n, nd

     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t) :: Element

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
!
     REAL(KIND=dp) :: ddBasisddx(nd,3,3)
     REAL(KIND=dp) :: Basis(nd)
     REAL(KIND=dp) :: dBasisdx(nd,3),detJ

     REAL(KIND=dp) :: UX(n), UY(n), UZ(n), Velo(3), dVelodx(3,3), Energy(n), &
                      Dissipation(n), Distance(n), Density(n), Viscosity(n)

     REAL(KIND=dp) :: A(2,2),M(2,2),Prod,div,ProdTensor(3,3)
     INTEGER :: i,j,c,p,q,t,dim,NBasis
     REAL(KIND=dp) :: LoadatIp(2),Cmu,Rho,mu,Tmu,Effmu(2)

     REAL(KIND=dp) :: s,u,v,w, K,Omega,Strain(3,3), Vorticity(3,3), dist

     REAL(KIND=dp) :: StrainMeasure,VorticityMeasure,X,Y,Z,SigmaK, &
             SigmaO,Beta,CD,F1,F2,F3,F4,rGamma, GradK(3), GradO(3)

     REAL(KIND=dp) :: Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3),SqrtMetric

     LOGICAL :: stat
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

!------------------------------------------------------------------------------

     dim = CoordinateSystemDimension()

     Viscosity(1:n) = GetReal( Material, 'Viscosity' )
     Density(1:n) = GetReal( Material, 'Density' )

     CALL GetScalarLocalSolution( UX, 'Velocity 1' )
     CALL GetScalarLocalSolution( UY, 'Velocity 2' )
     CALL GetScalarLocalSolution( UZ, 'Velocity 3' )

     CALL GetScalarLocalSolution( Energy, 'Kinetic energy' )
     CALL GetScalarLocalSolution( Dissipation, 'Kinetic Dissipation' )

     FORCE = 0.0D0
     STIFF = 0.0D0
     MASS  = 0.0D0

     NBasis = nd

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
       Velo = 0.0_dp
       Velo(1) = SUM( UX(1:n)*Basis(1:n) )
       Velo(2) = SUM( UY(1:n)*Basis(1:n) )
       Velo(3) = SUM( UZ(1:n)*Basis(1:n) )

       dVelodx = 0.0_dp
       DO i=1,dim
         dVelodx(1,i) = SUM( UX(1:n)*dBasisdx(1:n,i) )
         dVelodx(2,i) = SUM( UY(1:n)*dBasisdx(1:n,i) )
         dVelodx(3,i) = SUM( UZ(1:n)*dBasisdx(1:n,i) )
       END DO

       IF ( CurrentCoordinateSystem() == Cartesian ) THEN
         Strain  = 0.5_dp * (dVelodx + TRANSPOSE(dVelodx))
         StrainMeasure = MAX( SQRT(2 * SUM(Strain * Strain)), 1.0d-10 )

         Vorticity = 0.5_dp * (dVelodx - TRANSPOSE(dVelodx))
         VorticityMeasure = SQRT(2 * SUM(Vorticity * Vorticity))
       ELSE
         StrainMeasure = SQRT(SecondInvariant( Velo,dVelodx,Metric,Symb )/2)
       END IF

!------------------------------------------------------------------------------

       K = SUM( Energy(1:n) * Basis(1:n) )
       Omega = SUM( Dissipation(1:n) * Basis(1:n) )

       mu   = SUM( Viscosity(1:n) * Basis(1:n) )
       rho  = SUM( Density(1:n) * Basis(1:n) )

       Beta   = 0.075_dp
       SigmaK = 2.000_dp
       SigmaO = 2.000_dp
       rGamma = 5._dp/9._dp

       Tmu = rho*K/MAX(Omega,1.0d-10)
       Effmu(1) = mu + Tmu / SigmaK
       Effmu(2) = mu + Tmu / SigmaO

       Prod = 2*Tmu*SUM(Strain*dVelodx)
!------------------------------------------------------------------------------
!      Loop over basis functions of both unknowns and weights
!------------------------------------------------------------------------------
       DO p=1,NBasis
       DO q=1,NBasis
          M = 0.0d0
          A = 0.0d0

          M(1,1) = rho * Basis(q) * Basis(p)
          M(2,2) = rho * Basis(q) * Basis(p)

          A(1,1) = A(1,1) + rho * 0.09_dp * Omega * Basis(q) * Basis(p)
          A(2,2) = A(2,2) + rho * Beta * Omega * Basis(q) * Basis(p)
!------------------------------------------------------------------------------
!         The diffusion term
!------------------------------------------------------------------------------
          IF ( CurrentCoordinateSystem() == Cartesian ) THEN
             DO i=1,dim
               A(1,1) = A(1,1) + Effmu(1) * dBasisdx(q,i) * dBasisdx(p,i)
               A(2,2) = A(2,2) + Effmu(2) * dBasisdx(q,i) * dBasisdx(p,i)
             END DO
          ELSE
             DO i=1,dim
               DO j=1,dim
                  A(1,1) = A(1,1) + Metric(i,j) * Effmu(1) * &
                       dBasisdx(q,i) * dBasisdx(p,j)

                  A(2,2) = A(2,2) + Metric(i,j) * Effmu(2) * &
                       dBasisdx(q,i) * dBasisdx(p,j)
               END DO
             END DO
          END IF

!------------------------------------------------------------------------------
!           The convection term
!------------------------------------------------------------------------------
          DO i=1,dim
            A(1,1) = A(1,1) + rho * Velo(i) * dBasisdx(q,i) * Basis(p)
            A(2,2) = A(2,2) + rho * Velo(i) * dBasisdx(q,i) * Basis(p)
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
        LoadAtIP(1) = Prod
        LoadAtIP(2) = rGamma * Prod * Omega / K

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
!> Wall law for the k-Omega turbulence model.
!------------------------------------------------------------------------------
   SUBROUTINE OmegaWall( Element,n )
!------------------------------------------------------------------------------
     TYPE(Element_t), TARGET :: Element
     INTEGER :: n
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: omega_wall,dist,mu(32),rho(32),x0(n),y0(n),z0(n),x,y,z
     INTEGER :: i,j,np
     TYPE(Element_t), POINTER :: Parent
!------------------------------------------------------------------------------
     Parent => Element % BoundaryInfo % Left
     IF ( .NOT. ASSOCIATED(Parent) ) &
       Parent => Element % BoundaryInfo % Right
     IF ( .NOT. ASSOCIATED(Parent) ) RETURN

     np = GetElementNOFNodes(Parent)

     rho(1:np)= GetReal( GetMaterial(Parent), 'Density', UElement=Parent )
     mu(1:np) = GetReal( GetMaterial(Parent), 'Viscosity', UElement=Parent )

     x0(1:n) = Model % Nodes % x(Element % NodeIndexes)      
     y0(1:n) = Model % Nodes % y(Element % NodeIndexes)      
     z0(1:n) = Model % Nodes % z(Element % NodeIndexes)      

     omega_wall = 1.d10
     DO i=1,np
       j = Parent % NodeIndexes(i)
       IF ( ANY( j==Element % NodeIndexes(1:n) ) ) CYCLE

       x = Model % Nodes % x(j)
       y = Model % Nodes % y(j)
       z = Model % Nodes % z(j)

       dist = MINVAL( (x-x0(1:n))**2 + (y-y0(1:n))**2 + (z-z0(1:n))**2 )
       IF ( dist < AEPS ) CYCLE

!      omega_wall = 2*mu(i)/0.09_dp/rho(i)/dist
       omega_wall = 6*mu(i)/rho(i)/0.075_dp/dist

       j = 2*Solver % Variable % Perm(j)
       !Solver % Matrix % RHS(j) = omega_wall
       !CALL ZeroRow( Solver % Matrix, j )
       !CALL SetMatrixElement( Solver % Matrix, j,j, 1.0_dp )

       CALL UpdateDirichletDof( Solver % Matrix, j, omega_wall )
     END DO

!    DO i=1,n
!      j = 2*Solver % Variable % Perm(Element % NodeIndexes(i))
!      Solver % Matrix % RHS(j) = 10*omega_wall
!      CALL ZeroRow( Solver % Matrix, j )
!      CALL SetMatrixElement( Solver % Matrix, j,j, 1.0_dp )
!    END DO

!------------------------------------------------------------------------------
   END SUBROUTINE OmegaWall
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  END SUBROUTINE KOmega
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Initialization for the primary solver: KOmega
!> \ingroup Solvers
!------------------------------------------------------------------------------
   SUBROUTINE KOmega_Init( Model,Solver,dt,TransientSimulation )
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

     ! Everything below is specific to a p-element bubble ("Element =
     ! p:.. b:.."); the legacy "Stabilization Method = Bubbles" path (no
     ! "Element" override) doesn't go through GetElementNOFBDOFs'
     ! Solver % GlobalBubbles branch at all, so touching the list there would
     ! be a no-op at best and, via bandwidth optimization/mesh-level bubble
     ! DOF bookkeeping that DOES consult Solver % GlobalBubbles regardless of
     ! which bubble path a solver actually uses, a real (if tiny) unintended
     ! perturbation at worst -- see the matching comment and diffuser_v2f
     ! regression in KESolver_Init, which caught this the same way.
     str = ListGetString( SolverParams,'Element', Found )
     IF ( Found ) THEN
       IF ( INDEX( str, 'b:' ) > 0 ) THEN
         ! Left in the global system a bubble mode is a free per-element
         ! unknown driven by strongly nonlinear K/Omega reaction terms, with
         ! no neighboring element to diffuse against and no floor. Condense
         ! it out locally by default instead, like KESolver_Init,
         ! V2F_LDM_Init, SSTKOmega_Init and IncompressibleNSVec already do
         ! for their own bubbles; CondensatePTransient below makes that
         ! choice work for transient runs too. ListAddNew, so an explicit
         ! sif setting still wins.
         CALL ListAddNewLogical(SolverParams, 'Bubbles in Global System', .FALSE.)

         ! The recovery of a transient condensed bubble (see bx/bxprev and
         ! CondensatePTransient in KOmega, mirroring IncompressibleNSVec's
         ! own bx/bxprev, and the identical logic in KESolver_Init) needs at
         ! least TWO solves within one timestep. Only relevant where a
         ! bubble is actually condensed out.
         IF ( TransientSimulation .AND. &
              .NOT. ListGetLogical(SolverParams,'Bubbles in Global System',Found) ) THEN
           CALL ListAddNewInteger(SolverParams, 'Nonlinear System Min Iterations', 2)
           CALL ListAddNewInteger(SolverParams, 'Nonlinear System Max Iterations', 2)
         END IF
       END IF
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE KOmega_Init
!------------------------------------------------------------------------------
