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
! *  Original Date: 10 Nov 1997
! *
! ****************************************************************************/

!------------------------------------------------------------------------------
!> Solver for the (RC)SST K-omega-turbulence model.
!> \ingroup Solvers
!------------------------------------------------------------------------------
   SUBROUTINE SSTKOmega( Model,Solver,dt,TransientSimulation )
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
     INTEGER :: i,j,k,n,iter,t,body_id,eq_id,istat,LocalNodes,bf_id,DOFs

     TYPE(Nodes_t)   :: ElementNodes
     TYPE(Element_t),POINTER :: Element

     REAL(KIND=dp) :: RelativeChange,Norm,PrevNorm
     LOGICAL :: Stabilize = .TRUE.,NewtonLinearization = .FALSE.,gotIt

     LOGICAL :: AllocationsDone = .FALSE.

     TYPE(Variable_t), POINTER :: FlowSol, KE

     INTEGER, POINTER :: KinPerm(:)

     INTEGER :: NewtonIter,NonlinearIter
     REAL(KIND=dp) :: NewtonTol

     REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), &
       STIFF(:,:), LOAD(:,:),FORCE(:), LocalKinEnergy(:), TimeForce(:)

     TYPE(ValueList_t), POINTER :: BC, Equation, Material

     SAVE MASS,STIFF,LOAD,FORCE, ElementNodes,AllocationsDone,TimeForce

#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: at,at0, KMax, EMax, KVal, EVal
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
                 FORCE( 2*DOFs*N ), TimeForce( 2*DOFs*N ), STAT=istat )

       IF ( istat /= 0 ) THEN
         CALL Fatal( 'SSTKOmega', 'Memory allocation error.' )
       END IF

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

       CALL Info( 'SSTKOmega', ' ', Level=4 )
       CALL Info( 'SSTKOmega', ' ', Level=4 )
       CALL Info( 'SSTKOmega', &
          '-------------------------------------', Level=4 )
       WRITE( Message, * ) 'SSTKomega iteration: ', iter
       CALL Info( 'SSTKOmega', Message, Level=4 )
       CALL Info( 'SSTKOmega', &
          '-------------------------------------', Level=4 )
       CALL Info( 'SSTKOmega', ' ', Level=4 )
       CALL Info( 'SSTKOmega', 'Starting Assembly...', Level=4 )

       CALL DefaultInitialize()

!------------------------------------------------------------------------------
!      Bulk elements
!------------------------------------------------------------------------------
       body_id = -1
       DO t=1,Solver % NumberOfActiveElements

         IF ( RealTime() - at0 > 1.0 ) THEN
           WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
            (Solver % NumberOfActiveElements-t) / &
               (1.0*Solver % NumberOfActiveElements)), ' % done'

           CALL Info( 'SSTKOmega', Message, Level=5 )
           at0 =RealTime()
         END IF
!------------------------------------------------------------------------------
!        Check if this element belongs to a body where kinetic energy
!        should be calculated
!------------------------------------------------------------------------------
         Element => GetActiveElement(t)
         Material => GetMaterial()

         n = GetElementNOFNodes()
         CALL GetElementNodes( ElementNodes )
!------------------------------------------------------------------------------
!        Get element local matrices, and RHS vectors
!------------------------------------------------------------------------------
         CALL LocalMatrix( MASS,STIFF,FORCE,LOAD,Element,n,ElementNodes )
         TimeForce = 0.0_dp
         IF ( TransientSimulation ) THEN
            CALL Default1stOrderTime( MASS, STIFF, FORCE )
         END IF
         CALL Condensate( DOFs*N, STIFF, FORCE, TimeForce )
!------------------------------------------------------------------------------
!        Update global matrices from local matrices
!------------------------------------------------------------------------------
         CALL DefaultUpdateEquations( STIFF, FORCE )

!------------------------------------------------------------------------------
      END DO     !  Bulk elements
      CALL Info( 'SSTKOmega', 'Assembly done', Level=4 )

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

        IF ( GetLogical( BC, 'Wall Law',gotIt ) ) THEN
          CALL KomegaWallLaw(n)
        END IF
      END DO

      CALL DefaultDirichletBCs()
!------------------------------------------------------------------------------
      CALL Info( 'SSTKOmega', 'Set boundaries done', Level=4 )
!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
      PrevNorm = Norm
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
      CALL Info( 'SSTKOmega', Message, Level = 4 )

      RelativeChange = Solver % Variable % NonlinChange
      WRITE( Message,* ) 'Relative Change : ',RelativeChange
      CALL Info( 'SSTKOmega', Message, Level = 4 )

      IF ( Solver % Variable % NonlinConverged == 1 ) EXIT
!------------------------------------------------------------------------------
    END DO
!------------------------------------------------------------------------------

CONTAINS

!------------------------------------------------------------------------------
   SUBROUTINE KomegaWallLaw(n)
!------------------------------------------------------------------------------
     INTEGER :: n
     LOGICAL :: GotIt
     INTEGER :: i,j,k,DOF
     REAL(KIND=dp) :: Density(n),Viscosity(n),SurfaceRoughness(n),LayerThickness(n), &
             U(n), V(n), W(n), Kin, Eps, Omega

     Density(1:n)   = GetParentMatProp( 'Density' )
     Viscosity(1:n) = GetParentMatProp( 'Viscosity' )

     SurfaceRoughness(1:n) = GetReal( BC, 'Surface Roughness', gotIt )
     LayerThickness(1:n)   = GetReal( BC, 'Boundary Layer Thickness' )

     CALL GetScalarLocalSolution(U, 'Velocity 1')
     CALL GetScalarLocalSolution(V, 'Velocity 2')
     CALL GetScalarLocalSolution(W, 'Velocity 3')

     DOFs = Solver % Variable % DOFs
     DO j=1,n
      CALL KEWall( Kin, Eps, Omega, SQRT(U(j)**2+V(j)**2+W(j)**2), &
       LayerThickness(j), SurfaceRoughness(j), Viscosity(j), &
         Density(j) )

       k = DOFs*(Solver % Variable % Perm(Element % NodeIndexes(j))-1)
       !Solver % Matrix % RHS(k+1) = Kin
       !CALL ZeroRow( Solver % Matrix,k+1 )
       !CALL SetMatrixElement( Solver % Matrix,k+1,k+1,1.0d0 )

       CALL UpdateDirichletDof( Solver % Matrix, k+1, Kin )
       CALL UpdateDirichletDof( Solver % Matrix, k+2, Omega )
       
       !Solver % Matrix % RHS(k+2) = Omega
       !CALL ZeroRow( Solver % Matrix,k+2 )
       !CALL SetMatrixElement( Solver % Matrix,k+2,k+2,1.0d0 )
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE KomegaWallLaw
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE LocalMatrix( MASS,STIFF,FORCE, LOAD, Element,n,Nodes )
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
     USE MaterialModels

     IMPLICIT NONE

     REAL(KIND=dp), DIMENSION(:)   :: FORCE
     REAL(KIND=dp), DIMENSION(:,:) :: MASS,STIFF,LOAD

     INTEGER :: n

     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t) :: Element

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
!
     REAL(KIND=dp) :: ddBasisddx(2*n,3,3)
     REAL(KIND=dp) :: Basis(2*n)
     REAL(KIND=dp) :: dBasisdx(2*n,3),detJ

     REAL(KIND=dp) :: UX(n), UY(n), UZ(n), Velo(3), dVelodx(3,3), Energy(n), &
                      Dissipation(n), Distance(n), Density(n), Viscosity(n)

     REAL(KIND=dp) :: A(2,2),M(2,2),Prod,div,ProdK,ProdO,ProdTensor(3,3),Ident(3,3)
     INTEGER :: i,j,c,p,q,t,dim,NBasis
     REAL(KIND=dp) :: LoadatIp(2),Cmu,Rho,mu,Tmu,Effmu(2)

     REAL(KIND=dp) :: s,u,v,w, K,Omega,Strain(3,3), Vorticity(3,3), dist, &
            Mach_number_sq, Sound_speed_sq

     REAL(KIND=dp) :: StrainMeasure,VorticityMeasure,X,Y,Z,SigmaK,  &
             SigmaO,Beta,CD,F1,F2,F3,F4,rGamma, GradK(3), GradO(3), &
             Gravity(3), rho_g, Pr_rho(n), Pr, c3_omega(n), c3, Pressure(n), &
             ReferencePressure, SpecificHeatRatio

     REAL(KIND=dp), POINTER :: gWork(:,:)

     REAL(KIND=dp) :: Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3),SqrtMetric

     LOGICAL :: stat, Bubbles
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

!------------------------------------------------------------------------------

     dim = CoordinateSystemDimension()

     gWork => ListGetConstRealArray( Model % Constants,'Gravity',GotIt)
     IF ( GotIt ) THEN
       Gravity = gWork(1:3,1)*gWork(4,1)
     ELSE
       Gravity    =  0.00_dp
       Gravity(2) = -9.81_dp
     END IF

     Viscosity(1:n) = GetReal( Material, 'Viscosity' )
     CALL ElementDensity( Density, n )

     Pr_rho(1:n) = GetReal( Material, 'Turbulent Prandtl Number', stat )
     IF ( .NOT. stat ) Pr_rho(1:n) = 0.85_dp

     c3_omega(1:n) = GetReal( Material, 'Dissipation buoyancy coefficient', stat )
     IF ( .NOT. stat ) c3_omega(1:n) = 0.0_dp

     SpecificheatRatio = GetCReal( Material, 'Specific Heat Ratio', stat )
     CALL getScalarLocalSolution( Pressure, 'Pressure' )
     ReferencePressure = GetCReal( Material, 'Reference Pressure', stat )

     CALL GetScalarLocalSolution( UX, 'Velocity 1' )
     CALL GetScalarLocalSolution( UY, 'Velocity 2' )
     CALL GetScalarLocalSolution( UZ, 'Velocity 3' )

     CALL GetScalarLocalSolution( Energy, 'Kinetic energy' )
     CALL GetScalarLocalSolution( Distance, 'Wall Distance' )
     CALL GetScalarLocalSolution( Dissipation, 'Kinetic Dissipation' )

     FORCE = 0.0D0
     STIFF = 0.0D0
     MASS  = 0.0D0

     NBasis = 2*n
     Bubbles = .TRUE.

     Ident = 0._dp
     DO i=1,3
       Ident(i,i) = 1._dp
     END DO
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IF ( Bubbles ) THEN
        IntegStuff = GaussPoints( element, element % Type % GaussPoints2 )
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
         StrainMeasure = MAX( SQRT(2*SUM(Strain*Strain)), 1.0d-10 )

         Vorticity = 0.5_dp * (dVelodx - TRANSPOSE(dVelodx))
         VorticityMeasure = SQRT(2*SUM(Vorticity*Vorticity))
       ELSE
         StrainMeasure = SQRT(SecondInvariant( Velo,dVelodx,Metric,Symb )/2)
       END IF

!------------------------------------------------------------------------------

       K = SUM( Energy(1:n) * Basis(1:n) )
       Omega = SUM( Dissipation(1:n) * Basis(1:n) )

       DO i=1,dim
         GradK(i) = SUM( dBasisdx(1:n,i) * Energy(1:n) )
         GradO(i) = SUM( dBasisdx(1:n,i) * Dissipation(1:n) )
       END DO

       mu   = SUM( Viscosity(1:n) * Basis(1:n) )
       rho  = SUM( Density(1:n) * Basis(1:n) )
       rho_g = 0._dp
       DO i=1,dim
         rho_g = rho_g + SUM(Density(1:n) * dBasisdx(1:n,i)) * Gravity(i)
       END DO
       dist = SUM( Distance(1:n) * Basis(1:n) )

       Sound_speed_sq = (SUM(Basis(1:n)*Pressure(1:n))+ReferencePressure) * &
                     SpecificHeatRatio / rho

       Mach_number_sq = 0._dp
       IF ( Sound_speed_sq > 0._dp ) Mach_number_sq = K / Sound_speed_sq

       Pr = SUM( Basis(1:n) * Pr_rho(1:n) )
       c3 = SUM( Basis(1:n) * c3_omega(1:n) )

       CD = MAX(2*rho/1.168_dp/Omega*SUM(GradK(1:dim)*GradO(1:dim)),1.d-10)

       F1 = SQRT(K) / 0.09_dp / Omega / dist
       F1 = MAX( F1, 500 * mu / rho / Omega / dist**2 )
       F1 = MIN( F1, 4 * rho / 1.168_dp* K / CD / dist**2 )
       F1 = TANH( F1**4 )

       F2 = MAX( 2*SQRT(K)/0.09_dp/Omega/dist,500*mu/rho/omega/dist**2 )
       F2 = TANH(F2**2)

       ! don't use SST near rough walls; NOT IN USE
!      F3 = 1-TANH((150*mu/rho/omega/dist**2)**4)
       F3 = 1

       ! Rotation and curvature correction of Schweighofer & Hellsten; NOT IN USE
!      F4 = VorticityMeasure/StrainMeasure*(VorticityMeasure/StrainMeasure-1)
!      F4 = 1._dp / (1+3.6_dp*F4)
       F4 = 1

       Beta   = 0.075_dp*F1 + 0.0828_dp*(1-F1) 
       SigmaK = 1.176_dp*F1 + 1.0000_dp*(1-F1)
       SigmaO = 2.000_dp*F1 + 1.1680_dp*(1-F1)

!      rGamma = Beta/0.09_dp - 0.41_dp**2/SigmaO/SQRT(0.09_dp)
       rGamma = 5._dp/9._dp * F1 + 0.44_dp * (1-F1)

       Tmu = 0.31_dp*rho*K/MAX(0.31_dp*Omega,VorticityMeasure*F2*F3)
       Effmu(1) = mu + Tmu / SigmaK
       Effmu(2) = mu + Tmu / SigmaO

       ProdK = SUM((2*Tmu*Strain-2/3._dp*rho*Ident*K)*dVelodx) - &
             Tmu * rho_g / (rho * Pr) - 2*rho*0.09_dp*K*Omega*Mach_number_sq

       ProdO = SUM((2*rho*Strain-2/3._dp*rho*Ident*Omega)*dVelodx) - &
             c3 * rho_g / Pr
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
          A(2,2) = A(2,2) + rho * F4 * Beta * Omega * Basis(q) * Basis(p)
!------------------------------------------------------------------------------
!         The diffusion term
!------------------------------------------------------------------------------
          IF ( CurrentCoordinateSystem() == Cartesian ) THEN
             DO i=1,dim
               A(1,1) = A(1,1) + Effmu(1) * dBasisdx(q,i) * dBasisdx(p,i)
               A(2,2) = A(2,2) + Effmu(2) * dBasisdx(q,i) * dBasisdx(p,i)
               A(2,2) = A(2,2) - 2*rho*(1-F1)/1.168_dp/Omega*GradK(i)*dBasisdx(q,i)*Basis(p)
             END DO
          ELSE
             DO i=1,dim
               DO j=1,dim
                  A(1,1) = A(1,1) + Metric(i,j) * Effmu(1) * &
                       dBasisdx(q,i) * dBasisdx(p,i)

                  A(2,2) = A(2,2) + Metric(i,j) * Effmu(2) * &
                       dBasisdx(q,i) * dBasisdx(p,i)

                  A(2,2) = A(2,2) - 2*rho*(1-F1)/1.168_dp/Omega*Metric(i,j)* &
                       GradK(i)*dBasisdx(q,i)*Basis(p)
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
        LoadAtIP(1) = ProdK
        LoadAtIP(2) = rGamma*ProdO

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
   SUBROUTINE OmegaWall( Element,n )
!------------------------------------------------------------------------------
     TYPE(Element_t), POINTER :: Element
     INTEGER :: n
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Distance(32), omega_wall, dist, mu(32), rho(32)
     INTEGER :: i,j,np
     TYPE(Element_t), POINTER :: Parent
!------------------------------------------------------------------------------
     Parent => Element % BoundaryInfo % Left
     IF ( .NOT. ASSOCIATED(Parent) ) THEN
       Parent => Element % BoundaryInfo % Right
     ELSE
       IF ( .NOT. ALL(Solver % Variable % Perm(Parent % NodeIndexes)>0) ) &
         Parent => Element % BoundaryInfo % Right
     END IF
     IF(.NOT.ASSOCIATED(Parent))RETURN

     np = GetElementNOFNodes(Parent)

     rho(1:np)= GetReal( GetMaterial(Parent), 'Density', UElement=Parent )
     mu(1:np) = GetReal( GetMaterial(Parent), 'Viscosity', UElement=Parent )

     CALL GetScalarLocalSolution( Distance, 'Wall distance', Parent )

     omega_wall = 1.d10
     DO i=1,np
       j = Parent % NodeIndexes(i)
       IF ( Distance(i) < AEPS ) CYCLE
       IF ( ANY( j==Element % NodeIndexes(1:n) ) ) CYCLE

!      omega_wall = 2*mu(i)/0.09_dp/rho(i)/Distance(i)**2
       omega_wall = 6*mu(i)/rho(i)/0.075_dp/Distance(i)**2

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
  END SUBROUTINE SSTKOmega
!------------------------------------------------------------------------------
