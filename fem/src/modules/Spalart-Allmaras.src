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
! *  Original Date: 16 Nov 1997
! *
! ****************************************************************************/

!------------------------------------------------------------------------------
!> Solver for the Spalart-Allmaras-turbulence model.
!> \ingroup Solvers
!------------------------------------------------------------------------------
   SUBROUTINE SpalartAllmaras( Model,Solver,dt,TransientSimulation )
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

     REAL(KIND=dp) :: RelativeChange,Norm
     LOGICAL :: Stabilize = .TRUE.,gotIt

     LOGICAL :: AllocationsDone = .FALSE.

     TYPE(Variable_t), POINTER :: FlowSol, KE

     INTEGER, POINTER :: KinPerm(:)

     INTEGER :: NonlinearIter

     REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), &
       STIFF(:,:), LOAD(:,:),FORCE(:), LocalKinEnergy(:), TimeForce(:)

     TYPE(ValueList_t), POINTER :: BC, Equation, Material

     SAVE MASS,STIFF,LOAD,FORCE, ElementNodes,AllocationsDone,TimeForce

     REAL(KIND=dp) :: at,at0,CPUTime,RealTime, KMax, EMax, KVal, EVal
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
         CALL Fatal( 'SpalartAllmaras', 'Memory allocation error.' )
       END IF

       AllocationsDone = .TRUE.
     END IF

!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------

     NonlinearIter = ListGetInteger( Solver % Values, &
         'Nonlinear System Max Iterations',GotIt )

     IF ( .NOT.GotIt ) NonlinearIter = 1

!------------------------------------------------------------------------------
      DO i=1,Model % NumberOFBCs
        BC => Model % BCs(i) % Values
        IF ( GetLogical( BC, 'Noslip wall BC', gotit ) ) THEN
          CALL ListAddConstReal( BC, 'Turbulent Viscosity', 0.0_dp )
        END IF
      END DO
!------------------------------------------------------------------------------

     DO iter=1,NonlinearIter

       at  = CPUTime()
       at0 = RealTime()

       CALL Info( 'SpalartAllmaras', ' ', Level=4 )
       CALL Info( 'SpalartAllmaras', ' ', Level=4 )
       CALL Info( 'SpalartAllmaras', &
          '-------------------------------------', Level=4 )
       WRITE( Message, * ) 'Spalart-Allmaras iteration: ', iter
       CALL Info( 'SpalartAllmaras', Message, Level=4 )
       CALL Info( 'SpalartAllmaras', &
          '-------------------------------------', Level=4 )
       CALL Info( 'SpalartAllmaras', ' ', Level=4 )
       CALL Info( 'SpalartAllmaras', 'Starting Assembly...', Level=4 )

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

           CALL Info( 'SpalartAllmaras', Message, Level=5 )
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
      CALL Info( 'SpalartAllmaras', 'Assembly done', Level=4 )

!------------------------------------------------------------------------------
      CALL DefaultFinishAssembly()
      CALL DefaultDirichletBCs()
!------------------------------------------------------------------------------
      CALL Info( 'SpalartAllmaras', 'Set boundaries done', Level=4 )
!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
      Norm = DefaultSolve()
!------------------------------------------------------------------------------
!     Kinetic Energy Solution should be positive
!------------------------------------------------------------------------------
      n = SIZE( Solver % Variable % Values)
      DO i=1,SIZE(Solver % Variable % Perm)
         k = Solver % Variable % Perm(i)
         IF ( k <= 0 ) CYCLE
         Kval = Solver % Variable % Values(k)
         Solver % Variable % Values(k) = MAX( KVal, 1.0d-12 )
      END DO

!------------------------------------------------------------------------------

      IF ( Solver % Variable % NonlinConverged == 1 ) EXIT
!------------------------------------------------------------------------------
    END DO
!------------------------------------------------------------------------------

CONTAINS

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
!******************************************************************************
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

     REAL(KIND=dp) :: UX(n), UY(n), UZ(n), Velo(3), dVelodx(3,3),Tviscosity(n), &
                      Distance(n), Density(n), Viscosity(n)

     REAL(KIND=dp) :: A,M,Prod,div
     INTEGER :: i,j,c,p,q,t,dim,NBasis
     REAL(KIND=dp) :: LoadatIp,Cmu,Rho,mu,Tmu,Effmu

     REAL(KIND=dp) :: s,u,v,w,Strain(3,3), Vorticity(3,3), dist

     REAL(KIND=dp) :: StrainMeasure,VorticityMeasure,X,Y,Z,Sigma, &
        GradTmu(3), Cw1,Cw2,Cw3,fw,fw1,fw2,Cb1,Cb2,Cb3,St,Xi,Cv1,g,r

     REAL(KIND=dp) :: Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3),SqrtMetric

     LOGICAL :: stat, Bubbles
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

!------------------------------------------------------------------------------

     dim = CoordinateSystemDimension()

     Viscosity(1:n) = GetReal( Material, 'Viscosity' )
     Density(1:n) = GetReal( Material, 'Density' )

     CALL GetScalarLocalSolution( UX, 'Velocity 1' )
     CALL GetScalarLocalSolution( UY, 'Velocity 2' )
     CALL GetScalarLocalSolution( UZ, 'Velocity 3' )

     CALL GetScalarLocalSolution( TViscosity )
     CALL GetScalarLocalSolution( Distance, 'Wall Distance' )

     FORCE = 0.0_dp
     STIFF = 0.0_dp
     MASS  = 0.0_dp

     NBasis = 2*n
     Bubbles = .TRUE.

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
         StrainMeasure = SQRT(2 * SUM(Strain * Strain))

         Vorticity = 0.5_dp * (dVelodx - TRANSPOSE(dVelodx))
         VorticityMeasure = SQRT(2 * SUM(Vorticity * Vorticity))
       ELSE
         StrainMeasure = SQRT(SecondInvariant( Velo,dVelodx,Metric,Symb )/2)
       END IF

!------------------------------------------------------------------------------

       Tmu = SUM( Tviscosity(1:n) * Basis(1:n) )
       DO i=1,dim
         GradTmu(i) = SUM( Tviscosity(1:n) * dBasisdx(1:n,i) )
       END DO

       mu   = SUM( Viscosity(1:n) * Basis(1:n) )
       rho  = SUM( Density(1:n) * Basis(1:n) )
       dist = SUM( Distance(1:n) * Basis(1:n) )

       Cb1 = 0.1335_dp
       Cb2 = 0.6220_dp
       Cv1 = 7.1_dp
       Sigma = 2._dp/3._dp

       ! Rotation and curvature correction of Schweighofer & Helsten; NOT IN USE
       ! -----------------------------------------------------------------------
!      r = VorticityMeasure/StrainMeasure*(VorticityMeasure/StrainMeasure-1)
!      r = 1._dp / (1+3.6_dp*r)
       r = 1._dp
       Cw1 = r*(Cb1/0.41_dp**2 + (1+Cb2)/Sigma)

       Cw2 = 0.3_dp
       Cw3 = 2.0_dp

       Xi  = Tmu/(mu/rho)
       fw1 = Xi**3 / (Xi**3 + Cv1**3)
       fw2 = 1 - Xi / ( 1+Xi*fw1 )

       St = VorticityMeasure + 2 * MIN(0.0_dp, StrainMeasure-VorticityMeasure)
       St = St + Tmu / dist**2 / 0.41_dp**2 * fw2

       r  = Tmu / St / 0.41_dp**2 / dist**2
       g  = r + Cw2 * (r**6-r)
       fw = g*((1+Cw3**6)/(g**6+Cw3**6))**(1._dp/6._dp)

       Effmu = (mu + rho*Tmu)/Sigma
!------------------------------------------------------------------------------
!      Loop over basis functions of both unknowns and weights
!------------------------------------------------------------------------------
       DO p=1,NBasis
       DO q=1,NBasis
          M = 0.0d0
          A = 0.0d0

          M = rho * Basis(q) * Basis(p)
          A = A - rho * Cb1 * St * Basis(q) * Basis(p)/4
          A = A + rho * Cw1 * fw * Tmu / dist**2 * Basis(q) * Basis(p)
!------------------------------------------------------------------------------
!         The diffusion term
!------------------------------------------------------------------------------
          IF ( CurrentCoordinateSystem() == Cartesian ) THEN
             DO i=1,dim
               A = A + Effmu * dBasisdx(q,i) * dBasisdx(p,i) 
             END DO
          ELSE
             DO i=1,dim
               DO j=1,dim
                  A = A + Metric(i,j) * Effmu * dBasisdx(q,i) * dBasisdx(p,i)
               END DO
             END DO
          END IF

!------------------------------------------------------------------------------
!           The convection term
!------------------------------------------------------------------------------
          DO i=1,dim
            A = A + rho * (Velo(i)-Cb2*GradTmu(i)/Sigma) * dBasisdx(q,i) * Basis(p)
          END DO

          MASS(p,q)  = MASS(p,q)  + s*M
          STIFF(p,q) = STIFF(p,q) + s*A
        END DO
        END DO

        ! Load at the integration point:
        !-------------------------------
        LoadAtIp = 3._dp*rho * Cb1 * St * Tmu/4

!------------------------------------------------------------------------------
        DO p=1,NBasis
          FORCE(p) = FORCE(p)+s*LoadAtIp*Basis(p)
        END DO
      END DO
!------------------------------------------------------------------------------
   END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  END SUBROUTINE SpalartAllmaras
!------------------------------------------------------------------------------
