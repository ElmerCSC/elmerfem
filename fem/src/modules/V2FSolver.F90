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
! *  Module containing a solver for the V2-F-turbulence model.
! *
! ******************************************************************************
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
!> Solve the V^2-F turbulence model equations. Is this some other variant of it?
!> \ingroup Solvers
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
   SUBROUTINE V2F_LDM( Model,Solver,dt,TransientSimulation )
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
     INTEGER :: i,j,t,n,nd,iter,istat,LocalNodes,DOFs,body_id

     TYPE(Nodes_t)   :: ElementNodes
     TYPE(Element_t),POINTER :: Element

     REAL(KIND=dp) :: RelativeChange,Norm,S,C

     LOGICAL :: Stabilize = .TRUE.,NewtonLinearization = .FALSE.,gotIt
!
     LOGICAL :: AllocationsDone = .FALSE.

     CHARACTER(LEN=MAX_NAME_LEN) :: KEModel, V2FModel

     INTEGER :: NSDOFs,NewtonIter,NonlinearIter
     REAL(KIND=dp) :: NewtonTol, Clip, V2FCp

     REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), &
       STIFF(:,:), LOAD(:,:),FORCE(:),U(:),V(:),W(:), &
           Density(:),Viscosity(:),KinDis(:), &
                 KinEne(:),KESigmaK(:),KESigmaE(:),KECmu(:),KEC1(:),&
                   KEC2(:), TimeForce(:), V2(:), f(:), &
                    V2FCnu(:), V2FC1(:), V2FC2(:), V2FCL(:), V2FCT(:), V2FSigma(:)

     TYPE(ValueList_t), POINTER :: BC, Equation, Material

     SAVE U,V,W,MASS,STIFF,LOAD,FORCE, ElementNodes,Density,&
         AllocationsDone,Viscosity,LocalNodes, &
           KinDis,KinEne,KESigmaK,KESigmaE,KECmu, &
             TimeForce, KEC1, KEC2, &
               V2, F, V2FCnu, V2FC1, V2FC2, V2FCL, V2FCT, V2FSigma
     REAL(KIND=dp) :: at,at0,KMax, EMax, KVal, EVal
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
     IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

     LocalNodes = COUNT( Solver % Variable % Perm > 0 )
     IF ( LocalNodes <= 0 ) RETURN

     DOFs = Solver % Variable % DOFs
!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone ) THEN
       N = SOlver % Mesh % MaxElementDOFs

       ALLOCATE( U( N ), V( N ), W( N ),  &
                 Density( N ), &
                 Viscosity(N), &
                 KEC1(N),      &
                 KEC2(N),      &
                 KESigmaK(N),  &
                 KESigmaE(N),  &
                 KECmu(N),     &
                 V2FCnu(N),    &
                 V2FC1(N),     &
                 V2FC2(N),     &
                 V2FCL(N),     &
                 V2FCT(N),     &
                 V2FSigma(N),  &
                 KinEne( N ),  &
                 KinDis( N ),  &
                 V2( N ),      &
                 F( N ),       &
                 MASS( 2*DOFs*N,2*DOFs*N ), &
                 STIFF( 2*DOFs*N,2*DOFs*N ),LOAD( DOFs,N ), &
                 FORCE( 2*DOFs*N ), TimeForce( 2*DOFs*N ), STAT=istat )

       IF ( istat /= 0 ) THEN
         CALL Fatal( 'V2-F-Solver', 'Memory allocation error.' )
       END IF

       AllocationsDone = .TRUE.
     END IF

!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------
!    Stabilize = ListGetLogical( Solver % Values,'Stabilize',GotIt )
     Stabilize = .FALSE.

     NewtonTol = ListGetConstReal( Solver % Values, &
        'Nonlinear System Newton After Tolerance',gotIt )

     NewtonIter = ListGetInteger( Solver % Values, &
        'Nonlinear System Newton After Iterations',gotIt )

     NonlinearIter = ListGetInteger( Solver % Values, &
         'Nonlinear System Max Iterations',GotIt )

     IF ( .NOT.GotIt ) NonlinearIter = 1

!------------------------------------------------------------------------------

     DO iter=1,NonlinearIter

       at  = CPUTime()
       at0 = RealTime()

       CALL Info( 'V2-F-Solver', ' ', Level=4 )
       CALL Info( 'V2-F-Solver', ' ', Level=4 )
       CALL Info( 'V2-F-Solver', &
          '-------------------------------------', Level=4 )
       WRITE( Message, * ) 'V2-F iteration: ', iter
       CALL Info( 'V2-F-Solver', Message, Level=4 )
       CALL Info( 'V2-F-Solver', &
          '-------------------------------------', Level=4 )
       CALL Info( 'V2-F-Solver', ' ', Level=4 )
       CALL Info( 'V2-F-Solver', 'Starting Assembly...', Level=4 )

       CALL DefaultInitialize()

!------------------------------------------------------------------------------
!      Bulk elements
!------------------------------------------------------------------------------
       body_id = -1
       CALL StartAdvanceOutput( 'V2FSolver','Assembly:')
       DO t=1,Solver % NumberOfActiveElements

!------------------------------------------------------------------------------

         CALL AdvanceOutput(t,Solver % NumberOFActiveElements)

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
         END IF
!------------------------------------------------------------------------------
         CALL GetElementNodes( ElementNodes )
         n  = GetElementNOFNodes()
         nd = GetElementNOFDOFs()
         nd = 2*n
!------------------------------------------------------------------------------
         V2FModel = GetString( Material, 'V2-F Model', GotIt )
         IF ( .NOT. GotIt ) V2FModel = 'ldm'

         V2FSigma(1:n) = GetReal( Material, 'V2-F Sigma', GotIt )
         IF ( .NOT. GotIt ) THEN
           V2FSigma(1:n) = 1.0d0
           CALL ListAddConstReal( Material, 'V2-F Sigma', V2FSigma(1) )
         END IF

         KECmu(1:n) = ListGetConstReal( Material, 'KE Cmu', GotIt )
         IF ( .NOT. GotIt ) THEN
            KECmu(1:n) = 0.22_dp 
            CALL ListAddConstReal( Material, 'KE Cmu', KECmu(1) )
         END IF

         V2FC1(1:n) = GetReal( Material, 'V2-F C1', GotIt )
         IF ( .NOT. GotIt ) THEN
            V2FC1(1:n) = 1.4_dp
            CALL ListAddConstReal( Material, 'V2-F C1', V2FC1(1) )
         END IF

         V2FC2(1:n) = GetReal( Material, 'V2-F C2', GotIt )
         IF ( .NOT. GotIt ) THEN
            V2FC2(1:n) = 0.3_dp
            CALL ListAddConstReal( Material, 'V2-F C2', V2FC2(1) )
         END IF

         V2FCT(1:n) = GetReal( Material, 'V2-F CT', GotIt )
         IF ( .NOT. GotIt ) THEN
            V2FCT(1:n) = 6.0_dp
            CALL ListAddConstReal( Material, 'V2-F CT', V2FCT(1) )
         END IF

         V2FCL(1:n) = GetReal( Material, 'V2-F CL', GotIt )
         IF ( .NOT. GotIt ) THEN
            V2FCL(1:n) = 0.23_dp
            CALL ListAddConstReal( Material, 'V2-F CL', V2FCL(1) )
         END IF

         V2FCnu(1:n) = GetReal( Material, 'V2-F Cnu', GotIt )
         IF ( .NOT. GotIt ) THEN
            V2FCnu(1:n) = 70.0d0
            CALL ListAddConstReal( Material, 'V2-F Cnu', V2FCnu(1) )
         END IF

         V2FCp = GetCReal( Material, 'V2-F Cp', GotIt )
         IF ( .NOT. GotIt ) THEN
            CALL ListAddConstReal( Material, 'V2-F Cp', 0.05_dp )
         END IF
!------------------------------------------------------------------------------
         Density(1:n)   = GetReal( Material,'Density' )
         Viscosity(1:n) = GetReal( Material,'Viscosity' )
!------------------------------------------------------------------------------
         CALL GetScalarLocalSolution( V2, 'V2' )
         CALL GetScalarLocalSolution( F, 'F' )
         CALL GetScalarLocalSolution( KinEne, 'Kinetic Energy' )
         CALL GetScalarLocalSolution( KinDis, 'Kinetic Dissipation' )

         CALL GetScalarLocalSolution( U, 'Velocity 1' )
         CALL GetScalarLocalSolution( V, 'Velocity 2' )
         CALL GetScalarLocalSolution( W, 'Velocity 3' )
!------------------------------------------------------------------------------
!        Get element local matrices, and RHS vectors
!------------------------------------------------------------------------------
         CALL LocalMatrix( MASS,STIFF,FORCE,LOAD, &
             U,V,W,Element,n,nd,ElementNodes )
!------------------------------------------------------------------------------
         TimeForce = 0.0_dp
         IF ( TransientSimulation ) CALL Default1stOrderTime(MASS,STIFF,FORCE)
!------------------------------------------------------------------------------
!        Update global matrices from local matrices
!------------------------------------------------------------------------------
         CALL Condensate( DOFs*n, STIFF, FORCE, TimeForce )
         CALL DefaultUpdateEquations( STIFF, FORCE )
!------------------------------------------------------------------------------
      END DO     !  Bulk elements
      CALL DefaultFinishBulkAssembly()
      CALL Info( 'V2-F-Solver', 'Assembly done', Level=4 )

!------------------------------------------------------------------------------

      DO i=1,Model % NumberOFBCs
        BC => Model % BCs(i) % Values
        IF ( GetLogical(  BC, 'Noslip wall BC', gotit ) ) THEN
          CALL ListAddConstReal( BC, 'F',  0.0_dp )
          CALL ListAddConstReal( BC, 'V2', 0.0_dp )
        END IF
      END DO

      CALL DefaultFinishAssembly()
!------------------------------------------------------------------------------
!     Dirichlet boundary conditions
!------------------------------------------------------------------------------
      CALL DefaultDirichletBCs()
!------------------------------------------------------------------------------
      CALL Info( 'V2-F-Solver', 'Set boundaries done', Level=4 )
!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
      Norm = DefaultSolve()

      Solver % Variable % Values(1::2) = &
        MAX( Solver % Variable % Values(1::2), 1.0d-9 )

!------------------------------------------------------------------------------
      WRITE( Message,* ) 'Result Norm   : ',Norm
      CALL Info( 'V2-F-Solver', Message, Level = 4 )

       RelativeChange = Solver % Variable % NonlinChange
      WRITE( Message,* ) 'Relative Change : ',RelativeChange
      CALL Info( 'V2-F-Solver', Message, Level = 4 )

      IF ( RelativeChange < NewtonTol .OR. &
              iter > NewtonIter ) NewtonLinearization = .TRUE.

      IF ( Solver % Variable % NonlinConverged == 1 ) EXIT
!------------------------------------------------------------------------------
    END DO
!------------------------------------------------------------------------------

CONTAINS

!------------------------------------------------------------------------------
   SUBROUTINE LocalMatrix( MASS,STIFF,FORCE, &
            LOAD,UX,UY,UZ,Element,n,nd,Nodes )
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
!  REAL(KIND=dp) :: UX(:),UY(:),UZ(:)
!     INPUT: Nodal values of velocity components from previous iteration
!           used only if coefficient of the convection term (C1) is nonzero
!
!  LOGICAL :: Stabilize
!     INPUT: Should stabilzation be used ? Used only if coefficient of the
!            convection term (C1) is nonzero
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

     REAL(KIND=dp), DIMENSION(:)   :: FORCE,UX,UY,UZ
     REAL(KIND=dp), DIMENSION(:,:) :: MASS,STIFF,LOAD

     LOGICAL :: Stabilize

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

     REAL(KIND=dp) :: A(2,2),M(2,2)
     REAL(KIND=dp) :: LoadatIp(2),Rho,mu,Tmu,EffVisc
     INTEGER :: i,j,c,p,q,t,dim,N_Integ,NBasis

     REAL(KIND=dp) :: s,u,v,w, K,E,Eta,LV2,Lf,Strain(3,3),ProdTensor(3,3), &
                     Vorticity(3,3), nu, Cv, Prod,div

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     REAL(KIND=dp) :: C1,C2,CT,CL,Cnu,Cmu,TimeScale,LengthScale2,aparm(nd)
     REAL(KIND=dp) :: SecInv,X,Y,Z,Re_T
     REAL(KIND=dp) :: Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3),SqrtMetric

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     LOGICAL :: stat,Convection, Bubbles

!------------------------------------------------------------------------------

     dim = CoordinateSystemDimension()

     FORCE = 0.0_dp
     STIFF = 0.0_dp
     MASS  = 0.0_dp

     Bubbles = .TRUE.
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IF ( Bubbles ) THEN
        IntegStuff = GaussPoints( element, element % TYPE % GaussPoints2 )
     ELSE
        IntegStuff = GaussPoints( element )
     END IF

     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!    Now we start integrating
!------------------------------------------------------------------------------
     DO t=1,N_Integ
       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)
!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,detJ, &
               Basis,dBasisdx,Bubbles=.TRUE. )
!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
       s = detJ * S_Integ(t)
       IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
         X = SUM( Nodes % x(1:n)*Basis(1:n) )
         Y = SUM( Nodes % y(1:n)*Basis(1:n) )
         Z = SUM( nodes % z(1:n)*Basis(1:n) )
         CALL CoordinateSystemInfo(Metric,SqrtMetric,Symb,dSymb,X,Y,Z)

         s = s * SqrtMetric
       END IF

!      Velocity from previous iteration at the integration point
!------------------------------------------------------------------------------
       Velo = 0.0d0
       Velo(1) = SUM( UX(1:n)*Basis(1:n) )
       Velo(2) = SUM( UY(1:n)*Basis(1:n) )
       Velo(3) = SUM( UZ(1:n)*Basis(1:n) )

       dVelodx = 0.0d0
       DO i=1,dim
         dVelodx(1,i) = SUM( UX(1:n)*dBasisdx(1:n,i) )
         dVelodx(2,i) = SUM( UY(1:n)*dBasisdx(1:n,i) )
         dVelodx(3,i) = SUM( UZ(1:n)*dBasisdx(1:n,i) )
       END DO

       Strain    = 0.5_dp * ( dVelodx + TRANSPOSE(dVelodx) )
       Vorticity = 0.5_dp * ( dVelodx - TRANSPOSE(dVelodx) )

       IF ( CurrentCoordinateSystem() == Cartesian ) THEN
          Secinv = 2*SUM(Strain * Strain)
       ELSE
          SecInv = SecondInvariant(Velo,dVelodx,Metric,Symb)/2
       END IF

!------------------------------------------------------------------------------
!      Coefficient of the convection and time derivative terms
!      at the integration point
!------------------------------------------------------------------------------
       K   = SUM( KinEne(1:n) * Basis(1:n) )
       E   = SUM( KinDis(1:n) * Basis(1:n) )
       LF  = SUM( F(1:n)  * Basis(1:n) )
       LV2 = SUM( V2(1:n) * Basis(1:n) )

       Cnu = SUM( V2FCnu(1:n) * Basis(1:n) )
       Cmu = SUM( KECMu(1:n)  * Basis(1:n) )
       CT  = SUM( V2FCT(1:n)  * Basis(1:n) )
       CL  = SUM( V2FCL(1:n)  * Basis(1:n) )
       C1  = SUM( V2FC1(1:n)  * Basis(1:n) )
       C2  = SUM( V2FC2(1:n)  * Basis(1:n) )

       rho = SUM( Basis(1:n) * Density(1:n) )
       mu  = SUM( Basis(1:n) * Viscosity(1:n) )

       Re_T = K**2 / ((mu/Rho)*E)
       TimeScale = MAX( K/E, CT*SQRT((mu/rho)/E))
       Lengthscale2 = CL**2 * MAX(K**3 /E**2, Cnu**2*SQRT((mu/rho)**3/E)) 

       Tmu = rho*Cmu*LV2*TimeScale
       EffVisc = mu + Tmu / SUM(V2FSigma(1:n)*Basis(1:n))

!      div = Strain(1,1) + Strain(2,2) + Strain(3,3)
!      ProdTensor = 2*Tmu*Strain
!      DO i=1,3
!        ProdTensor(i,i) = ProdTensor(i,i) - 2.0_dp*(Tmu*div+K)/3.0_dp
!      END DO
!      Prod = SUM(ProdTensor * dVelodx) / Rho
       Prod = Tmu * SecInv / Rho

!------------------------------------------------------------------------------
!      Loop Over basis functions of both unknowns and weights
!------------------------------------------------------------------------------
       DO p=1,nd
       DO q=1,nd
!------------------------------------------------------------------------------
!         The diffusive-convective equation without stabilization
!------------------------------------------------------------------------------
          M = 0.0d0
          A = 0.0d0

          M(1,1) = rho * Basis(q) * Basis(p)

          A(1,1) = A(1,1) + 6 * Rho / TimeScale  * Basis(q) * Basis(p)
          A(1,2) = A(1,2) - Rho * K * Basis(q) * Basis(p)

          A(2,1) = A(2,1) + (C1-6) / K / TimeScale * Basis(q) * Basis(p)
          A(2,2) = A(2,2) + Basis(q) * Basis(p)

          ! The diffusion term:
          !--------------------
          IF ( CurrentCoordinateSystem() == Cartesian ) THEN
             DO i=1,dim
               A(1,1) = A(1,1) + EffVisc * dBasisdx(q,i) * dBasisdx(p,i)
               A(2,2) = A(2,2) + LengthScale2*dBasisdx(q,i) * dBasisdx(p,i)
             END DO
          ELSE
             DO i=1,dim
               DO j=1,dim
                  A(1,1) = A(1,1) + Metric(i,j) * EffVisc * &
                       dBasisdx(q,i) * dBasisdx(p,i)

                  A(2,2) = A(2,2) + Metric(i,j) * LengthScale2 * &
                       dBasisdx(q,i) * dBasisdx(p,i)
               END DO
             END DO
          END IF
          ! The convection term:
          !---------------------
          DO i=1,dim
            A(1,1) = A(1,1) + rho * Velo(i) * dBasisdx(q,i) * Basis(p)
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
       LoadAtIP = 0.0_dp
       LoadAtIP(2) = ((C1-1)*2/3.0_dp+C2*Prod/E)/TimeScale
!------------------------------------------------------------------------------
        DO p=1,nd
           FORCE(2*(p-1)+1) = FORCE(2*(p-1)+1) + s*LoadAtIp(1)*Basis(p)
           FORCE(2*(p-1)+2) = FORCE(2*(p-1)+2) + s*LoadAtIp(2)*Basis(p)
        END DO
      END DO
!------------------------------------------------------------------------------
   END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  END SUBROUTINE V2F_LDM
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  Solve the V^2-F turbulence model equations.
!> \ingroup Solvers
!------------------------------------------------------------------------------
   SUBROUTINE V2F( Model,Solver,dt,TransientSimulation )
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
     INTEGER :: i,j,t,n,nd,iter,istat,LocalNodes,DOFs,body_id

     TYPE(Nodes_t)   :: ElementNodes
     TYPE(Element_t),POINTER :: Element

     REAL(KIND=dp) :: RelativeChange,Norm,S,C

     LOGICAL :: Stabilize = .TRUE.,NewtonLinearization = .FALSE.,gotIt
!
     LOGICAL :: AllocationsDone = .FALSE.

     CHARACTER(LEN=MAX_NAME_LEN) :: KEModel, V2FModel

     INTEGER :: NSDOFs,NewtonIter,NonlinearIter
     REAL(KIND=dp) :: NewtonTol, Clip, V2FCp

     REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), &
       STIFF(:,:), LOAD(:,:),FORCE(:),U(:),V(:),W(:), &
           Density(:),Viscosity(:),KinDis(:), &
                 KinEne(:),KESigmaK(:),KESigmaE(:),KECmu(:),KEC1(:),&
                   KEC2(:), TimeForce(:), V2(:), f(:), &
                    V2FCnu(:), V2FC1(:), V2FC2(:), V2FCL(:), V2FCT(:), V2FSigma(:)

     TYPE(ValueList_t), POINTER :: BC, Equation, Material

     SAVE U,V,W,MASS,STIFF,LOAD,FORCE, ElementNodes,Density,&
         AllocationsDone,Viscosity,LocalNodes, &
           KinDis,KinEne,KESigmaK,KESigmaE,KECmu, &
             TimeForce, KEC1, KEC2, &
               V2, F, V2FCnu, V2FC1, V2FC2, V2FCL, V2FCT, V2FSigma
     REAL(KIND=dp) :: at,at0, KMax, EMax, KVal, EVal
!------------------------------------------------------------------------------
     CHARACTER(LEN=MAX_NAME_LEN) :: VersionID = "$Id$"

!------------------------------------------------------------------------------
!    Check if version number output is requested
!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone ) THEN
        IF ( ListGetLogical( GetSimulation(), 'Output Version Numbers', GotIt ) ) THEN
           CALL Info( 'V2-F-Solver', 'V2-F Solver version:', Level = 0 ) 
           CALL Info( 'V2-F-Solver', VersionID, Level = 0 ) 
           CALL Info( 'V2-F-Solver', ' ', Level = 0 ) 
        END IF
     END IF

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
     IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

     LocalNodes = COUNT( Solver % Variable % Perm > 0 )
     IF ( LocalNodes <= 0 ) RETURN

     DOFs = Solver % Variable % DOFs
!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone ) THEN
       N = SOlver % Mesh % MaxElementDOFs

       ALLOCATE( U( N ), V( N ), W( N ),  &
                 Density( N ), &
                 Viscosity(N), &
                 KEC1(N),      &
                 KEC2(N),      &
                 KESigmaK(N),  &
                 KESigmaE(N),  &
                 KECmu(N),     &
                 V2FCnu(N),    &
                 V2FC1(N),     &
                 V2FC2(N),     &
                 V2FCL(N),     &
                 V2FCT(N),     &
                 V2FSigma(N),  &
                 KinEne( N ),  &
                 KinDis( N ),  &
                 V2( N ),      &
                 F( N ),       &
                 MASS( 2*DOFs*N,2*DOFs*N ), &
                 STIFF( 2*DOFs*N,2*DOFs*N ),LOAD( DOFs,N ), &
                 FORCE( 2*DOFs*N ), TimeForce(2*DOFs*N), STAT=istat )

       IF ( istat /= 0 ) THEN
         CALL Fatal( 'V2-F-Solver', 'Memory allocation error.' )
       END IF

       AllocationsDone = .TRUE.
     END IF

!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------
!    Stabilize = ListGetLogical( Solver % Values,'Stabilize',GotIt )
     Stabilize = .FALSE.

     NewtonTol = ListGetConstReal( Solver % Values, &
        'Nonlinear System Newton After Tolerance',gotIt )

     NewtonIter = ListGetInteger( Solver % Values, &
        'Nonlinear System Newton After Iterations',gotIt )

     NonlinearIter = ListGetInteger( Solver % Values, &
         'Nonlinear System Max Iterations',GotIt )

     IF ( .NOT.GotIt ) NonlinearIter = 1

!------------------------------------------------------------------------------

     DO iter=1,NonlinearIter

       at  = CPUTime()
       at0 = RealTime()

       CALL Info( 'V2-F-Solver', ' ', Level=4 )
       CALL Info( 'V2-F-Solver', ' ', Level=4 )
       CALL Info( 'V2-F-Solver', &
          '-------------------------------------', Level=4 )
       WRITE( Message, * ) 'V2-F iteration: ', iter
       CALL Info( 'V2-F-Solver', Message, Level=4 )
       CALL Info( 'V2-F-Solver', &
          '-------------------------------------', Level=4 )
       CALL Info( 'V2-F-Solver', ' ', Level=4 )
       CALL Info( 'V2-F-Solver', 'Starting Assembly...', Level=4 )

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

           CALL Info( 'V2-F-Solver', Message, Level=5 )
           at0 =RealTime()
         END IF
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
         END IF
!------------------------------------------------------------------------------
         CALL GetElementNodes( ElementNodes )
         n  = GetElementNOFNodes()
         nd = GetElementNOFDOFs()
         nd = 2*n
!------------------------------------------------------------------------------
         V2FModel = GetString( Material, 'V2-F Model', GotIt )
         IF ( .NOT. GotIt ) V2FModel = 'ldm'

         V2FSigma(1:n) = GetReal( Material, 'V2-F Sigma', GotIt )
         IF ( .NOT. GotIt ) THEN
           V2FSigma(1:n) = 1.0d0
           CALL ListAddConstReal( Material, 'V2-F Sigma', V2FSigma(1) )
         END IF

         KECmu(1:n) = ListGetConstReal( Material, 'KE Cmu', GotIt )
         IF ( .NOT. GotIt ) THEN
            KECmu(1:n) = 0.19_dp 
            CALL ListAddConstReal( Material, 'KE Cmu', KECmu(1) )
         END IF

         V2FC1(1:n) = GetReal( Material, 'V2-F C1', GotIt )
         IF ( .NOT. GotIt ) THEN
            V2FC1(1:n) = 1.4_dp
            CALL ListAddConstReal( Material, 'V2-F C1', V2FC1(1) )
         END IF

         V2FC2(1:n) = GetReal( Material, 'V2-F C2', GotIt )
         IF ( .NOT. GotIt ) THEN
            V2FC2(1:n) = 0.3_dp
            CALL ListAddConstReal( Material, 'V2-F C2', V2FC2(1) )
         END IF

         V2FCT(1:n) = GetReal( Material, 'V2-F CT', GotIt )
         IF ( .NOT. GotIt ) THEN
            V2FCT(1:n) = 6.0_dp
            CALL ListAddConstReal( Material, 'V2-F CT', V2FCT(1) )
         END IF

         V2FCL(1:n) = GetReal( Material, 'V2-F CL', GotIt )
         IF ( .NOT. GotIt ) THEN
            V2FCL(1:n) = 0.3_dp
            CALL ListAddConstReal( Material, 'V2-F CL', V2FCL(1) )
         END IF

         V2FCnu(1:n) = GetReal( Material, 'V2-F Cnu', GotIt )
         IF ( .NOT. GotIt ) THEN
            V2FCnu(1:n) = 70.0d0
            CALL ListAddConstReal( Material, 'V2-F Cnu', V2FCnu(1) )
         END IF

         V2FCp = GetCReal( Material, 'V2-F Cp', GotIt )
         IF ( .NOT. GotIt ) THEN
            CALL ListAddConstReal( Material, 'V2-F Cp', 0.045_dp )
         END IF
!------------------------------------------------------------------------------
         Density(1:n)   = GetReal( Material,'Density' )
         Viscosity(1:n) = GetReal( Material,'Viscosity' )
!------------------------------------------------------------------------------
         CALL GetScalarLocalSolution( V2, 'V2' )
         CALL GetScalarLocalSolution( F, 'F' )
         CALL GetScalarLocalSolution( KinEne, 'Kinetic Energy' )
         CALL GetScalarLocalSolution( KinDis, 'Kinetic Dissipation' )

         CALL GetScalarLocalSolution( U, 'Velocity 1' )
         CALL GetScalarLocalSolution( V, 'Velocity 2' )
         CALL GetScalarLocalSolution( W, 'Velocity 3' )
!------------------------------------------------------------------------------
!        Get element local matrices, and RHS vectors
!------------------------------------------------------------------------------
         CALL LocalMatrix( MASS,STIFF,FORCE,LOAD, &
             U,V,W,Element,n,nd,ElementNodes )
!------------------------------------------------------------------------------
         TimeForce = 0.0_dp
         IF ( TransientSimulation ) CALL Default1stOrderTime(MASS,STIFF,FORCE)
!------------------------------------------------------------------------------
!        Update global matrices from local matrices
!------------------------------------------------------------------------------
         CALL Condensate( DOFs*n, STIFF, FORCE, TimeForce )
         CALL DefaultUpdateEquations( STIFF, FORCE )
!------------------------------------------------------------------------------
      END DO     !  Bulk elements
      CALL DefaultFinishBulkAssembly()
      CALL Info( 'V2-F-Solver', 'Assembly done', Level=4 )

!------------------------------------------------------------------------------

      CALL DefaultFinishAssembly()
!------------------------------------------------------------------------------
!     Dirichlet boundary conditions
!------------------------------------------------------------------------------
      CALL DefaultDirichletBCs()
!------------------------------------------------------------------------------
      CALL Info( 'V2-F-Solver', 'Set boundaries done', Level=4 )
!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
      Norm = DefaultSolve()

      Solver % Variable % Values(1::3) = &
        MAX( Solver % Variable % Values(1::3), 1.0d-9 )

!------------------------------------------------------------------------------

      RelativeChange = Solver % Variable % NonlinChange

      WRITE( Message,* ) 'Result Norm   : ',Norm
      CALL Info( 'V2-F-Solver', Message, Level = 4 )
      WRITE( Message,* ) 'Relative Change : ',RelativeChange
      CALL Info( 'V2-F-Solver', Message, Level = 4 )

      IF ( RelativeChange < NewtonTol .OR. &
              iter > NewtonIter ) NewtonLinearization = .TRUE.

      IF ( Solver % Variable % NonlinConverged == 1 ) EXIT
!------------------------------------------------------------------------------
    END DO
!------------------------------------------------------------------------------

CONTAINS

!------------------------------------------------------------------------------
   SUBROUTINE LocalMatrix( MASS,STIFF,FORCE, &
            LOAD,UX,UY,UZ,Element,n,nd,Nodes )
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
!  REAL(KIND=dp) :: UX(:),UY(:),UZ(:)
!     INPUT: Nodal values of velocity components from previous iteration
!           used only if coefficient of the convection term (C1) is nonzero
!
!  LOGICAL :: Stabilize
!     INPUT: Should stabilzation be used ? Used only if coefficient of the
!            convection term (C1) is nonzero
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

     REAL(KIND=dp), DIMENSION(:)   :: FORCE,UX,UY,UZ
     REAL(KIND=dp), DIMENSION(:,:) :: MASS,STIFF,LOAD

     LOGICAL :: Stabilize

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

     REAL(KIND=dp) :: A(3,3),M(3,3)
     REAL(KIND=dp) :: LoadatIp(2),Rho,mu,Tmu,EffVisc
     INTEGER :: i,j,c,p,q,t,dim,N_Integ,NBasis

     REAL(KIND=dp) :: s,u,v,w, K,E,Eta,LV2,Lf,Strain(3,3),ProdTensor(3,3), &
                     Vorticity(3,3), nu, Cv, Prod,div

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     REAL(KIND=dp) :: C1,C2,CT,CL,Cnu,Cmu,TimeScale,LengthScale2,aparm(nd)
     REAL(KIND=dp) :: SecInv,X,Y,Z,Re_T
     REAL(KIND=dp) :: Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3),SqrtMetric

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     LOGICAL :: stat,Convection, Bubbles

!------------------------------------------------------------------------------

     dim = CoordinateSystemDimension()

     FORCE = 0.0_dp
     STIFF = 0.0_dp
     MASS  = 0.0_dp

     Bubbles = .TRUE.
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IF ( Bubbles ) THEN
        IntegStuff = GaussPoints( element, element % TYPE % GaussPoints2 )
     ELSE
        IntegStuff = GaussPoints( element )
     END IF

     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!    Now we start integrating
!------------------------------------------------------------------------------
     DO t=1,N_Integ
       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)
!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,detJ, &
               Basis,dBasisdx,Bubbles=.TRUE. )
!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
       s = detJ * S_Integ(t)
       IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
         X = SUM( Nodes % x(1:n)*Basis(1:n) )
         Y = SUM( Nodes % y(1:n)*Basis(1:n) )
         Z = SUM( nodes % z(1:n)*Basis(1:n) )
         CALL CoordinateSystemInfo(Metric,SqrtMetric,Symb,dSymb,X,Y,Z)

         s = s * SqrtMetric
       END IF

!      Velocity from previous iteration at the integration point
!------------------------------------------------------------------------------
       Velo = 0.0d0
       Velo(1) = SUM( UX(1:n)*Basis(1:n) )
       Velo(2) = SUM( UY(1:n)*Basis(1:n) )
       Velo(3) = SUM( UZ(1:n)*Basis(1:n) )

       dVelodx = 0.0d0
       DO i=1,dim
         dVelodx(1,i) = SUM( UX(1:n)*dBasisdx(1:n,i) )
         dVelodx(2,i) = SUM( UY(1:n)*dBasisdx(1:n,i) )
         dVelodx(3,i) = SUM( UZ(1:n)*dBasisdx(1:n,i) )
       END DO

       Strain = 0.5_dp * ( dVelodx + TRANSPOSE(dVelodx) )
       IF ( CurrentCoordinateSystem() == Cartesian ) THEN
          Secinv = 2*SUM(Strain * Strain)
       ELSE
          SecInv = SecondInvariant(Velo,dVelodx,Metric,Symb)/2
       END IF

!------------------------------------------------------------------------------
!      Coefficient of the convection and time derivative terms
!      at the integration point
!------------------------------------------------------------------------------
       K   = SUM( KinEne(1:n) * Basis(1:n) )
       E   = SUM( KinDis(1:n) * Basis(1:n) )
       LF  = SUM( F(1:n)  * Basis(1:n) )
       LV2 = SUM( V2(1:n) * Basis(1:n) )

       Cnu = SUM( V2FCnu(1:n) * Basis(1:n) )
       Cmu = SUM( KECMu(1:n)  * Basis(1:n) )
       CT  = SUM( V2FCT(1:n)  * Basis(1:n) )
       CL  = SUM( V2FCL(1:n)  * Basis(1:n) )
       C1  = SUM( V2FC1(1:n)  * Basis(1:n) )
       C2  = SUM( V2FC2(1:n)  * Basis(1:n) )

       rho = SUM( Basis(1:n) * Density(1:n) )
       mu  = SUM( Basis(1:n) * Viscosity(1:n) )

       TimeScale = MAX( K/E, CT*SQRT((mu/rho)/E))
       Lengthscale2 = CL**2 * MAX(K**3 /E**2, Cnu**2*SQRT((mu/rho)**3/E)) 

       Tmu = rho*Cmu*LV2*TimeScale
       EffVisc = mu + Tmu / SUM(V2FSigma(1:n)*Basis(1:n))

!      div = Strain(1,1) + Strain(2,2) + Strain(3,3)
!      ProdTensor = 2*Tmu*Strain
!      DO i=1,3
!        ProdTensor(i,i) = ProdTensor(i,i) - 2.0_dp*(Tmu*div+K)/3.0_dp
!      END DO
!      Prod = SUM(ProdTensor * dVelodx) / Rho
       Prod = Tmu * SecInv / Rho

!------------------------------------------------------------------------------
!      Loop Over basis functions of both unknowns and weights
!------------------------------------------------------------------------------
       DO p=1,nd
       DO q=1,nd
!------------------------------------------------------------------------------
!         The diffusive-convective equation without stabilization
!------------------------------------------------------------------------------
          M = 0.0d0
          A = 0.0d0

          M(1,1) = rho * Basis(q) * Basis(p)

          A(1,1) = A(1,1) + 6 * Rho / TimeScale  * Basis(q) * Basis(p)
          A(1,2) = A(1,2) - Rho * K * Basis(q) * Basis(p)

          A(2,1) = A(2,1) + (C1-6) / K / TimeScale * Basis(q) * Basis(p)
          A(2,2) = A(2,2) + Basis(q) * Basis(p)

          A(3,3) = Basis(q) * Basis(p)
          A(3,1) = A(3,1) - E / K**2  * Basis(q) * Basis(p)

          ! The diffusion term:
          !--------------------
          IF ( CurrentCoordinateSystem() == Cartesian ) THEN
             DO i=1,dim
               A(1,1) = A(1,1) + EffVisc*dBasisdx(q,i)*dBasisdx(p,i)
               A(2,2) = A(2,2) + LengthScale2*dBasisdx(q,i)*dBasisdx(p,i)
               A(2,3) = A(2,3) - 5*LengthScale2*dBasisdx(q,i)*dBasisdx(p,i)
             END DO
          ELSE
             DO i=1,dim
               DO j=1,dim
                  A(1,1) = A(1,1) + Metric(i,j) * EffVisc * &
                       dBasisdx(q,i) * dBasisdx(p,i)

                  A(2,2) = A(2,2) + Metric(i,j) * LengthScale2 * &
                       dBasisdx(q,i) * dBasisdx(p,i)

                  A(2,3) = A(2,3) - Metric(i,j) * 5*LengthScale2 * &
                       dBasisdx(q,i) * dBasisdx(p,i)
               END DO
             END DO
          END IF
          ! The convection term:
          !---------------------
          DO i=1,dim
            A(1,1) = A(1,1) + rho * Velo(i) * dBasisdx(q,i) * Basis(p)
          END DO

          DO i=1,3
             DO j=1,3
               STIFF(3*(p-1)+i,3*(q-1)+j) = STIFF(3*(p-1)+i,3*(q-1)+j)+s*A(i,j)
               MASS(3*(p-1)+i,3*(q-1)+j)  = MASS(3*(p-1)+i,3*(q-1)+j) +s*M(i,j)
             END DO
          END DO
       END DO
       END DO

       ! Load at the integration point:
       !-------------------------------
       LoadAtIP = 0.0_dp
       LoadAtIP(2) = ((C1-1)*2/3.0_dp+C2*Prod/E)/TimeScale

!------------------------------------------------------------------------------
        DO p=1,nd
           FORCE(3*(p-1)+2) = FORCE(3*(p-1)+2) + s*LoadAtIp(2)*Basis(p)
        END DO
      END DO
!------------------------------------------------------------------------------
   END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  END SUBROUTINE V2F
!------------------------------------------------------------------------------
