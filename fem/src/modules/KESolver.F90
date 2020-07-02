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
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Matrix_t),POINTER  :: StiffMatrix
     INTEGER :: i,j,k,l,n,nd,nb,t,iter,istat,LocalNodes,DOFs

     TYPE(Nodes_t)   :: ElementNodes
     TYPE(Element_t),POINTER :: Element

     REAL(KIND=dp) :: RelativeChange,Norm,PrevNorm

     INTEGER, POINTER :: NodeIndexes(:)
     LOGICAL :: gotIt

     LOGICAL :: AllocationsDone = .FALSE., Bubbles


     TYPE(Variable_t), POINTER :: FlowSol, KE, KEkin, KEdis

     INTEGER, POINTER :: FlowPerm(:),KinPerm(:)

     INTEGER :: NSDOFs,NonlinearIter,NoActive
 
     REAL(KIND=dp), POINTER :: KEpsilon(:),&
         FlowSolution(:), ForceVector(:)

     REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), &
         STIFF(:,:),LayerThickness(:), &
         LOAD(:,:),FORCE(:),U(:),V(:),W(:), &
         Density(:),Viscosity(:),EffectiveVisc(:,:),Work(:),  &
         TurbulentViscosity(:),LocalDissipation(:), &
         LocalKinEnergy(:),LocalKEratio(:), &
         C0(:,:), SurfaceRoughness(:), TimeForce(:),LocalV2(:)
     
     TYPE(ValueList_t), POINTER :: BC

     SAVE U,V,W,MASS,STIFF,LOAD,FORCE, &
         ElementNodes,LayerThickness,Density,&
         AllocationsDone,Viscosity,LocalNodes,Work,TurbulentViscosity, &
         LocalDissipation,LocalKinEnergy,LocalKEratio, C0, &
         SurfaceRoughness, TimeForce, EffectiveVisc, LocalV2

     REAL(KIND=dp), POINTER :: SecInv(:)
     SAVE SecInv

     REAL(KIND=dp) :: at,at0,KMax, EMax, KVal, EVal

     TYPE(ValueList_t), POINTER :: SolverParams
     CHARACTER(*), PARAMETER :: Caller = 'KESolver'

     INTEGER :: KEComp, Ncomp
     TYPE(Variable_t), POINTER :: KEratio
     LOGICAL :: UseKEratio, LimitKE
     
!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
     IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

     SolverParams => GetSolverParams()
     
     KEComp = ListGetInteger( SolverParams,'KE component',GotIt)
     IF( KEComp == 0 ) THEN
       Ncomp = 2
     ELSE
       Ncomp = 1
     END IF
     
     KE => Solver % Variable
     IF ( ASSOCIATED( KE ) ) THEN
       DOFs     =  KE % DOFs
       KinPerm  => KE % Perm
       KEpsilon => KE % Values
     END IF

     KEkin => VariableGet( Solver % Mesh % Variables, 'Kinetic Energy' )
     KEdis => VariableGet( Solver % Mesh % Variables, 'Kinetic Dissipation' )
     KEratio => VariableGet( Solver % Mesh % Variables, 'KE ratio' )
     UseKEratio = ASSOCIATED( KEratio )      
     LimitKE = .NOT. UseKEratio
     
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
           LocalKEratio( N ), &
           LocalKinEnergy( N ),     &
           LocalDissipation( N ),&
           LocalV2(N), &
           MASS( 2*DOFs*N,2*DOFs*N ), &
           STIFF( 2*DOFs*N,2*DOFs*N ),LOAD( DOFs,N ), &
           FORCE( 2*DOFs*N ), TimeForce( 2*DOFs*N ), STAT=istat )

       IF ( istat /= 0 ) THEN
         CALL Fatal( Caller, 'Memory allocation error.' )
       END IF
                
       NULLIFY(SecInv)
       AllocationsDone = .TRUE.
     END IF

!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------
     NonlinearIter = ListGetInteger( SolverParams, &
         'Nonlinear System Max Iterations',GotIt )
     IF ( .NOT.GotIt ) NonlinearIter = 1

     Bubbles = GetString( SolverParams, &
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

       CALL Info( Caller, ' ', Level=5 )
       CALL Info( Caller,'-------------------------------------', Level=4 )
       IF( KEComp == 0 ) THEN
         CALL Info( Caller,'KEpsilon iteration: '//TRIM(I2S(iter)), Level=4 )
       ELSE IF( KEComp == 1 ) THEN
         CALL Info( Caller,'Kinetic energy iteration: '//TRIM(I2S(iter)), Level=4 )
       ELSE
         CALL Info( Caller,'Dissipation iteration: '//TRIM(I2S(iter)), Level=4 )
       END IF
       CALL Info( Caller,'-------------------------------------', Level=5 )
       CALL Info( Caller, ' ', Level=5 )
       CALL Info( Caller, 'Starting Assembly...', Level=5 )

       CALL DefaultInitialize()

       IF( UseKEratio ) THEN
         CALL ComputeKEratio()
       END IF
       
!------------------------------------------------------------------------------
!      Bulk elements
!------------------------------------------------------------------------------
       CALL StartAdvanceOutput(Caller,'Assembly' )
       
       NoActive = GetNOFActive()
       
       DO t=1,NoActive

!------------------------------------------------------------------------------

         CALL AdvanceOutput(t,NoActive)

!------------------------------------------------------------------------------
!        Check if this element belongs to a body where kinetic energy
!        should be calculated
!------------------------------------------------------------------------------
         Element => GetActiveElement(t)
         
!------------------------------------------------------------------------------
         CALL GetElementNodes( ElementNodes )
         n  = GetElementNOFNodes()
         nd = GetElementNOFDOFs()
         IF ( Bubbles ) nd=2*n
         nb = GetElementNOFBDOFs()
         NodeIndexes => Element % NodeIndexes

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
      CALL Info( Caller, 'Assembly done', Level=5 )

      IF(ListGetLogicalAnyBC(Model,'Epsilon Wall BC') .OR. &
          ListGetLogicalAnyBC(Model,'Noslip Wall BC') ) THEN
        DO t = 1, Solver % Mesh % NumberOfBoundaryElements
          Element => GetBoundaryElement(t)
          n = GetElementNOFNodes()

          BC => GetBC()
          IF( KEComp == 1 ) CYCLE
          IF ( ASSOCIATED( BC ) ) THEN
            IF ( GetLogical( BC, 'Epsilon Wall BC', gotIt ) .OR. &
                 GetLogical( BC, 'Noslip Wall BC',  gotIt ) ) THEN
              DO i=1,n
                j = KinPerm(Element % NodeIndexes(i))
                CALL ZeroRow( Solver % Matrix, Dofs*j )
                Solver % Matrix % RHS(Dofs*j) = 0.0_dp
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

              k = KinPerm(NodeIndexes(j))
              IF( KEcomp == 0 ) THEN             
                CALL UpdateDirichletDof( StiffMatrix, 2*k-1, Work(1) )
                CALL UpdateDirichletDof( StiffMatrix, 2*k, Work(2) )
              ELSE 
                CALL UpdateDirichletDof( StiffMatrix, k, Work(KEcomp) )
              END IF
            END DO
          END IF

          IF ( GetLogical( BC, 'Epsilon Wall BC', gotIt ) .OR. &
              GetLogical( BC, 'Noslip Wall BC',  gotIt ) ) THEN
            CALL EpsilonWall( Element, n, STIFF, FORCE )
          END IF
        END IF
      END DO
!------------------------------------------------------------------------------

      CALL DefaultFinishBoundaryAssembly()
      CALL DefaultFinishAssembly()
      
!------------------------------------------------------------------------------
!     Dirichlet boundary conditions
!------------------------------------------------------------------------------
      CALL DefaultDirichletBCs()
!------------------------------------------------------------------------------
      CALL Info( Caller, 'Set boundaries done', Level=5 )
!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
      PrevNorm = Norm

      Norm = DefaultSolve()
!------------------------------------------------------------------------------
!      Kinetic Energy Solution should be positive
!------------------------------------------------------------------------------

      IF( LimitKE ) THEN
        CALL LimitKEfields()
      END IF
        
!------------------------------------------------------------------------------
      RelativeChange = Solver % Variable % NonlinChange

      IF ( Solver % Variable % NonlinConverged == 1 ) EXIT
!------------------------------------------------------------------------------
    END DO
!------------------------------------------------------------------------------

CONTAINS

  SUBROUTINE ComputeKEratio()
    REAL(KIND=dp) :: Cmu, TmuMin, Tmu, Lmax, Lstar        
    TYPE(ValueList_t), POINTER :: Material

    Material => GetMaterial(GetActiveElement(1))

    Cmu = ListGetCReal( Material, 'KE Cmu', GotIt )
    TmuMin = GetCReal( Material, 'Minimum Turbulent Viscosity')
    Lmax = GetCReal( Material, 'Maximum Turbulent Mixing Length')

    n = Solver % Mesh % NumberOfNodes

    DO k=1,SIZE( KEkin % Values )
      Kval = MAX( KEkin % Values(k), 0.0_dp ) 
      Eval = MAX( KEdis % Values(k), 0.0_dp )

      Lstar = Lmax
      IF( Cmu * Kval**1.5 < Eval * Lmax ) THEN
        Lstar = Cmu * Kval**1.5 / Eval
      END IF
      
      Tmu  = MAX( Lstar * SQRT(Kval), TmuMin )
           
      KEratio % Values(k) = Cmu * Kval / Tmu
    END DO

  END SUBROUTINE ComputeKEratio

  
  SUBROUTINE LimitKEfields()
    REAL(KIND=dp) :: Clip, Cmu        
    TYPE(ValueList_t), POINTER :: Material

    Material => GetMaterial(GetActiveElement(1))

    Clip = GetConstReal( Material, 'KE Clip', GotIt )
    Cmu = ListGetCReal( Material, 'KE Cmu', GotIt )

    n = Solver % Mesh % NumberOfNodes
    Kmax = MAXVAL( KEkin % Values )
    Emax = MAXVAL( KEdis % Values ) 
    DO k=1,SIZE(KEkin % Values)
      Kval = KEkin % Values(k)
      Eval = KEdis % Values(k)
      
      IF ( KVal < Clip*Kmax ) Kval = Clip*KMax

      IF ( Eval < Clip*EMax ) THEN
        KVal = Clip*EMax
        Eval = MAX(Density(1)*Cmu*KVal**2/Viscosity(1),Clip*EMax)
      END IF

      IF( KEcomp /= 2 ) THEN
        KEkin % Values(k) = MAX( KVal, 1.0d-10 )
      END IF
      IF( KEcomp /= 1 ) THEN
        KEdis % Values(k) = MAX( EVal, 1.0d-10 )
      END IF
    END DO
    
  END SUBROUTINE LimitKEfields
    
   
  
!------------------------------------------------------------------------------
!>  Return element local matrices and RSH vector for KE model.
!------------------------------------------------------------------------------
   SUBROUTINE LocalMatrix( MASS,STIFF,FORCE, &
         LOAD, UX,UY,UZ,Element,n,nd,Nodes )
!------------------------------------------------------------------------------
     USE MaterialModels

     IMPLICIT NONE

     REAL(KIND=dp), DIMENSION(:,:) :: MASS  !< time derivative coefficient matrix
     REAL(KIND=dp), DIMENSION(:,:) :: STIFF !< rest of the equation coefficients
     REAL(KIND=dp), DIMENSION(:)   :: FORCE !< RHS vector
     REAL(KIND=dp), DIMENSION(:)   :: UX,UY,UZ !< Nodal values of velocity components from previous iteration.
     REAL(KIND=dp), DIMENSION(:,:) :: LOAD
     INTEGER :: n               !< Number of element nodes
     INTEGER :: nd              !< Number of dofs 
     TYPE(Nodes_t) :: Nodes     !< Element node coordinates
     TYPE(Element_t) :: Element !< Structure describing the element

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
!
     REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),detJ
     REAL(KIND=dp) :: Velo(3),dVelodx(3,3)
     REAL(KIND=dp) :: A(2,2),M(2,2),ProdK,ProdE
     
     INTEGER :: i,j,p,q,t,dim,NBasis
     REAL(KIND=dp) :: LoadatIp(2),Cmu,Rho,mu,Tmu,Effmu(2),TimeScale,LC1,LC2,LC3
     REAL(KIND=dp) :: s,u,v,w, K,E,Eta,Strain(3,3),&
         alpha,oldalpha,dalpha,err,ww,olderr,derr

     TYPE(GaussIntegrationPoints_t), TARGET :: IP

     REAL(KIND=dp) :: Pressure(n), Pref, Sound_speed_sq, Mach_number_sq
     REAL(KIND=dp) :: Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3), &
              SqrtMetric,Gravity(3),Pr,rho_g
     REAL(KIND=dp) :: C0(2),C1,CT,C2(2),SecInv,X,Y,Z,LV2,LCT,SigmaK,SigmaE,CvRat

     REAL(KIND=dp), POINTER :: gWork(:,:)

     LOGICAL :: stat,UseRNGModel,GotCvRat,GotGrav

     TYPE(ValueList_t), POINTER :: Material
     INTEGER :: body_id = -1
     CHARACTER(LEN=MAX_NAME_LEN) :: KEModel, V2FModel
     REAL(KIND=dp) :: V2FCp, EperK, TmuMin

     
!     REAL(KIND=dp) :: div,ProdTensor(3,3),nu

     SAVE Material, body_id, V2FCp, KEModel, V2FModel, UseRNGModel, &
         Gravity, GotGrav, Pr, LC3, LCT, Cmu, LC1, LC2, SigmaE, &
         SigmaK, cvRat, GotCvRat, Pref, TmuMin
     
     
!------------------------------------------------------------------------------

     dim = CoordinateSystemDimension()

     FORCE = 0.0D0
     STIFF = 0.0D0
     MASS  = 0.0D0

     NBasis = nd
     IF ( Bubbles ) NBasis = 2*n

     IF ( Element % BodyId /= body_id ) THEN
       body_id = Element % BodyId

       Material => GetMaterial()
              
       KEModel = GetString( Material, 'KE Model', GotIt )
       IF ( .NOT. GotIt ) KEModel = 'standard'
       
       V2FModel = GetString( Material, 'V2-F Model', GotIt )
     
       UseRNGModel = ( KEModel == 'rng' )
       
       gWork => ListGetConstRealArray( Model % Constants,'Gravity',GotGrav)
       IF ( GotGrav ) THEN
         Gravity = gWork(1:3,1)*gWork(4,1)
       ELSE
         Gravity    =  0.00_dp
         Gravity(2) = -9.81_dp
       END IF       
       ! Even historically gravity has not been used!
       GotGrav = .FALSE.     
       
       Pr = ListGetCReal( Material, 'Turbulent Prandtl Number', stat )
       IF ( .NOT. stat ) Pr = 0.85_dp
       
       LC3 = ListGetCReal( Material, 'Dissipation buoyancy coefficient', stat )
       
       IF ( KEModel == 'v2-f' ) THEN
         LCT = ListGetCReal( Material, 'V2-F CT', GotIt )
         IF ( .NOT. GotIt ) THEN
           LCT = 6.0_dp
           CALL ListAddConstReal( Material, 'V2-F CT', LCT )
         END IF

         V2FCp = ListGetCReal( Material, 'V2-F Cp', GotIt )
         IF ( .NOT. GotIt ) THEN
           V2FCp = 0.05_dp
           CALL ListAddConstReal( Material, 'V2-F Cp', V2FCp )
         END IF
       END IF

       Cmu = ListGetCReal( Material, 'KE Cmu', GotIt )
       IF ( .NOT. GotIt ) THEN
         SELECT CASE( KEModel )
         CASE( 'standard' )
           Cmu = 0.09_dp
         CASE( 'v2-f')
           Cmu = 0.22_dp
         CASE( 'rng' )
           Cmu = 0.0845_dp
         CASE DEFAULT
           CALL Fatal( Caller, 'Unknown K-Epsilon model' )
         END SELECT
         CALL ListAddConstReal( Material, 'KE Cmu', Cmu )
       END IF

       LC1 = ListGetCReal( Material, 'KE C1', GotIt )
       IF ( .NOT. GotIt ) THEN
         SELECT CASE( KEModel )
         CASE( 'standard' )
           LC1 = 1.44_dp
         CASE( 'v2-f' )
           LC1 = 1.4_dp
         CASE( 'rng' )
           LC1 = 1.42_dp
         CASE DEFAULT
           CALL Fatal( Caller, 'Unknown K-Epsilon model' )
         END SELECT
         CALL ListAddConstReal( Material, 'KE C1', LC1 )
       END IF

       LC2 = ListGetCReal( Material, 'KE C2', GotIt )
       IF ( .NOT. GotIt ) THEN
         SELECT CASE( KEModel )
         CASE( 'standard' )
           LC2 = 1.92_dp
         CASE( 'v2-f' )
           LC2 = 1.9_dp
         CASE( 'rng' )
           LC2 = 1.68_dp
         CASE DEFAULT
           CALL Fatal( Caller, 'Unknown K-Epsilon model' )
         END SELECT
         CALL ListAddConstReal( Material, 'KE C2', LC2 )
       END IF

       SigmaE = ListGetCReal( Material, 'KE SigmaE', GotIt )
       IF ( .NOT. GotIt ) THEN
         SELECT CASE( KEModel )
         CASE( 'standard','v2-f' )
           SigmaE = 1.3_dp
         CASE( 'rng' )
           SigmaE = 1.0_dp
         CASE DEFAULT
           CALL Fatal( Caller, 'Unknown K-Epsilon model' )
         END SELECT
         CALL ListAddConstReal( Material, 'KE SigmaE', SigmaE )
       END IF

       SigmaK = ListGetCReal( Material, 'KE SigmaK', GotIt )
       IF ( .NOT. GotIt ) THEN
         SigmaK = 1.0d0
         CALL ListAddConstReal( Material, 'KE SigmaK', SigmaK )
       END IF

       cvRat = GetCReal( Material, 'Specific Heat Ratio', GotCvRat )       
       Pref = GetCReal( Material, 'Reference Pressure', stat )

       TmuMin = GetCReal( Material, 'Minimum Turbulent Viscosity',GotIt)          
     END IF
            
     Viscosity(1:n) = GetReal( Material,'Viscosity' )
     IF( GotCvRat ) THEN
       CALL getScalarLocalSolution( Pressure, 'Pressure' )
     END IF

     ! This may include compressibility models etc. 
     !Density(1:n)   = GetReal( Material,'Density' )
     CALL ElementDensity( Density, n )

     CALL GetScalarLocalSolution( LocalV2, 'V2' )
     CALL GetScalarLocalSolution( LocalKinEnergy, 'Kinetic Energy' )
     CALL GetScalarLocalSolution( LocalDissipation, 'Kinetic Dissipation' )

     IF( UseKEratio ) THEN
       CALL GetScalarLocalSolution( LocalKEratio, 'KE ratio' )
     END IF
     
     CALL GetScalarLocalSolution( Ux, 'Velocity 1' )
     CALL GetScalarLocalSolution( Uy, 'Velocity 2' )
     CALL GetScalarLocalSolution( Uz, 'Velocity 3' )
       
   
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IF ( Bubbles ) THEN
       IP = GaussPoints( element, element % TYPE % GaussPoints2 )
     ELSE
       IP = GaussPoints( element )
     END IF

!------------------------------------------------------------------------------
!    Now we start integrating
!------------------------------------------------------------------------------
     DO t=1,IP % n
       u = IP % u(t)
       v = IP % v(t)
       w = IP % w(t)
!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,detJ, &
             Basis,dBasisdx,Bubbles=Bubbles )
!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
       s = detJ * IP % s(t)
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

       IF( UseKEratio ) THEN
         EperK = SUM( LocalKEratio(1:n) * Basis(1:n) )
       END IF
       
       IF( GotGrav ) THEN       
         rho_g = 0._dp
         DO i=1,dim
           rho_g = rho_g + SUM(Density(1:n) * dBasisdx(1:n,i)) * Gravity(i)
         END DO
       END IF
         
       mu  = SUM( Viscosity(1:n) * Basis(1:n) )
       rho = SUM( Density(1:n) * Basis(1:n) )

       Mach_number_sq = 0._dp
       IF( GotCvRat ) THEN
         Sound_speed_sq = (SUM(Basis(1:n)*Pressure(1:n))+Pref) * CvRat / rho
         IF ( Sound_speed_sq > 0.0_dp) Mach_number_sq = K/Sound_speed_sq
       END IF
         
       IF ( KEModel=='v2-f' ) THEN
         LV2 = SUM( LocalV2(1:n) * Basis(1:n) )
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
         IF( E > 0.0_dp ) THEN
           Tmu  = MAX( Rho * Cmu*K**2  / E, TmuMin )
         ELSE
           Tmu = TmuMin
         END IF
       END IF

       IF( GotGrav ) THEN              
         ProdK = Tmu * (SecInv-rho_g/(rho*Pr)) / rho
         ProdE = Tmu * (SecInv-LC3*rho_g/(rho*Pr)) / rho
       ELSE
         ProdK = Tmu * SecInv / rho
         ProdE = Tmu * SecInv / rho
       END IF
         
       Effmu(1) = mu + Tmu / SigmaK
       Effmu(2) = mu + Tmu / SigmaE

       C0(1) = Rho
       IF ( KEModel == 'v2-f' ) THEN
         C0(2) = Rho * LC2 / TimeScale
       ELSE
         IF( UseKEratio ) THEN
           C0(2) = Rho * LC2 * EperK
         ELSE
           C0(2) = Rho * LC2 * E / K
         END IF
       END IF

! moved to ---> rhs
!      C0(2) = C0(2) + Cmu*Rho*Eta**3*(1-Eta/4.38d0) / &
!                 (1.0d0 + 0.012d0*Eta**3) * E / K

       C1 = Rho
       CT = Rho
!------------------------------------------------------------------------------
!      Coefficient of the diffusion term &
!------------------------------------------------------------------------------
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

         C2(1) = Alpha * Effmu(1)
         C2(2) = Alpha * Effmu(2)
       ELSE
         C2(1) = Effmu(1)
         C2(2) = Effmu(2)        
       END IF
         
       ! Load at the integration point:
       !-------------------------------
       IF( GotCvRat ) THEN
         LoadAtIP(1) = rho*(ProdK-2*E*Mach_number_sq)
       ELSE
         LoadAtIp(1) = rho*ProdK
       END IF
       
       IF ( KEModel=='v2-f' ) THEN
         LoadAtIP(2) = Rho*LC1*ProdE/TimeScale
       ELSE
         IF( UseKEratio ) THEN
           LoadAtIP(2) = Rho*LC1*ProdE*EperK
         ELSE
           LoadAtIP(2) = Rho*LC1*ProdE*E/K
         END IF
       END IF
       
       IF ( UseRNGModel ) THEN
         IF( UseKEratio ) THEN
           LoadatIP(2) = LoadatIP(2) - Cmu*Rho*Eta**3*(1-Eta/4.38d0) / &
               (1.0d0 + 0.012d0*Eta**3) * E*EperK
         ELSE
           LoadatIP(2) = LoadatIP(2) - Cmu*Rho*Eta**3*(1-Eta/4.38d0) / &
               (1.0d0 + 0.012d0*Eta**3) * E**2 / K
         END IF
       END IF
         
         
!------------------------------------------------------------------------------
!       Loop over basis functions of both unknowns and weights
!------------------------------------------------------------------------------
       DO p=1,NBasis

         IF( KEcomp == 0 ) THEN
           FORCE(2*(p-1)+1) = FORCE(2*(p-1)+1)+s*LoadAtIp(1)*Basis(p)
           FORCE(2*(p-1)+2) = FORCE(2*(p-1)+2)+s*LoadAtIp(2)*Basis(p)
         ELSE 
           FORCE(p) = FORCE(p)+s*LoadAtIp(KEComp)*Basis(p)
         END IF
         
         DO q=1,NBasis
!------------------------------------------------------------------------------
!         The diffusive-convective equation without stabilization
!------------------------------------------------------------------------------
           M = 0.0d0
           A = 0.0d0
                      
           M(1,1) = CT * Basis(q) * Basis(p)
           M(2,2) = CT * Basis(q) * Basis(p)

           IF( UseKEratio ) THEN
             A(1,1) = C0(1) * EperK * Basis(q) * Basis(p) 
           ELSE
             A(1,2) = C0(1) * Basis(q) * Basis(p)
           END IF
             
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
           
           IF( KEComp == 0 ) THEN          
             DO i=1,2
               DO j=1,2
                 STIFF(2*(p-1)+i,2*(q-1)+j) = STIFF(2*(p-1)+i,2*(q-1)+j) + s*A(i,j)
                 MASS(2*(p-1)+i,2*(q-1)+j)  = MASS(2*(p-1)+i,2*(q-1)+j)  + s*M(i,j)
               END DO
             END DO
           ELSE
             STIFF(p,q) = STIFF(p,q) + s * A(KEComp,KEComp)
             MASS(p,q) = MASS(p,q) + s * M(KeComp,KEComp)
             IF( KEComp == 1 ) THEN
               FORCE(p) = FORCE(p) - s * A(1,2) * LocalDissipation(q)
             ELSE
               FORCE(p) = FORCE(p) - s * A(2,1) * LocalKinEnergy(q)
             END IF
           END IF
         END DO
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

     INTEGER :: i,j,p,q,t,np, dim     
     REAL(KIND=dp) :: s,u,v,w,E,K,Kder,X,Y,Z,Normal(3),mu,rho,Relax     
     TYPE(GaussIntegrationPoints_t), TARGET :: IP
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
     IP = GaussPoints( Element )

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

     DO t=1,IP % n
       u = IP % u(t)
       v = IP % v(t)
       w = IP % w(t)

       ! Basis function values & derivatives at the integration point:
       !--------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,detJ, BasisB )

       s = detJ * IP % s(t)
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
     INTEGER :: KEcomp 
!------------------------------------------------------------------------------
     SolverParams => GetSolverParams()

     KEComp = ListGetInteger( SolverParams,'KE component',Found ) 
     
     !CALL ListAddNewLogical( SolverParams,'Global Mass Matrix',.TRUE.) 
     
     str = GetString( SolverParams,'Variable', Found )

     IF( KEComp ==  0 ) THEN
       IF ( .NOT. Found ) str = 'K-eps'
       IF ( INDEX( str, '[' ) <= 0 ) THEN
         CALL ListAddString( SolverParams, 'Variable', &
             TRIM(str) // '[Kinetic Energy:1 Kinetic Dissipation:1]' )
       END IF
     ELSE IF( KEcomp == 1 ) THEN
       CALL ListAddString( SolverParams, 'Variable','Kinetic Energy' )
     ELSE IF( KEcomp == 2 ) THEN
       CALL ListAddString( SolverParams, 'Variable','Kinetic Dissipation' )
     ELSE
       CALL Fatal('KESolver_init','Invalid value for "KE Component": '//TRIM(I2S(KEcomp)))
     END IF
     
     IF( GetLogical( SolverParams,'USE KE ratio',Found) ) THEN
       CALL ListAddString( SolverParams,&
           NextFreeKeyword('Exported Variable',SolverParams),'KE ratio')
     END IF
 
!------------------------------------------------------------------------------
   END SUBROUTINE KESolver_Init
!------------------------------------------------------------------------------

