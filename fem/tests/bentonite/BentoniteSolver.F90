!/*****************************************************************************
! *
! *       ELMER, A Computational Fluid Dynamics Program.
! *
! *       Copyright 1st April 1995 - , CSC - IT Center for Science Ltd.,
! *                                    Finland.
! *
! ****************************************************************************/
!
!/*****************************************************************************
! *
! *****************************************************************************
! *
! *                    Author:       Juha Ruokolainen
! *
! *                 Address: CSC - IT Center for Science Ltd.
! *                        Keilaranta 14, P.O. BOX 405
! *                          02101 Espoo, Finland
! *                          Tel. +358 0 457 2723
! *                        Telefax: +358 0 457 2302
! *                      EMail: Juha.Ruokolainen@csc.fi
! *
! *                       Date: 04 Oct 2000
! *
! *                Modified by:
! *
! *       Date of modification:
! *
! ****************************************************************************/

!------------------------------------------------------------------------------
RECURSIVE SUBROUTINE BentoniteSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve the Bentonite equation!
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh, materials, BCs, etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear & nonlinear equation solver options
!
!  REAL(KIND=dp) :: dt,
!     INPUT: Timestep size for time dependent simulations
!
!  LOGICAL :: TransientSimulation
!     INPUT: Steady state or transient simulation
!
!******************************************************************************
   USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  TYPE(Nodes_t) :: ElementNodes

  INTEGER, POINTER :: NodeIndexes(:)

  LOGICAL :: AllocationsDone = .FALSE., Bubbles, stat

  INTEGER, POINTER :: BentonitePerm(:)
  REAL(KIND=dp), POINTER :: Bentonite(:)

  INTEGER :: iter, i, j, k, n, t, istat, eq, LocalNodes
  REAL(KIND=dp) :: Norm, PrevNorm, RelativeChange

  TYPE(ValueList_t), POINTER :: BC, SolverParams

  INTEGER :: NonlinearIter, DOFs, active
  REAL(KIND=dp) :: NonlinearTol,s

  REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), STIFF(:,:), &
        FORCE(:), SoundSpeed(:), Xi(:), Temp(:), Eta(:)

  SAVE STIFF, MASS, FORCE, ElementNodes, &
                 Xi, Temp, Eta, AllocationsDone

  TYPE(Mesh_t), POINTER :: Mesh

#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: at,at0,totat,st,totst,t1
#else
  REAL(KIND=dp) :: at,at0,totat,st,totst,t1,CPUTime,RealTime
#endif
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Get variables needed for solution
!------------------------------------------------------------------------------
  IF ( .NOT.ASSOCIATED( Solver % Matrix ) ) RETURN

  Mesh => GetMesh()

  Bentonite => Solver % Variable % Values
  DOFs = Solver % Variable % DOFs
  BentonitePerm => Solver % Variable % Perm

  LocalNodes = COUNT( BentonitePerm > 0 )
  IF ( LocalNodes <= 0 ) RETURN

  Norm = Solver % Variable % Norm
!------------------------------------------------------------------------------
! Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
  IF ( .NOT. AllocationsDone .OR. Mesh % Changed ) THEN
     N = Mesh % MaxElementNodes

     IF ( AllocationsDone ) THEN
        DEALLOCATE(                 &
             ElementNodes % x,      &
             ElementNodes % y,      &
             ElementNodes % z,      &
             FORCE, STIFF, MASS, Xi, Eta, Temp )
     END IF

     ALLOCATE(                   &
          ElementNodes % x( N ), &
          ElementNodes % y( N ), &
          ElementNodes % z( N ), &
          FORCE( DOFs*N ),       &
          STIFF(DOFs*N,DOFs*N ), &
          MASS( DOFs*N,DOFs*N ),  Xi(N), Temp(N), Eta(N), STAT=istat )

     IF ( istat /= 0 ) THEN
        CALL Fatal( 'BentoniteSolve', 'Memory allocation error.' )
     END IF

!
!    Compute initial value for Eta:
!    ------------------------------
     DO i=1,Mesh % NumberOFnodes
        IF ( BentonitePerm(i) > 0 ) THEN
           Xi(1)   = Bentonite( DOFs*(BentonitePerm(i)-1) + 1 )
           Temp(1) = Bentonite( DOFs*(BentonitePerm(i)-1) + 2 )
           Eta(1)  = Bentonite( DOFs*(BentonitePerm(i)-1) + 3 )
           CALL CompEta( Xi(1), Temp(1), Eta(1) )
           Bentonite( DOFs*(BentonitePerm(i)-1) + 3 ) = Eta(1)
        END IF
     END DO

     AllocationsDone = .TRUE.
  END IF
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Do some additional initialization, and go for it
!------------------------------------------------------------------------------
  SolverParams => GetSolverParams()

  NonlinearTol = GetConstReal( SolverParams, &
       'Nonlinear System Convergence Tolerance',stat )

  NonlinearIter = GetInteger( SolverParams, &
       'Nonlinear System Max Iterations', stat )

  IF ( .NOT.stat ) NonlinearIter = 1
!------------------------------------------------------------------------------
! Iterate over any nonlinearity of material or source
!------------------------------------------------------------------------------
  totat = 0.0d0
  totst = 0.0d0

  DO iter=1,NonlinearIter
!------------------------------------------------------------------------------
     at  = CPUTime()
     at0 = RealTime()

     CALL Info( 'BentoniteSolve', ' ', Level=4 )
     CALL Info( 'BentoniteSolve', '-------------------------------------', Level=4 )
     WRITE( Message, * ) 'Bentonite iteration', iter
     CALL Info( 'BentoniteSolve', Message, Level=4 )
     CALL Info( 'BentoniteSolve', '-------------------------------------', Level=4 )
     CALL Info( 'BentoniteSolve', ' ', Level=4 )
     CALL Info( 'BentoniteSolve', 'Starting Assembly', Level=4 )

!
!    Do the bulk assembly:
!    ---------------------
     active = GetNOFActive()
     CALL DefaultInitialize()
!------------------------------------------------------------------------------
     DO t=1,active
!------------------------------------------------------------------------------
        IF ( RealTime() - at0 > 1.0 ) THEN
          WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
           (Solver % NumberOfActiveElements-t) / &
              (1.0*Solver % NumberOfActiveElements)), ' % done'
          CALL Info( 'BentoniteSolve', Message, Level=5 )
                      
          at0 = RealTime()
        END IF
!------------------------------------------------------------------------------
        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes

        ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
        ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
        ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)


        ! Get elementwise nodal values of the current solution:
        !------------------------------------------------------
        Xi(1:n)   = Bentonite( DOFs * (BentonitePerm(NodeIndexes)-1) + 1 )
        Temp(1:n) = Bentonite( DOFs * (BentonitePerm(NodeIndexes)-1) + 2 )
        Eta(1:n)  = Bentonite( DOFs * (BentonitePerm(NodeIndexes)-1) + 3 )

        ! Get element local matrix and rhs vector:
        !-----------------------------------------
        CALL LocalMatrix( MASS, STIFF, FORCE, &
               Xi, Temp, Eta, Element, n, ElementNodes )

        IF ( TransientSimulation ) CALL Default1stOrderTime(MASS,STIFF,FORCE)
        CALL DefaultUpdateEquations( STIFF, FORCE )
!------------------------------------------------------------------------------
     END DO
!------------------------------------------------------------------------------

!
!    Neumann & Newton BCs:
!    ---------------------
!------------------------------------------------------------------------------
     DO t = 1,Mesh % NumberOfBoundaryElements
!------------------------------------------------------------------------------
        Element => GetBoundaryElement(t)
        IF ( .NOT. ActiveBoundaryElement() ) CYCLE

        BC => GetBC()
        IF ( .NOT.GetLogical( BC, 'Robin Condition', stat ) ) CYCLE
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes

        ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
        ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
        ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)

        ! Get elementwise nodal values of the current solution:
        !------------------------------------------------------
        Xi(1:n)   = Bentonite( DOFs * (BentonitePerm(NodeIndexes)-1) + 1 )
        Temp(1:n) = Bentonite( DOFs * (BentonitePerm(NodeIndexes)-1) + 2 )
        Eta(1:n)  = Bentonite( DOFs * (BentonitePerm(NodeIndexes)-1) + 3 )

        ! Get element local matrix and rhs vector:
        !-----------------------------------------
        CALL LocalMatrixBoundary(  STIFF, FORCE, &
             Xi, Temp, Eta, Element, n, ElementNodes )

        !Update global matrix and rhs vector from local matrix & vector:
        !---------------------------------------------------------------
        IF ( TransientSimulation ) THEN
           ! NOTE: This will replace STIFF and FORCE with the
           !      combined information...
           !--------------------------------------------------
           MASS = 0.0d0
           CALL Default1stOrderTime( MASS, STIFF, FORCE )
        END IF
        CALL DefaultUpdateEquations( STIFF, FORCE )
     END DO
!------------------------------------------------------------------------------

     CALL DefaultFinishAssembly()
     CALL DefaultDirichletBCs()
     
     ! Solve the system and we are done:
     !----------------------------------
     at = CPUTime() - at
     st = CPUTime()

     PrevNorm = Norm
     Norm = DefaultSolve()

     st = CPUTIme()-st
     totat = totat + at
     totst = totst + st
     WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Assembly: (s)', at, totat
     CALL Info( 'BentoniteSolve', Message, Level=4 )
     WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Solve:    (s)', st, totst
     CALL Info( 'BentoniteSolve', Message, Level=4 )

!------------------------------------------------------------------------------
     IF ( PrevNorm + Norm /= 0.0d0 ) THEN
        RelativeChange = 2*ABS(PrevNorm - Norm) / (PrevNorm + Norm)
     ELSE
        RelativeChange = 0.0d0
     END IF

     CALL Info( 'BentoniteSolve', ' ', Level=4 )
     WRITE( Message, * ) 'Result Norm    : ',Norm
     CALL Info( 'BentoniteSolve', Message, Level=4 )
     WRITE( Message, * ) 'Relative Change: ',RelativeChange
     CALL Info( 'BentoniteSolve', Message, Level=4 )

     IF ( RelativeChange < NonlinearTol ) EXIT
!------------------------------------------------------------------------------
  END DO ! of nonlinear iteration
!------------------------------------------------------------------------------

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE CompEta( Xi, Temp, Eta )
!------------------------------------------------------------------------------
    IMPLICIT NONE

    REAL(KIND=dp) :: Xi,Temp,Eta,nu_0,P_0,rho_1,M_v,R,k_u,mu_1,a,D_0,nD,c_vp,c_lp
    REAL(KIND=dp) :: Eta_0,T_0,L_0,rho_dry,c_eff,A1,A2,dx,x_0,Delta_eff,Xi_0,Xi_max

    nu_0 = 0.45d0
    P_0 = 1.013d5
    rho_1 = 998.0d0
    M_v = 0.018d0
    R = 8.314351d0
    k_u = 2.0d-21
    mu_1 = 1.0d-3
    a = 0.4d0
    D_0 = 0.216d-4
    nD = 1.8d0
    c_vp = 1.866d3
    c_lp = 4.182d3
    Eta_0 = 0.023d0
    Xi_0 = 0.54d0
    Xi_max = 0.999d0
    T_0 = 293.15d0
    L_0 = 2454.3d3 - (c_vp-c_lp)*T_0
    rho_dry = 1690.0d0

    c_eff = nu_0*Xi*rho_1*c_lp + rho_dry*(1.38*T_0 + 355.6);

    A1 = 0.57d0
    A2 = 1.28D0
    x_0 = 0.65d0
    dx = 0.1d0
    Delta_eff = A2 + (A1-A2)/(1+EXP((Xi-x_0)/dx))

    Eta = EXP( M_v/(R*Temp)*(L_0*(Temp-T_0)/T_0 + &
          (c_vp-c_lp)*Temp*LOG(Temp/T_0)) - a*(1/Xi**2 - 1/Xi_max**2)) * Eta_0
!------------------------------------------------------------------------------
  END SUBROUTINE CompEta
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  MASS, STIFF, FORCE, NodalXi, &
         NodalT, NodalEta, Element, n, Nodes )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: MASS(:,:),STIFF(:,:), FORCE(:)
    REAL(KIND=dp) :: NodalXi(:),NodalT(:),NodalEta(:)
    INTEGER :: n
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3)
    REAL(KIND=dp) :: detJ,U,V,W,S
    REAL(KIND=dp) :: LMASS(3,3), LSTIFF(3,3), LFORCE(3)
    LOGICAL :: Stat
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
    INTEGER :: i,j,k,l,p,q,t,dim, NBasis, CoordSys

    REAL(KIND=dp) :: x,y,z,SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3)
    REAL(KIND=dp) :: Grad, Base

    REAL(KIND=dp) :: a0,a1,a2,a3,a4,a5,b0,b1,b2,b3,c0,c1,c2,F, SMALL=1.0D-12
    REAL(KIND=dp) :: Xi,Temp,Eta,nu_0,P_0,rho_1,M_v,R,k_u,mu_1,a,D,c_vp,c_lp
    REAL(KIND=dp) :: Eta_0,T_0,L_0,rho_dry,c_eff, Ar1,Ar2,dx,x_0,Delta_eff,Xi_0, area
    REAL(KIND=dp) :: Xi_max, nD, D_0
!------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()
    CoordSys = CurrentCoordinateSystem()

    Metric = 0.0d0
    Metric(1,1) = 1.0d0
    Metric(2,2) = 1.0d0
    Metric(3,3) = 1.0d0

    MASS  = 0.0d0
    STIFF = 0.0d0
    FORCE = 0.0d0
!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------
    NBasis = n
    IntegStuff = GaussPoints( Element )
!------------------------------------------------------------------------------
    DO t=1,IntegStuff % n
       U = IntegStuff % u(t)
       V = IntegStuff % v(t)
       W = IntegStuff % w(t)
       S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, U, V, W, detJ, &
               Basis, dBasisdx )

       s = s * detJ
       IF ( CoordSys /= Cartesian ) THEN
          X = SUM( Nodes % x(1:n) * Basis(1:n) )
          Y = SUM( Nodes % y(1:n) * Basis(1:n) )
          Z = SUM( Nodes % z(1:n) * Basis(1:n) )
          CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,X,Y,Z )
          s = s * SqrtMetric
       END IF
!------------------------------------------------------------------------------
!      Current solution at the integration point
!------------------------------------------------------------------------------
       Xi = SUM( Basis(1:n) * NodalXi(1:n) )
       Temp = SUM( Basis(1:n) * NodalT(1:n) )
       Eta  = MAX( SUM( Basis(1:n) * NodalEta(1:n) ), SMALL )
!------------------------------------------------------------------------------
!      The source term and the coefficient of the time derivative and 
!      diffusion terms at the integration point
!------------------------------------------------------------------------------

       nu_0    = 0.45d0
       P_0     = 1.013d5
       rho_1   = 998.0d0
       M_v     = 0.018d0
       R       = 8.314351d0
       k_u     = 2.0d-21
       mu_1    = 1.0d-3
       a       = 0.4d0
       D_0     = 0.216d-4
       nD      = 1.8d0
       c_vp    = 1.866d3
       c_lp    = 4.182d3
       Eta_0   = 0.023d0
       Xi_0    = 0.54d0
       Xi_max  = 0.999d0
       T_0     = 293.15d0
       L_0     = 2454.3d3 - (c_vp-c_lp)*T_0
       rho_dry = 1690.0d0

       c_eff = nu_0*Xi*rho_1*c_lp + rho_dry*(1.38*T + 355.6);

       Ar1  = 0.57d0
       Ar2  = 1.28D0
       x_0  = 0.65d0
       dx   = 0.1d0
       Delta_eff = Ar2 + (Ar1-Ar2)/(1+EXP((Xi-x_0)/dx))

!
!      1st (Xi) equation coefficients:
!      -------------------------------
       a0 =  nu_0 * (rho_1-Eta*M_v/(R*Temp)*P_0)                        ! dXi  / dt
       a1 = -nu_0 * (1-Xi)*Eta*M_v/(R*Temp**2) * P_0                    ! dT   / dt
       a2 =  nu_0 * (1-Xi) * M_v/(R*Temp) * P_0                         ! dEta / dt
       a3 =  rho_1 * k_u / mu_1 * 2*a* rho_1 * R/M_v * Temp                          ! div(grad(Xi))
       a4 = -rho_1 * k_u / mu_1 * 2*a* rho_1 * R/M_v * (Xi/Xi_max)*(Xi_max-Xi)       ! div(grad(T))
       a5 =  M_v / (R*Temp) * P_0 * nu_0 * (1-Xi) * D_0*(Temp/273)**nD   ! div(grad(Eta))

!
!      2nd (T) equation coefficients:
!      ------------------------------
       b0 = -(L_0+(c_vp-c_lp)*Temp)*nu_0*rho_1                          ! dXi / dt
       b1 = c_eff                                                       ! dT  / dt
       b2 = -(L_0+(c_vp-c_lp)*Temp)*rho_1*k_u/mu_1*2*a*rho_1*R*Temp/M_v ! div(grad(Xi))
       b3 = Delta_eff                                                   ! div(grad(T))

!
!      3rd (Eta) equation coefficients (after newton lin.):
!      ----------------------------------------------------
       c0 = -2*a / Xi**3                                                ! Xi
       c1 = -(M_v*L_0/(R*Temp**2) + M_v*(c_vp-c_lp)/(R*Temp))           ! T
       c2 = 1.0d0 / Eta                                                 ! Eta

!
!      rhs for 3rd eq. (with newton lin. terms):
!      -----------------------------------------
       F = 0.0d0
       F = F - a*(1/Xi**2 - 1/Xi_max**2) - 2*a/Xi**2
       F = F + M_v*L_0/(R*T_0) - 2*M_v*L_0/(R*Temp) 
       F = F + M_v*(c_vp-c_lp)/R*LOG(Temp/T_0) - M_v*(c_vp-c_lp)/R
       F = F - LOG(Eta/Eta_0) + 1.0d0
       F = s * F

       DO p=1,NBasis
          DO q=1,NBasis
             LSTIFF = 0.0d0
             LMASS  = 0.0d0

             Base = s * Basis(p) * Basis(q)
             Grad = s * SUM( dBasisdx(q,1:dim) * dBasisdx(p,1:dim) )

             !
             ! 1st (Xi) equation:
             ! ------------------
             LMASS(1,1)  = LMASS(1,1)  + a0 * Base  ! dXi  / dt
             LMASS(1,2)  = LMASS(1,2)  + a1 * Base  ! dT   / dt
             LMASS(1,3)  = LMASS(1,3)  + a2 * Base  ! dEta / dt

             LSTIFF(1,1) = LSTIFF(1,1) + a3 * Grad  ! div(grad(Xi))
             LSTIFF(1,2) = LSTIFF(1,2) + a4 * Grad  ! div(grad(T))
             LSTIFF(1,3) = LSTIFF(1,3) + a5 * Grad  ! div(grad(Eta))

             !
             ! 2nd (T) equation:
             ! -----------------
             LMASS(2,1)  = LMASS(2,1)  + b0 * Base  ! dXi / dt
             LMASS(2,2)  = LMASS(2,2)  + b1 * Base  ! dT  / dt

             LSTIFF(2,1) = LSTIFF(2,1) + b2 * Grad  ! div(grad(Xi))
             LSTIFF(2,2) = LSTIFF(2,2) + b3 * Grad  ! div(grad(T))

             !
             ! 3rd (Eta) equation:
             ! -------------------
             LSTIFF(3,1) = LSTIFF(3,1) + c0 * Base  ! Xi
             LSTIFF(3,2) = LSTIFF(3,2) + c1 * Base  ! T
             LSTIFF(3,3) = LSTIFF(3,3) + c2 * Base  ! Eta

             !
             ! add nodal matrices to element matrices:
             ! ---------------------------------------
             DO i=1,DOFs
                k = DOFs*(p-1) + i
                DO j=1,DOFs
                   l = DOFs*(q-1) + j
                   MASS(k,l)  = MASS(k,l)  + LMASS(i,j)
                   STIFF(k,l) = STIFF(k,l) + LSTIFF(i,j)
                END DO
             END DO
          END DO

          !
          ! rhs (for 3rd eq):
          ! -----------------
          LFORCE = 0.0d0
          LFORCE(3) = LFORCE(3) + F * Basis(p)

          !
          ! add nodal vector to element vector:
          ! -----------------------------------
          DO i=1,DOFs
             k = DOFs*(p-1) + i
             FORCE(k) = FORCE(k) + LFORCE(i)
          END DO
       END DO
!------------------------------------------------------------------------------
    END DO
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBoundary(  STIFF, FORCE, &
         NodalXi, NodalT, NodalEta, Element, n, Nodes )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:),FORCE(:)
    REAL(KIND=dp) :: NodalXi(:), NodalT(:), NodalEta(:)
    INTEGER :: n
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: detJ,U,V,W,S,LFORCE(3), LSTIFF(3,3)
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),X,Y,Z
    REAL(KIND=dp) :: Normal(3)
    LOGICAL :: Stat
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
    INTEGER :: i,j,k,l,p,q,dim,CoordSys, NBasis, Integ

    REAL(KIND=dp) :: H_0, T_1, H, R_2, R_1, T_h, T_r, A, L_sat, T, Text
!------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()
    CoordSys = CurrentCoordinateSystem()

!------------------------------------------------------------------------------
!   Initialize stuff
!------------------------------------------------------------------------------
    STIFF = 0.0d0
    FORCE = 0.0d0

    T_1   = 286.58d0
    L_sat = 1.28d0

    H_0   = 1.9d0

    R_1   = 0.485d0
    R_2   = 1.14d0
    T_h   = 373.15d0
    T_r   = 285.15d0
    A = L_sat / ( ( R_2 * LOG(R_2 / R_1)) )

!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------
    NBasis = n
    IntegStuff = GaussPoints( Element )
!------------------------------------------------------------------------------
    DO Integ=1,IntegStuff % n
       U = IntegStuff % u(Integ)
       V = IntegStuff % v(Integ)
       W = IntegStuff % w(Integ)
       S = IntegStuff % s(Integ)
!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, U, V, W, detJ, &
               Basis, dBasisdx )

       s = s * detJ
       IF ( CoordSys /= Cartesian ) THEN
          x = SUM( Nodes % x(1:n) * Basis(1:n) )
          y = SUM( Nodes % y(1:n) * Basis(1:n) )
          z = SUM( Nodes % z(1:n) * Basis(1:n) )
          s = s * CoordinateSqrtMetric( x,y,z )
       END IF

!------------------------------------------------------------------------------
!      Current solution at the integration point
!------------------------------------------------------------------------------
       T = SUM( Basis(1:n) * NodalT(1:n) )

!      IF ( T <= T_1 ) THEN
          H = H_0
          Text = T_r
!      ELSE
!        H = A*(T_h-T)/(T-T_r)
!        Text = T_r
!      END IF

       DO p=1,NBasis
          DO q=1,NBasis
             LSTIFF = 0.0d0
             LSTIFF(2,2) = LSTIFF(2,2) + s * H * Basis(q) * Basis(p)

             !
             ! add nodal matrices to element matrices:
             ! ---------------------------------------
             DO i=1,DOFs
                k = DOFs*(p-1) + i
                DO j=1,DOFs
                   l = DOFs*(q-1) + j
                   STIFF(k,l) = STIFF(k,l) + LSTIFF(i,j)
                END DO
             END DO
          END DO

          LFORCE = 0.0d0
          LFORCE(2) = LFORCE(2) + s * H*Text * Basis(p)

          !
          ! add nodal vector to element vector:
          ! -----------------------------------
          DO i=1,DOFs
             k = DOFs*(p-1) + i
             FORCE(k) = FORCE(k) + LFORCE(i)
          END DO
       END DO
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBoundary
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE BentoniteSolver
!------------------------------------------------------------------------------
