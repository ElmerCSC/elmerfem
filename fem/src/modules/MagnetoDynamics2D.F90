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
! *  Module for solving magnetic vector potential in cartesian and
! *  cylindrically symmetric 2D case. In both cases the vector potential
! *  is reduced to a single component. 
! *
! *  Authors: Juha Ruokolainen, Mika Malinen, Peter R�back
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
SUBROUTINE MagnetoDynamics2D_Init( Model,Solver,dt,TransientSimulation )
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

  Params => GetSolverParams()
  CALL ListAddInteger( Params, 'Variable Dofs',1 )
  IF( .NOT. ListCheckPresent(  Params,'Variable') ) THEN
    CALL ListAddString( Params,'Variable','Potential')
  END IF

  IF(.NOT. ListCheckPresent( Params,'Apply Mortar BCs') ) THEN
    CALL ListAddLogical( Params,'Apply Mortar BCs',.TRUE.)
  END IF
  
!------------------------------------------------------------------------------
END SUBROUTINE MagnetoDynamics2D_Init
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solver the magnetic vector potential in cartesian 2D case.
!> The solver may take into account rotating boundary conditions.
!> Also optionally compute moments and inertia. 
!------------------------------------------------------------------------------
SUBROUTINE MagnetoDynamics2D( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE CircuitUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: TransientSimulation !< Steady state or transient simulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  LOGICAL :: AllocationsDone = .FALSE., Found
  TYPE(Element_t), POINTER :: Element

  REAL(KIND=dp) :: Norm
  INTEGER :: i,j,k,n, nb, nd, t, istat, Active, NonlinIter, iter

  TYPE(ValueList_t), POINTER :: BC
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:), &
               NEWX(:), NEWY(:), POT(:)

  TYPE(Mesh_t),   POINTER :: Mesh

  LOGICAL :: NewtonRaphson = .FALSE., CSymmetry
  INTEGER :: CoupledIter
  TYPE(Variable_t), POINTER :: IterV, CoordVar

  TYPE(Matrix_t),POINTER::CM

!------------------------------------------------------------------------------

  CALL Info( 'MagnetoDynamics2D','------------------------------------------------', Level=4 )
  CALL Info( 'MagnetoDynamics2D', 'Solving equation for magnetic vector potential', Level=4 )
  CALL Info( 'MagnetoDynamics2D','------------------------------------------------', Level=4 )

  CALL DefaultStart()
  
  ! Allocate some permanent storage, this is done first time only:
  ! --------------------------------------------------------------
  NULLIFY(BC)
  Mesh => GetMesh()

  IF(GetCoupledIter()>1) NewtonRaphson = .TRUE.

  NonlinIter = GetInteger(GetSolverParams(), &
           'Nonlinear System Max Iterations',Found)
  IF(.NOT.Found) NonlinIter = 1

  CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )

  DO iter = 1,NonlinIter
    IF(Iter > 1) NewtonRaphson=.TRUE.
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
!$omp parallel do private(Element, n, nd, BC,Found)
    DO t=1,active
      Element => GetBoundaryElement(t)
      BC=>GetBC( Element )
      IF(.NOT.ASSOCIATED(BC)) CYCLE

      IF(GetLogical(BC,'Infinity BC',Found)) THEN
         n  = GetElementNOFNodes( Element )
         nd = GetElementNOFDOFs( Element )
         CALL LocalMatrixBC(Element, n, nd)
      ELSE IF(GetLogical(BC,'Air Gap',Found)) THEN
         n  = GetElementNOFNodes( Element )
         nd = GetElementNOFDOFs( Element )
         CALL LocalMatrixAirGapBC(Element, BC, n, nd)
      END IF
    END DO
!$omp end parallel do

    CALL DefaultFinishAssembly()

    CALL SetMagneticFluxDensityBC()
    CALL DefaultDirichletBCs()
    Norm = DefaultSolve()
 
    IF( Solver % Variable % NonlinConverged > 0 ) EXIT
  END DO

  ! For cylindrical symmetry the model lumping has not been implemented
  IF( .NOT. CSymmetry ) THEN
    CALL CalculateLumped(Model % NumberOfBodyForces)
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
  
CONTAINS

!------------------------------------------------------------------------------
 SUBROUTINE CalculateLumped(nbf)
!------------------------------------------------------------------------------
   INTEGER::nbf
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: torq,TorqArea,a(nbf),u(nbf),IMoment,IA, &
       rinner,router,ctorq
   INTEGER :: i,bfid,n,nd
   LOGICAL :: Found
   TYPE(ValueList_t),POINTER::Params
!------------------------------------------------------------------------------

   CALL Info('MagnetoDynamics2D','Calculating lumped parameters',Level=8)
   
   U=0._dp
   a=0._dp
   torq=0._dp
   TorqArea=0._dp
   IMoment=0._dp
   IA=0
   DO i=1,GetNOFActive()
     Element => GetActiveElement(i)
     nd = GetElementNOFDOFs(Element)
     n  = GetElementNOFNodes(Element)

     CALL Torque(Torq,TorqArea,Element,n,nd)

     Params=>GetBodyForce(Element)
     IF(ASSOCIATED(Params)) THEN
       bfid=GetBodyForceId(Element)
       IF(GetLogical(Params,'Calculate Potential',Found)) &
         CALL Potential(u(bfid),a(bfid),Element,n,nd)
     END IF

     Params=>GetBodyParams(Element)
     IF(ASSOCIATED(Params)) THEN
       IF(GetLogical(Params,'Calculate Inertial Moment',Found)) &
         CALL InertialMoment(IMoment,IA,Element,n,nd)
     END IF
   END DO

   DO i=1,nbf
     a(i) = ParallelReduction(a(i))
     u(i) = ParallelReduction(u(i))
   END DO
   IMoment = ParallelReduction(IMoment)
   IA = ParallelReduction(IA)

   Torq = ParallelReduction(Torq)
   WRITE(Message,'(A,ES15.4)') 'Air gap initial torque:', Torq
   CALL Info('MagnetoDynamics2D',Message,Level=8)

   TorqArea = ParallelReduction(TorqArea)
   rinner = ListGetCRealAnyBody( Model,'r inner',Found )
   router = ListGetCRealAnyBody( Model,'r outer',Found )
   IF (TorqArea /= 0) THEN
      Ctorq = PI*(router**2-rinner**2) / TorqArea
   ELSE
      Ctorq = 0.0_dp
   END IF
   WRITE(Message,'(A,ES15.4)') 'Air gap correction:', cTorq
   CALL Info('MagnetoDynamics2D',Message,Level=8)
   Torq = Ctorq * Torq
   
   DO i=1,nbf
     IF(a(i)>0) THEN
       CALL ListAddConstReal(Model % Simulation,'res: Potential / bodyforce ' &
                     //TRIM(i2s(i)),u(i)/a(i))
       CALL ListAddConstReal(Model % Simulation,'res: area / bodyforce ' &
                     //TRIM(i2s(i)),a(i))
     END IF
   END DO
   CALL ListAddConstReal(Model % Simulation,'res: air gap torque', Torq)
   CALL ListAddConstReal(Model % Simulation,'res: inertial volume', IA)
   CALL ListAddConstReal(Model % Simulation,'res: inertial moment', IMoment)
   

   WRITE(Message,'(A,ES15.4)') 'Air gap torque:', Torq
   CALL Info('MagnetoDynamics2D',Message,Level=7)
   WRITE(Message,'(A,ES15.4)') 'Inertial volume:', IA
   CALL Info('MagnetoDynamics2D',Message,Level=7)
   WRITE(Message,'(A,ES15.4)') 'Inertial moment:', Imoment
   CALL Info('MagnetoDynamics2D',Message,Level=7)
   
!------------------------------------------------------------------------------
 END SUBROUTINE CalculateLumped
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE InertialMoment(U,A,Element,n,nd)
!------------------------------------------------------------------------------
    INTEGER :: n,nd
    REAL(KIND=dp)::U,a
    TYPE(Element_t)::Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd), DetJ,x,y,r,Density(n)
    INTEGER :: t
    LOGICAL :: stat,Found
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
	!$OMP THREADPRIVATE(Nodes)

    Density(1:n) = GetReal(GetMaterial(),'Density',Found,Element)
    IF(.NOT.Found) RETURN

    CALL GetElementNodes( Nodes, Element )
  
    !Numerical integration:
    !----------------------
    IP = GaussPoints(Element)
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                IP % W(t), detJ, Basis )

      x = SUM(Nodes % x(1:nd)*Basis(1:nd))
      y = SUM(Nodes % y(1:nd)*Basis(1:nd))
      r = SQRT(x**2+y**2)
      A = A + IP % s(t)*detJ
      U = U + IP % s(t)*detJ*R*SUM(Density(1:n)*Basis(1:n))
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE InertialMoment
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE Torque(U,Area,Element,n,nd)
!------------------------------------------------------------------------------
    INTEGER :: n,nd
    REAL(KIND=dp)::U,Area
    TYPE(Element_t)::Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: dBasisdx(nd,3),Basis(nd), DetJ, &
             POT(nd),x,y,r,r0,r1,Br,Bp,Bx,By,B(3,nd)
    INTEGER :: t
    LOGICAL :: stat, Found
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
	!$OMP THREADPRIVATE(Nodes)

    r0 = GetCReal(GetBodyParams(),'r inner',Found)
    r1 = GetCReal(GetBodyParams(),'r outer',Found)
    IF (.NOT.Found) RETURN

    CALL GetElementNodes( Nodes, Element )

    x = SUM(Nodes % x(1:n))/n
    y = SUM(Nodes % y(1:n))/n
    r = SQRT(x**2+y**2)
    IF (r<r0.OR.r>r1) RETURN

    CALL GetLocalSolution(POT, UElement=Element)
  
    !Numerical integration:
    !----------------------
    IP = GaussPoints(Element)
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                IP % W(t), detJ, Basis, dBasisdx )

      x = SUM(Nodes % x(1:nd)*Basis(1:nd))
      y = SUM(Nodes % y(1:nd)*Basis(1:nd))
      r = SQRT(x**2+y**2)

      Bx =  SUM(POT*dBasisdx(:,2))
      By = -SUM(POT*dBasisdx(:,1))
      Br =  x/r*Bx + y/r*By
      Bp = -y/r*Bx + x/r*By
      U = U + IP % s(t)*detJ*r*Br*Bp/(PI*4.0d-7*(r1-r0))
      Area = Area + IP % s(t)*detJ
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE Torque
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE Potential( U, A, Element,n,nd)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: U,A
    INTEGER :: n, nd
    TYPE(Element_t) :: Element

    REAL(KIND=dp) :: Basis(nd), DetJ,POT(nd),pPOT(nd),dPOT(nd)
    INTEGER :: t
    LOGICAL :: stat
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
	!$OMP THREADPRIVATE(Nodes)

    CALL GetElementNodes( Nodes )

    CALL GetLocalSolution(POT,UElement=Element)
    CALL GetLocalSolution(pPOT,tstep=-1,UElement=Element)
    IF(Solver % Order<2.OR.GetTimeStep()<=2) THEN 
      dPot = (POT - pPOT)/dt
    ELSE
      dPot = 1.5_dp*POT - 2*pPOT
      CALL GetLocalSolution(pPOT,tstep=-2,UElement=Element)
      dPot = (dPOT + 0.5_dp*pPOT)/dt
    END IF


    !Numerical integration:
    !----------------------
    IP = GaussPoints(Element)
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                IP % W(t), detJ, Basis )
      A = A + IP % s(t) * detJ
      U = U + IP % s(t) * detJ * SUM(dPot(1:nd)*Basis(1:nd))
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE Potential
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(Element, n, nd)
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    INTEGER :: i,p,q,t,siz

    TYPE(GaussIntegrationPoints_t) :: IP

    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), &
      LOAD(nd),R(2,2,n),C(n), mu,muder,Babs,POT(nd), &
        JAC(nd,nd),Agrad(3),C_ip,M(2,n),M_ip(2),x

    LOGICAL :: Cubic, HBcurve, Found, Stat

    REAL(KIND=dp), POINTER :: Bval(:), Hval(:), Cval(:)
    REAL(KIND=dp), POINTER :: CubicCoeff(:) => NULL(), HB(:,:) => NULL()
    TYPE(ValueListEntry_t), POINTER :: Lst
    TYPE(ValueList_t), POINTER :: Material, BodyForce

    TYPE(Nodes_t), SAVE :: Nodes
!$omp threadprivate(Nodes, CubicCoeff, HB)
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilType
    LOGICAL :: CoilBody    
    TYPE(ValueList_t), POINTER :: CompParams

    REAL(KIND=dp) :: Bt(nd,2), Ht(nd,2)
    REAL(KIND=dp) :: nu_tensor(2,2)
    REAL(KIND=dp) :: B_ip(2), Alocal
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes,Element )
    STIFF = 0._dp
    JAC  = 0._dp
    FORCE = 0._dp
    IF(TransientSimulation) MASS = 0._dp

    Material => GetMaterial(Element)

    CALL GetConstRealArray( Material, HB, 'H-B curve', HBCurve )
    siz = 0
    Cval => NULL()
    IF ( HBCurve ) THEN
      siz = SIZE(HB,1)
      IF(siz>1) THEN
        Bval=>HB(:,1)
        Hval=>HB(:,2)
        Cubic = GetLogical( Material, 'Cubic spline for H-B curve',Found)
        IF (Cubic.AND..NOT.ASSOCIATED(CubicCoeff)) THEN
          ALLOCATE(CubicCoeff(siz))
          CALL CubicSpline(siz,Bval,Hval,CubicCoeff)
        END IF
        Cval=>CubicCoeff
        HBCurve = .TRUE.
      END IF
    END IF

    IF(siz<=1) THEN
      Lst => ListFind(Material,'H-B Curve',HBcurve)
      IF(HBcurve) THEN
        Cval => Lst % CubicCoeff
        Bval => Lst % TValues
        Hval => Lst % FValues(1,1,:)
      END IF
    END IF

    IF(HBcurve) THEN
      CALL GetLocalSolution(POT,UElement=Element,USolver=Solver)
      IF (.NOT. ASSOCIATED(Bval) ) CALL Fatal ('mgdyn2D','bval not associated')
      IF (.NOT. ASSOCIATED(Hval) ) CALL Fatal ('mgdyn2D','hval not associated')
    ELSE
      CALL GetReluctivity(Material,R,n,Element)
    END IF

    C = GetReal( Material, 'Electric Conductivity', Found, Element)

    M(1,:) = GetReal( Material, 'Magnetization 1', Found, Element)
    M(2,:) = GetReal( Material, 'Magnetization 2', Found, Element)

    Load = 0.0d0
    BodyForce => GetBodyForce(Element)
    IF ( ASSOCIATED(BodyForce) ) &
       Load(1:n) = GetReal(BodyForce, 'Current Density', Found, Element)

    !Numerical integration:
    !----------------------
    IP = GaussPoints(Element)
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
      LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )

      nu_tensor = 0.0_dp
      IF (HBcurve) THEN
        Agrad = 0.0_dp
        Agrad = MATMUL( POT,dBasisdx )
        Alocal = SUM( POT(1:n) * Basis(1:n) )
        ! Sign?
        ! -----
        B_ip(1) = -Agrad(2) 
        B_ip(2) = Agrad(1)
        IF( CSymmetry ) B_ip(2) = B_ip(2) + Alocal/x
        ! -----
        Babs = MAX( SQRT(SUM(B_ip**2)), 1.d-8 )
        mu = InterpolateCurve(Bval,Hval,Babs,CubicCoeff=Cval)/Babs
        muder = (DerivateCurve(Bval,Hval,Babs,CubicCoeff=Cval)-mu)/Babs
        nu_tensor(1,1) = mu ! Mu is really nu!!! too lazy to correct now...
        nu_tensor(2,2) = mu
      ELSE
        muder=0._dp
        DO p=1,2
          DO q=1,2
            nu_tensor(p,q) = SUM(Basis(1:n) * R(p,q,1:n))
          END DO
        END DO
     END IF

      CoilBody = .FALSE.
      CompParams => GetComponentParams( Element )
      CoilType = ''
      IF (ASSOCIATED(CompParams)) THEN
        CoilType = GetString(CompParams, 'Coil Type', Found)
        IF (Found) THEN
          SELECT CASE (CoilType)
          CASE ('stranded')
             CoilBody = .TRUE.
          CASE ('massive')
             CoilBody = .TRUE.
          CASE ('foil winding')
             CoilBody = .TRUE.
    !         CALL GetElementRotM(Element, RotM, n)
          CASE DEFAULT
             CALL Fatal ('MagnetoDynamics2DHarmonic', 'Non existent Coil Type Chosen!')
          END SELECT
        END IF
      END IF

      C_ip = SUM( Basis(1:n) * C(1:n) )
      M_ip = MATMUL( M,Basis(1:n) )


      ! Finally, the elemental matrix & vector:
      !----------------------------------------
      IF(TransientSimulation .AND. C_ip/=0._dp) THEN
        DO p=1,nd
          DO q=1,nd
            IF(CoilType /= 'stranded') MASS(p,q) = MASS(p,q) + IP % s(t) * detJ * C_ip * Basis(q)*Basis(p)
          END DO
        END DO
      END IF

      ! Is the sign correct?
      !---------------------
      Bt(:,1) = -dbasisdx(:,2)
      Bt(:,2) =  dbasisdx(:,1)
      IF ( CSymmetry ) Bt(:,2) = Bt(:,2) + Basis(:)/x

      DO p = 1,nd
        Ht(p,:) = MATMUL(nu_tensor, Bt(p,:))
      END DO

      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + IP % s(t) * DetJ * &
             MATMUL(Ht, TRANSPOSE(Bt))

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
              muder/babs*SUM(B_ip(:)*Bt(q,:))*SUM(B_ip*Bt(p,:))
          END DO
        END DO
      END IF

      FORCE(1:nd) = FORCE(1:nd) + IP % s(t) * DetJ * (LoadAtip * Basis(1:nd) + &
           (M_ip(1)*dBasisdx(1:nd,2)-M_ip(2)*dBasisdx(1:nd,1)))
    END DO

    IF (HBcurve .AND. NewtonRaphson) THEN
      STIFF = STIFF + JAC
      FORCE = FORCE + MATMUL(JAC,POT)
    END IF

    IF(TransientSimulation) THEN
      CALL Default1stOrderTime( MASS, STIFF, FORCE,UElement=Element, USolver=Solver )
    END IF
    CALL DefaultUpdateEquations( STIFF, FORCE,UElement=Element, USolver=Solver)

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC(Element, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    LOGICAL :: Stat
    INTEGER :: i,p,q,t
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), R(2,2,n), R_ip, &
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
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixAirGapBC(Element, BC, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    LOGICAL :: Stat, Found
    INTEGER :: i,p,q,t
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), R(n), R_ip, &
            Inf_ip,Coord(3),Normal(3),mu,u,v, AirGapLength(nd), &
            AirGapMu(nd), AirGapL

    TYPE(ValueList_t), POINTER :: Material
    TYPE(ValueList_t), POINTER :: BC

    TYPE(Element_t), POINTER :: Parent
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !$OMP THREADPRIVATE(Nodes)
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes, Element )
    STIFF = 0._dp
    FORCE = 0._dp

    AirGapLength=GetConstReal( BC, 'Air Gap Length', Found)
    if (.not. Found) CALL FATAL('LocalMatrixAirGapBC', 'Air Gap Length not found!')

    AirGapMu=GetConstReal( BC, 'Air Gap Relative Permeability', Found)
    if (.not. Found) AirGapMu=1d0

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
    REAL(KIND=dp) :: Bx(Solver % Mesh % MaxElementNodes), &
                      By(Solver % Mesh % MaxElementNodes)
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
        IF ( ListCheckPresent( BC, 'Magnetic Flux Density 1') .OR. &
             ListCheckPresent( BC, 'Magnetic Flux Density 2')      &
            ) THEN
          Bx = 0._dp
          By = 0._dp

          Bx(1:n) = GetReal(BC, 'Magnetic Flux Density 1', Found)
          IF (.NOT. Found) Bx = 0._dp
          By(1:n) = GetReal(BC, 'Magnetic Flux Density 2', Found)
          IF (.NOT. Found) By = 0._dp
          DO j = 1,n
            k = Element % NodeIndexes(j)
            x = Mesh % Nodes % x(k)
            y = Mesh % Nodes % y(k)
            k = Perm(k)
            !b(k) = y * Bx(j) - x * By(j)

            CALL UpdateDirichletDof( A, k, y * Bx(j) - x * By(j) )

            !CALL ZeroRow(A, k)
            !CALL AddToMatrixElement(A, k, k, 1._dp)
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
SUBROUTINE MagnetoDynamics2DHarmonic_Init0( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: TransientSimulation !< Steady state or transient simulation
!------------------------------------------------------------------------------
  IF( .NOT.ListCheckPresent( Solver % Values, 'Apply Mortar BCs') ) &
    CALL ListAddLogical( Solver % Values, 'Apply Mortar BCs', .TRUE.)

  IF( .NOT.ListCheckPresent( Solver % Values, 'Linear System Complex') ) &
    CALL ListAddLogical( Solver % Values, 'Linear System Complex', .TRUE.)
!------------------------------------------------------------------------------
END SUBROUTINE MagnetoDynamics2DHarmonic_Init0
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE MagnetoDynamics2DHarmonic_Init( Model,Solver,dt,TransientSimulation )
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

  Params => GetSolverParams(Solver)
  CALL ListAddInteger( Params, 'Variable Dofs',2 )
  IF( .NOT. ListCheckPresent(  Params,'Variable') ) THEN
    CALL ListAddString( Params,'Variable',&
        'Potential[Potential re:1 Potential im:1]')
  END IF

  IF(.NOT. ListCheckPresent( Params,'Apply Mortar BCs') ) THEN
    CALL ListAddLogical( Params,'Apply Mortar BCs',.TRUE.)
  END IF
  
  IF(.NOT. ListCheckPresent( Params,'Linear System Complex') ) THEN
    CALL ListAddLogical( Params,'Linear System Complex',.TRUE.)
  END IF

!------------------------------------------------------------------------------
END SUBROUTINE MagnetoDynamics2DHarmonic_Init
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solver the magnetic vector potential in cartesian 2D & complex case.
!> The solver may take into account rotating boundary conditions.
!> Also optionally compute moments and inertia. 
!------------------------------------------------------------------------------
SUBROUTINE MagnetoDynamics2DHarmonic( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE CircuitUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: TransientSimulation !< Steady state or transient simulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  LOGICAL :: AllocationsDone = .FALSE., Found
  TYPE(Element_t),POINTER :: Element

  REAL(KIND=dp) :: Norm
  INTEGER :: i,j,k,ip,jp,n, nb, nd, t, istat, Active, iter, NonlinIter

  TYPE(ValueList_t), POINTER :: BC
  TYPE(Mesh_t),   POINTER :: Mesh
  COMPLEX(KIND=dp), PARAMETER :: im=(0._dp,1._dp)

  LOGICAL, SAVE :: NewtonRaphson = .FALSE., CSymmetry
  INTEGER :: CoupledIter
  TYPE(Variable_t), POINTER :: IterV, CoordVar

  TYPE(Matrix_t),POINTER::CM

!------------------------------------------------------------------------------


  CALL Info( 'MagnetoDynamics2DHarmonic',&
      '------------------------------------------------', Level=4 )
  CALL Info( 'MagnetoDynamics2DHarmonic', &
      'Solving equation for magnetic vector potential', Level=4 )
  CALL Info( 'MagnetoDynamics2DHarmonic',&
      '------------------------------------------------', Level=4 )

  CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )

  ! Allocate some permanent storage, this is done first time only:
  ! --------------------------------------------------------------
  Mesh => GetMesh()
  NULLIFY(BC)

  IF(GetCoupledIter() > 1) NewtonRaphson=.TRUE.

  NonlinIter = GetInteger(GetSolverParams(), &
      'Nonlinear system max iterations',Found)
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

    Active = GetNOFBoundaryElements()
!$omp parallel do private(Element, n, nd, BC, Found)
    DO t=1,active
      Element => GetBoundaryElement(t)
      BC=>GetBC(Element)
      IF(.NOT.ASSOCIATED(BC)) CYCLE

      IF(GetLogical(BC,'Infinity BC',Found)) THEN
         n  = GetElementNOFNodes(Element)
         nd = GetElementNOFDOFs(Element)
         CALL LocalMatrixBC(  Element, n, nd )
      ELSE IF(GetLogical(BC,'Air Gap',Found)) THEN
         n  = GetElementNOFNodes( Element )
         nd = GetElementNOFDOFs( Element )
         CALL LocalMatrixAirGapBC(Element, BC, n, nd)
      END IF
    END DO
!$omp end parallel do

    CALL DefaultFinishAssembly()

    CALL SetMagneticFluxDensityBC()
    CALL DefaultDirichletBCs()
    Norm = DefaultSolve()
 
    IF( Solver % Variable % NonlinConverged == 1 ) EXIT
  END DO

   IF(.NOT. CSymmetry ) THEN
     CALL CalculateLumped(Model % NumberOfBodyForces)
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
   
CONTAINS

!------------------------------------------------------------------------------
 SUBROUTINE CalculateLumped(nbf)
!------------------------------------------------------------------------------
   INTEGER::nbf
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: a(nbf),IMoment,IA,TorqArea,rinner,router,ctorq, &
       xRe, xIm, torq
   COMPLEX(KIND=dp)::u(nbf)
   INTEGER :: i,bfid,n,nd
   TYPE(ValueList_t),POINTER::Params
!------------------------------------------------------------------------------

   U=0._dp; a=0._dp; torq=0._dp; TorqArea=0._dp; IMoment=0._dp;IA=0
   DO i=1,GetNOFActive()
     Element => GetActiveElement(i)
     nd = GetElementNOFDOFs(Element)
     n  = GetElementNOFNodes(Element)

     CALL Torque(Torq,TorqArea,Element,n,nd)

     Params=>GetBodyForce(Element)
     IF(ASSOCIATED(Params)) THEN
       bfid=GetBodyForceId(Element)
       IF(GetLogical(Params,'Calculate Potential',Found)) &
         CALL Potential(u(bfid),a(bfid),Element,n,nd)
     END IF

     Params=>GetBodyParams(Element)
     IF(ASSOCIATED(Params)) THEN
       IF(GetLogical(Params,'Calculate Inertial Moment',Found)) &
         CALL InertialMoment(IMoment,IA,Element,n,nd)
     END IF
   END DO

   DO i=1,nbf
     a(i) = ParallelReduction(a(i))
     xRe = REAL( u(i) ); xIm = AIMAG( u(i) )
     xRe = ParallelReduction(xRe)
     xIm = ParallelReduction(xIm)
     u(i) = CMPLX( xRe, xIm )
   END DO
   IMoment = ParallelReduction(IMoment)
   IA = ParallelReduction(IA)
   Torq = ParallelReduction(Torq)

   
   WRITE(Message,'(A,ES15.4)') 'Air gap initial torque:', Torq
   CALL Info('MagnetoDynamics2D',Message,Level=8)

   TorqArea = ParallelReduction(TorqArea)
   rinner = ListGetCRealAnyBody( Model,'r inner',Found )
   router = ListGetCRealAnyBody( Model,'r outer',Found )
   IF (TorqArea /= 0) THEN
      Ctorq = PI*(router**2-rinner**2) / TorqArea
   ELSE
      Ctorq = 0.0_dp
   END IF
   WRITE(Message,'(A,ES15.4)') 'Air gap correction:', cTorq
   CALL Info('MagnetoDynamics2D',Message,Level=8)
   Torq = Ctorq * Torq
   
   DO i=1,nbf
     IF(a(i)>0) THEN
       CALL ListAddConstReal(Model % Simulation,'res: Potential re / bodyforce ' &
                     //TRIM(i2s(i)),REAL(u(i))/a(i))
       CALL ListAddConstReal(Model % Simulation,'res: Potential im / bodyforce ' &
                     //TRIM(i2s(i)),AIMAG(u(i))/a(i))
     END IF
   END DO
   CALL ListAddConstReal(Model % Simulation,'res: air gap torque', Torq)
   CALL ListAddConstReal(Model % Simulation,'res: inertial volume', IA)
   CALL ListAddConstReal(Model % Simulation,'res: inertial moment', IMoment)

   WRITE(Message,'(A,ES15.4)') 'Air gap torque:', Torq
   CALL Info('MagnetoDynamics2D',Message,Level=7)
   WRITE(Message,'(A,ES15.4)') 'Inertial volume:', IA
   CALL Info('MagnetoDynamics2D',Message,Level=7)
   WRITE(Message,'(A,ES15.4)') 'Inertial moment:', Imoment
   CALL Info('MagnetoDynamics2D',Message,Level=7)
!------------------------------------------------------------------------------
 END SUBROUTINE CalculateLumped
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE InertialMoment(U,A,Element,n,nd)
!------------------------------------------------------------------------------
    INTEGER :: n,nd
    REAL(KIND=dp)::U,a
    TYPE(Element_t)::Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd), DetJ,x,y,r,Density(n)
    INTEGER :: t
    LOGICAL :: stat,Found
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
	!$OMP THREADPRIVATE(Nodes)

    Density(1:n) = GetReal(GetMaterial(),'Density',Found,UElement=Element)
    IF(.NOT.Found) RETURN

    CALL GetElementNodes( Nodes, Element )
  
    !Numerical integration:
    !----------------------
    IP = GaussPoints(Element)
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                IP % W(t), detJ, Basis )

      x = SUM(Nodes % x(1:nd)*Basis(1:nd))
      y = SUM(Nodes % y(1:nd)*Basis(1:nd))
      r = SQRT(x**2+y**2)
      A = A + IP % s(t)*detJ
      U = U + IP % s(t)*detJ*R*SUM(Density(1:n)*Basis(1:n))
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE InertialMoment
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE Torque(U,Area,Element,n,nd)
!------------------------------------------------------------------------------
    INTEGER :: n,nd
    REAL(KIND=dp) :: Area
    REAL(KIND=dp)::U
    TYPE(Element_t)::Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: dBasisdx(nd,3),Basis(nd), DetJ, &
             POT(2,nd),x,y,r,r0,r1
    COMPLEX(KIND=dp)::POTC(nd),Br,Bp,Bx,By
    REAL(KIND=dp)::BrRe,BpRe,BrIm,BpIm
    INTEGER :: t
    LOGICAL :: stat
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
	!$OMP THREADPRIVATE(Nodes)

    r0 = GetCReal(GetBodyParams(),'r inner',Found)
    r1 = GetCReal(GetBodyParams(),'r outer',Found)
    IF (.NOT.Found) RETURN

    CALL GetElementNodes( Nodes, Element )

    x = SUM(Nodes % x(1:n))/n
    y = SUM(Nodes % y(1:n))/n
    r = SQRT(x**2+y**2)
    IF (r<r0.OR.r>r1) RETURN

    CALL GetLocalSolution(POT, UElement=Element)
    POTC=POT(1,:)+im*POT(2,:)
  
    !Numerical integration:
    !----------------------
    IP = GaussPoints(Element)
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                IP % W(t), detJ, Basis, dBasisdx )

      x = SUM(Nodes % x(1:nd)*Basis(1:nd))
      y = SUM(Nodes % y(1:nd)*Basis(1:nd))
      r = SQRT(x**2+y**2)

      Bx =  SUM(POTC*dBasisdx(:,2))
      By = -SUM(POTC*dBasisdx(:,1))
      Br =  x/r*Bx + y/r*By
      Bp = -y/r*Bx + x/r*By
      BrRe = REAL( Br ); BrIm = AIMAG( Br )
      BpRe = REAL( Bp ); BpIm = AIMAG( Bp )


      U = U + IP % s(t)*detJ*r*(BrRe*BpRe+BrIm*BpIm)/(2*PI*4.0d-7*(r1-r0))
      Area = Area + IP % s(t)*detJ
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE Torque
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE Potential( U, A, Element,n,nd)
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp) :: U
    REAL(KIND=dp) :: A
    INTEGER :: n, nd
    TYPE(Element_t) :: Element

    REAL(KIND=dp) :: Basis(nd), DetJ,POT(2,nd),Omega
    COMPLEX(KIND=dp) ::  POTC(nd)
    INTEGER :: t
    LOGICAL :: stat
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP

    CALL GetElementNodes( Nodes, Element )

    CALL GetLocalSolution(POT, UElement=Element)
    POTC = POT(1,:) + im*POT(2,:)
    Omega = GetAngularFrequency(Found=Found)

    !Numerical integration:
    !----------------------
    IP = GaussPoints(Element)
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                IP % W(t), detJ, Basis )
      A = A + IP % s(t) * detJ
      U = U + IP % s(t) * detJ * im*Omega*SUM(POTC*Basis)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE Potential
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  Element, n, nd)
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,x
    INTEGER :: i,p,q,t,siz
    TYPE(GaussIntegrationPoints_t) :: IP
    COMPLEX(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LoadAtIp,&
      JAC(nd,nd),Agrad(3),Load(n),M(2,n),M_ip(2),POTC(nd), C(n), C_ip

    REAL(KIND=dp) :: POT(2,nd),Babs,mu,muder,Omega
    COMPLEX(KIND=dp) :: nu_tensor(2,2)
    COMPLEX(KIND=dp) :: R(2,2,n)       

    LOGICAL :: Cubic, HBcurve, Found, Stat, StrandedHomogenization

    REAL(KIND=dp), POINTER :: Bval(:), Hval(:), Cval(:), &
      CubicCoeff(:) => NULL(), HB(:,:) => NULL()
    TYPE(ValueListEntry_t), POINTER :: Lst
    TYPE(ValueList_t), POINTER :: Material,  BodyForce

    TYPE(Nodes_t), SAVE :: Nodes
    
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilType
    LOGICAL :: CoilBody    
    TYPE(ValueList_t), POINTER :: CompParams

    COMPLEX(KIND=dp) :: Bt(nd,2)
    COMPLEX(KIND=dp) :: Ht(nd,2) 
    COMPLEX(KIND=dp) :: B_ip(2), Alocal

    REAL(KIND=dp) :: nu_11(nd), nuim_11(nd), nu_22(nd), nuim_22(nd)
    REAL(KIND=dp) :: nu_val, nuim_val

    REAL(KIND=dp) :: foilthickness, coilthickness, nofturns, skindepth, mu0 
    COMPLEX(KIND=dp) :: FR
    LOGICAL :: InPlaneProximity = .FALSE.

    LOGICAL :: FoundIm

!$omp threadprivate(Nodes)
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes,Element )
    STIFF = 0._dp
    JAC  = 0._dp
    FORCE = 0._dp
    IF(TransientSimulation) MASS = 0._dp

    Material => GetMaterial(Element)

    Omega = GetAngularFrequency(Found=Found)
    InPlaneProximity = .FALSE.
   
    CoilBody = .FALSE.
    CompParams => GetComponentParams( Element )
    CoilType = ''
    StrandedHomogenization = .FALSE.
    IF (ASSOCIATED(CompParams)) THEN
      CoilType = GetString(CompParams, 'Coil Type', Found)
      IF (Found) THEN
        SELECT CASE (CoilType)
        CASE ('stranded')
           CoilBody = .TRUE.
           StrandedHomogenization = GetLogical(CompParams, 'Homogenization Model', Found)
           IF ( .NOT. Found ) StrandedHomogenization = .FALSE.
             
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
           CoilBody = .TRUE.
        CASE ('foil winding')
           CoilBody = .TRUE.
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
           CALL Fatal ('MagnetoDynamics2DHarmonic', 'Non existent Coil Type Chosen!')
        END SELECT
      END IF
    END IF

    CALL GetConstRealArray( Material, HB, 'H-B curve', HBCurve )
    siz = 0
    Cval => NULL()
    IF ( HBCurve ) THEN
      siz = SIZE(HB,1)
      IF(siz>1) THEN
        Bval=>HB(:,1)
        Hval=>HB(:,2)
        Cubic = GetLogical( Material, 'Cubic spline for H-B curve',Found)
        IF (Cubic.AND..NOT.ASSOCIATED(CubicCoeff)) THEN
          ALLOCATE(CubicCoeff(siz))
          CALL CubicSpline(siz,Bval,Hval,CubicCoeff)
        END IF
        Cval=>CubicCoeff
        HBCurve = .TRUE.
      END IF
    END IF

    IF(siz<=1) THEN
      Lst => ListFind(Material,'H-B Curve',HBcurve)
      IF(HBcurve) THEN
        Cval => Lst % CubicCoeff
        Bval => Lst % TValues
        Hval => Lst % FValues(1,1,:)
      END IF
    END IF

    Lst => ListFind(Material,'H-B Curve',HBcurve)
    IF(HBcurve) THEN
      CALL GetLocalSolution(POT,UElement=Element)
      POTC=POT(1,:)+im*POT(2,:)
      IF (.NOT. ASSOCIATED(Bval) ) CALL Fatal ('mgdyn2D','bval not associated')
      IF (.NOT. ASSOCIATED(Hval) ) CALL Fatal ('mgdyn2D','hval not associated')
    ELSE IF (.NOT. StrandedHomogenization) THEN 
      CALL GetReluctivity(Material,R,n,Element)
    END IF
 
    C = GetReal( Material, 'Electric Conductivity', Found, Element)
    C = C + im * GetReal( Material, 'Electric Conductivity im', Found, Element)

    M(1,:) = GetReal( Material, 'Magnetization 1', Found, Element)
    M(1,:) = M(1,:) + im*GetReal( Material, 'Magnetization 1 im', Found, Element)

    M(2,:) = GetReal( Material, 'Magnetization 2', Found, Element)
    M(2,:) = M(2,:) + im*GetReal( Material, 'Magnetization 2 im', Found, Element)

    Load = 0.0d0
    BodyForce => GetBodyForce(Element)
    IF ( ASSOCIATED(BodyForce) ) THEN
       Load(1:n) = GetReal( BodyForce, 'Current Density', Found, Element )
       Load(1:n) = Load(1:n) + im*GetReal( BodyForce, 'Current Density im', Found, Element )
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
        mu = InterpolateCurve(Bval,Hval,Babs,CubicCoeff=Cval)/Babs
        muder = (DerivateCurve(Bval,Hval,Babs,CubicCoeff=Cval)-mu)/Babs
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

      DO p=1,nd
        DO q=1,nd
          IF(CoilType /= 'stranded') STIFF(p,q) = STIFF(p,q) + &
                                   IP % s(t) * detJ * im * omega * C_ip * Basis(q)*Basis(p)
        END DO
      END DO

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

      IF (HBcurve.AND.NewtonRaphson) THEN
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

      FORCE(1:nd) = FORCE(1:nd) + IP % s(t) * DetJ * (LoadAtip * Basis(1:nd) + &
           (M_ip(1)*dBasisdx(1:nd,2)-M_ip(2)*dBasisdx(1:nd,1)))
    END DO

    IF (HBcurve.AND.NewtonRaphson) THEN
      STIFF = STIFF + JAC
      FORCE = FORCE + MATMUL(JAC,POTC)
    END IF

    IF(TransientSimulation) THEN
      CALL Default1stOrderTime( MASS, STIFF, FORCE, UElement=Element )
    END IF
    CALL DefaultUpdateEquations( STIFF, FORCE, UElement=Element )

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC(Element, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    LOGICAL :: Stat
    INTEGER :: i,p,q,t
    TYPE(GaussIntegrationPoints_t) :: IP

    REAL(KIND=dp) :: R_ip, &
            Inf_ip,Coord(3),Normal(3),mu,u,v
    
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
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixAirGapBC(Element, BC, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    LOGICAL :: Stat, Found
    INTEGER :: i,p,q,t
    TYPE(GaussIntegrationPoints_t) :: IP
    COMPLEX(KIND=dp) :: STIFF(nd,nd), FORCE(nd)
    REAL(KIND=dp) :: R(n), R_ip, &
            Coord(3),Normal(3),mu,u,v, AirGapLength(nd), &
            AirGapMu(nd), AirGapL

    TYPE(ValueList_t), POINTER :: Material
    TYPE(ValueList_t), POINTER :: BC

    TYPE(Element_t), POINTER :: Parent
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !$OMP THREADPRIVATE(Nodes)
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes, Element )
    STIFF = 0._dp
    FORCE = 0._dp

    AirGapLength=GetConstReal( BC, 'Air Gap Length', Found)
    if (.not. Found) CALL FATAL('LocalMatrixAirGapBC', 'Air Gap Length not found!')

    AirGapMu=GetConstReal( BC, 'Air Gap Relative Permeability', Found)
    if (.not. Found) AirGapMu=1d0

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
    REAL(KIND=dp) :: Bx(Solver % Mesh % MaxElementNodes), &
                      Bxim(Solver % Mesh % MaxElementNodes), &
                      By(Solver % Mesh % MaxElementNodes), &
                      Byim(Solver % Mesh % MaxElementNodes)
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
        IF ( ListCheckPresent( BC, 'Magnetic Flux Density 1') .OR. &
             ListCheckPresent( BC, 'Magnetic Flux Density 1 im') .OR. &
             ListCheckPresent( BC, 'Magnetic Flux Density 2') .OR. &
             ListCheckPresent( BC, 'Magnetic Flux Density 2 im') &
            ) THEN
          Bx = 0._dp
          Bxim = 0._dp
          By = 0._dp
          Byim = 0._dp

          Bx(1:n) = GetReal(BC, 'Magnetic Flux Density 1', Found)
          IF (.NOT. Found) Bx = 0._dp
          Bxim(1:n) = GetReal(BC, 'Magnetic Flux Density 1 im', Found)
          IF (.NOT. Found) Bxim = 0._dp
          By(1:n) = GetReal(BC, 'Magnetic Flux Density 2', Found)
          IF (.NOT. Found) By = 0._dp
          Byim(1:n) = GetReal(BC, 'Magnetic Flux Density 2 im', Found)
          IF (.NOT. Found) Byim = 0._dp
          DO j = 1,n
            k = Element % NodeIndexes(j)
            x = Mesh % Nodes % x(k)
            y = Mesh % Nodes % y(k)
            k = Perm(k)
            !b(2*k-1) = y * Bx(j) - x * By(j)
            !b(2*k) = y * Bxim(j) - x * Byim(j)

            CALL UpdateDirichletDof( A, 2*k-1, y * Bx(j) - x * By(j) )
            CALL UpdateDirichletDof( A, 2*k, y * Bxim(j) - x * Byim(j) )

            !CALL ZeroRow(A, 2*k-1)
            !CALL ZeroRow(A, 2*k)
            !CALL AddToCmplxMatrixElement(A, 2*k-1, 2*k-1, 1._dp, 0._dp)
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
  IF( .NOT. ListCheckPresent( SolverParams,'Variable') ) THEN
    CALL ListAddString( SolverParams, 'Variable','-nooutput bsolver_temp' )
  END IF
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


END SUBROUTINE Bsolver_init


!------------------------------------------------------------------------------
!> Given the vector potential computes its gradient i.e. the magnetic
!> field intensity.  
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
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName, CondName
  INTEGER :: i,j,k,dim,FluxDofs,firstmag,TotDofs
  LOGICAL :: ConstantBulkMatrix, ConstantBulkMatrixInUse
  LOGICAL :: GotIt, Visited = .FALSE.
  REAL(KIND=dp) :: Unorm, Totnorm, val
  REAL(KIND=dp), ALLOCATABLE, TARGET :: ForceVector(:,:)
  REAL(KIND=dp), POINTER :: SaveRHS(:)  
  TYPE(Variable_t), POINTER :: FluxSol, HeatingSol, JouleSol, AzSol
  LOGICAL ::  CSymmetry, LossEstimation, JouleHeating, ComplexPowerCompute,&
              AverageBCompute, BodyICompute, BodyVolumesCompute = .FALSE., &
              CirCompVolumesCompute = .FALSE., HomogenizationParamCompute, &
              LorentzForceCompute = .FALSE.
  TYPE(Matrix_t),POINTER::CM
  REAL(KIND=dp) :: Omega
  
  TYPE(Variable_t), POINTER :: CurrDensSol
  
  SAVE Visited

  CALL Info( 'BSolver', '-------------------------------------',Level=4 )
  CALL Info( 'BSolver', 'Computing the magnetic field density ',Level=4 )
  CALL Info( 'BSolver', '-------------------------------------',Level=4 )
  dim = CoordinateSystemDimension()

  CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )

!------------------------------------------------------------------------------
!  Check what needs to be computed
!------------------------------------------------------------------------------
  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
  IF ( COUNT( Solver % Variable % Perm > 0 ) <= 0 ) RETURN
  
  SolverParams => GetSolverParams()

  VarName = GetString(GetSolverParams(),'Target Variable',GotIt)
  IF(.NOT. GotIt) VarName = 'Potential'
  AzSol => VariableGet( Solver % Mesh % Variables, VarName ) 
  IF( .NOT. ASSOCIATED( AzSol ) ) THEN
    CALL Fatal('BSolver','Target field not present: '//TRIM(VarName) )
  END IF

  FluxSol => VariableGet(Solver % Mesh % Variables, 'B')
  IF( .NOT. ASSOCIATED( FluxSol ) ) THEN
    CALL Fatal('BSolver','Solution field not present: B' )
  END IF

  IF( FluxSol % Dofs / AzSol % Dofs /= 2 ) THEN
    WRITE( Message,'(A,I0,A)') 'B field should have ',2 * AzSol % Dofs,' dofs!'
    CALL Fatal('BSolver',Message)
  END IF


  FluxDofs = FluxSol % Dofs
  ConstantBulkMatrix = GetLogical( SolverParams, 'Constant Bulk Matrix', GotIt )
  ConstantBulkMatrixInUse = ConstantBulkMatrix .AND. &
      ASSOCIATED(Solver % Matrix % BulkValues)
  
  IF ( ConstantBulkMatrixInUse ) THEN
    Solver % Matrix % RHS = 0.0_dp
    Solver % Matrix % Values = Solver % Matrix % BulkValues        
  ELSE
    CALL DefaultInitialize()
  END IF
  
  TotDofs = FluxDofs
  JouleHeating = ListGetLogical( SolverParams, 'Calculate Joule Heating', GotIt )

  IF( JouleHeating ) THEN
    IF( FluxDofs /= 4 ) THEN
      CALL Fatal('BSolver','Joule heating can only be computed for complex problems!')
    ELSE
      TotDofs = TotDofs + 2
      HeatingSol => VariableGet(Solver % Mesh % Variables, 'Joule Heating')
      IF( .NOT. ASSOCIATED( HeatingSol ) ) THEN
        CALL Fatal('BSolver','Solution field not present: Joule Heating' )
      END IF
      JouleSol => VariableGet(Solver % Mesh % Variables, 'Joule Field')
      IF( .NOT. ASSOCIATED( JouleSol ) ) THEN
        CALL Fatal('BSolver','Solution field not present: Joule Field' )
      END IF
      TotDofs = TotDofs + 2
      CurrDensSol => VariableGet(Solver % Mesh % Variables, 'Current Density')
      IF( .NOT. ASSOCIATED( CurrDensSol ) ) THEN
        CALL Fatal('BSolver','Solution field not present: Current Density' )
      END IF
    END IF
  END IF

  !------------------------------------------------------------------------------
  ! In the case of time-harmonic analysis losses may be estimated in terms of B
  !------------------------------------------------------------------------------ 
  LossEstimation = GetLogical(SolverParams,'Loss Estimation',GotIt)
  IF( LossEstimation .AND. FluxDofs /= 4) THEN
    CALL Fatal( 'BSolver', 'Real solution, loss estimation omitted' )
  END IF

  HomogenizationParamCompute = GetLogical(SolverParams, 'Calculate Homogenization Parameters', GotIt)
  IF (.NOT. GotIt ) HomogenizationParamCompute = .FALSE.
  IF( HomogenizationParamCompute.AND. FluxDofs /= 4) THEN
    CALL Fatal( 'BSolver', 'Real solution, Calculate Homogenization Parameters omitted' )
  END IF

  ComplexPowerCompute = GetLogical(SolverParams,'Calculate Complex Power',GotIt)
  IF (.NOT. GotIt ) ComplexPowerCompute = .FALSE.
  IF( ComplexPowerCompute.AND. FluxDofs /= 4) THEN
    CALL Fatal( 'BSolver', 'Real solution, Complex Power omitted' )
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
  CALL DefaultFinishAssembly()
!        
!------------------------------------------------------------------------------     

   CALL DefaultDirichletBCs()
      
   TotNorm = 0.0_dp
   DO i=1,TotDofs
     Solver % Matrix % RHS => ForceVector(:,i)
     Solver % Variable % Values = 0
     UNorm = DefaultSolve()
     TotNorm = TotNorm + SUM(Solver % Variable % Values**2)
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
   CALL Info( 'BSolver', Message, Level=4 )


CONTAINS


!------------------------------------------------------------------------------
  SUBROUTINE BulkAssembly()
!------------------------------------------------------------------------------
       
    INTEGER :: elem,t,i,j,k,p,q,n,nd, Rank, BodyId
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: weight,coeff,detJ,BAtIp(8),PotAtIp(2),MuAtIp, &
        Omega,TotalHeating, DesiredHeating, HeatingCoeff
    COMPLEX(KIND=dp) :: CondAtIp
    REAL(KIND=dp) :: Freq, FreqPower, FieldPower, ComponentLoss(2), LossCoeff, &
        ValAtIp, ValAtIpim, TotalLoss, x
    LOGICAL :: Found, SetHeating
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
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilType, bodyNumber, XYNumber
    LOGICAL :: CoilBody, EddyLoss
    COMPLEX(KIND=dp) :: imag_value, imag_value2
    INTEGER :: IvarId, ReIndex, ImIndex, VvarDofs, VvarId
    REAL(KIND=DP) :: grads_coeff, nofturns
    REAL(KIND=DP) :: i_multiplier_re, i_multiplier_im, ModelDepth
    COMPLEX(KIND=dp) :: i_multiplier, Bx, By, Jz, LorentzForceDensX, &
                        LorentzForceDensY
    REAL(KIND=dp) :: ValueNorm

    INTEGER :: NofComponents=0, bid
    INTEGER, POINTER :: BodyIds(:)
    REAL(KIND=DP) :: Vol
    CHARACTER(LEN=MAX_NAME_LEN) :: CompNumber, OutputComp
    
    LOGICAL :: StrandedHomogenization, FoundIm

    REAL(KIND=dp), ALLOCATABLE :: sigma_33(:), sigmaim_33(:)

    LOGICAL :: LaminateModelPowerCompute=.FALSE., InPlaneProximity=.FALSE.
    REAL(KIND=dp) :: LaminatePowerDensity, BMagnAtIP, Fsk, Lambda, LaminateThickness, &
                     mu0=4d-7*PI, skindepth
    
    SAVE Nodes

    n = 2*MAX(Solver % Mesh % MaxElementDOFs,Solver % Mesh % MaxElementNodes)
    ALLOCATE( STIFF(n,n), FORCE(Totdofs,n) )
    ALLOCATE( POT(2,n), Basis(n), dBasisdx(n,3), alpha(n) )
    ALLOCATE( Cond(n), mu(n), sigma_33(n), sigmaim_33(n) ) 
    LagrangeVar => VariableGet( Solver % Mesh % Variables,'LagrangeMultiplier')
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
        BodyAverageBCompute(i) = ListGetLogical(Model % Bodies(i) % Values, 'Compute Average Magnetic Flux Density', Found)
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
      
      CoilType = ''
      CompParams => GetComponentParams( Element )
      StrandedHomogenization = .FALSE.
      InPlaneProximity = .FALSE.
      LaminateModelPowerCompute = .FALSE.
      IF (ASSOCIATED(CompParams)) THEN    
        CoilType = GetString(CompParams, 'Coil Type', Found)
        
        SELECT CASE (CoilType)
        CASE ('stranded')
          CoilBody = .TRUE.
 
          IvarId = GetInteger (CompParams, 'Circuit Current Variable Id', Found)
          IF (.NOT. Found) CALL Fatal ('MagnetoDynamicsCalcFields', 'Circuit Current Variable Id not found!')
 
          N_j = GetConstReal (CompParams, 'Stranded Coil N_j', Found)
          IF (.NOT. Found) CALL Fatal ('MagnetoDynamicsCalcFields', 'Stranded Coil N_j not found!')
 
          nofturns = GetConstReal(CompParams, 'Number of Turns', Found)
          IF (.NOT. Found) CALL Fatal('MagnetoDynamicsCalcFields','Stranded Coil: Number of Turns not found!')
          
          i_multiplier_re = GetConstReal(CompParams, 'Current Multiplier re', Found)
          IF (.NOT. Found) i_multiplier_re = 0._dp
          
          i_multiplier_im = GetConstReal(CompParams, 'Current Multiplier im', Found)
          IF (.NOT. Found) i_multiplier_im = 0._dp
          
          i_multiplier = i_multiplier_re + im * i_multiplier_im

          StrandedHomogenization = GetLogical(CompParams, 'Homogenization Model', Found)
          IF ( .NOT. Found ) StrandedHomogenization = .FALSE.

          IF ( StrandedHomogenization ) THEN 
!            nu_11 = 0._dp
!            nuim_11 = 0._dp
!            nu_11 = GetReal(CompParams, 'nu 11', Found)
!            nuim_11 = GetReal(CompParams, 'nu 11 im', FoundIm)
!            IF ( .NOT. Found .AND. .NOT. FoundIm ) CALL Fatal ('MagnetoDynamicsCalcFields', &
!                                                      'Homogenization Model nu 11 not found!')
!            nu_22 = 0._dp
!            nuim_22 = 0._dp
!            nu_22 = GetReal(CompParams, 'nu 22', Found)
!            nuim_22 = GetReal(CompParams, 'nu 22 im', FoundIm)
!            IF ( .NOT. Found .AND. .NOT. FoundIm ) CALL Fatal ('MagnetoDynamicsCalcFields', &
!                                                      'Homogenization Model nu 22 not found!')
            sigma_33 = GetReal(CompParams, 'sigma 33', Found)
            IF ( .NOT. Found ) sigma_33 = 0._dp
            sigmaim_33 = GetReal(CompParams, 'sigma 33 im', FoundIm)
            IF ( .NOT. FoundIm ) sigmaim_33 = 0._dp
            IF ( .NOT. Found .AND. .NOT. FoundIm ) CALL Fatal ('MagnetoDynamicsCalcFields', &
                                                                 'Homogenization Model Sigma 33 not found!')
          END IF
 
        CASE ('massive')
          CoilBody = .TRUE.

          VvarId = GetInteger (CompParams, 'Circuit Voltage Variable Id', Found)
          IF (.NOT. Found) CALL Fatal ('MagnetoDynamicsCalcFields', 'Circuit Voltage Variable Id not found!')

        CASE ('foil winding')
          CoilBody = .TRUE.
          CALL GetLocalSolution(alpha,'Alpha')

          VvarId = GetInteger (CompParams, 'Circuit Voltage Variable Id', Found)
          IF (.NOT. Found) CALL Fatal ('MagnetoDynamicsCalcFields', 'Circuit Voltage Variable Id not found!')

          coilthickness = GetConstReal(CompParams, 'Coil Thickness', Found)
          IF (.NOT. Found) CALL Fatal('MagnetoDynamicsCalcFields','Foil Winding: Coil Thickness not found!')
 
          nofturns = GetConstReal(CompParams, 'Number of Turns', Found)
          IF (.NOT. Found) CALL Fatal('MagnetoDynamicsCalcFields','Foil Winding: Number of Turns not found!')
 
          VvarDofs = GetInteger (CompParams, 'Circuit Voltage Variable dofs', Found)
          IF (.NOT. Found) CALL Fatal ('MagnetoDynamicsCalcFields', 'Circuit Voltage Variable dofs not found!')
          InPlaneProximity = GetLogical(CompParams, 'Foil In Plane Proximity', Found)
          IF (InPlaneProximity) THEN
             LaminateThickness = coilthickness/nofturns
             LaminateModelPowerCompute = .TRUE.
          END IF
 
        CASE DEFAULT
          CALL Fatal ('BSolver', 'Non existent Coil Type Chosen!')
        END SELECT
      END IF

      
      ! Integrate local stresses:
      ! -------------------------
      IntegStuff = GaussPoints( Element )
      STIFF  = 0.0_dp
      FORCE  = 0.0_dp

      CALL GetLocalSolution( POT, VarName )

      IF( JouleHeating ) THEN
        BodyId = GetBody() 
        Material => GetMaterial()
        Cond(1:n) = GetReal( Material, 'Electric Conductivity', Found, Element)
      END IF

      IF( LossEstimation ) THEN
        BodyId = GetBody() 
        LossCoeff = ListGetFun( Material,'Fourier Loss Coefficient',Freq,Found )
        EddyLoss = .FALSE.
        IF (.NOT. Found) EddyLoss = .TRUE.
      END IF
      
      IF (BodyVolumesCompute) THEN
        BodyId = GetBody()
        BodyVolumes(BodyId) = 0._dp
      END IF

      IF (ComplexPowerCompute) THEN
        BodyId = GetBody()
        Material => GetMaterial()

        IF (StrandedHomogenization) CALL Fatal ('MagnetoDynamics2D','Calculate Complex Power for Stranded & 
                                                 Homogenization model is not implemented.')

        mu = GetReal( Material, 'Relative Permeability', Found)
        mu = mu * 4.d-7*PI
        IF ( .NOT. Found ) CALL Warn('BSolver', 'Relative Permeability not found!')
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
                                                         
          IF (CoilType /= 'stranded') THEN
            PotAtIp(1) =   Omega * SUM(POT(2,1:nd) * Basis(1:nd))
            PotAtIp(2) = - Omega * SUM(POT(1,1:nd) * Basis(1:nd))
          ELSE
            PotAtIp(1) = 0._dp
            PotAtIp(2) = 0._dp
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
          BMagnAtIP = SQRT(ABS(imag_value**2._dp) + ABS(imag_value2**2._dp))
        END IF
        
        IF (LorentzForceCompute) THEN
          BodyId = GetBody()
          ! Let's compute the JxB for all the bodies and 
          ! then we sum from these for the components which are outputed.

          Bx = CMPLX(BatIp(1), BatIp(3), KIND=dp)
          By = CMPLX(BatIp(2), BatIp(4), KIND=dp)
          Jz = CMPLX(BatIp(7), BatIp(8), KIND=dp)

          LorentzForceDensX = ModelDepth * Weight * By / Jz * ABS(Jz)**2._dp
          LorentzForceDensY = -ModelDepth * Weight * Bx / Jz * ABS(Jz)**2._dp
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
          cmplx_power = cmplx_power + im * ModelDepth * Weight * Omega/MuAtIp * (ABS(imag_value)**2._dp+ABS(imag_value2)**2._dp)

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
        Solver % Matrix % Rhs => SaveRhs
        CALL DefaultUpdateEquations( STIFF, FORCE(1,1:nd), &
            BulkUpdate=ConstantBulkMatrix )
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
      CALL Info('MagnetoDynamics2D',Message)
      CALL ListAddConstReal( Model % Simulation, 'res: Joule heating',TotalHeating)
      
      DesiredHeating = ListGetConstReal( SolverParams, &
          'Desired Heating Power',Found)        
      IF( Found .AND. TotalHeating > 0.0_dp ) THEN
        HeatingCoeff = DesiredHeating / TotalHeating

        WRITE(Message,'(A,ES15.4)') 'Joule coefficient: ',HeatingCoeff
        CALL Info('MagnetoDynamics2D',Message)
        CALL ListAddConstReal( Model % Simulation, 'res: Joule coefficient',HeatingCoeff)
      
        ForceVector(:,5) = HeatingCoeff * ForceVector(:,5) 
        ForceVector(:,6) = HeatingCoeff * ForceVector(:,6) 
      END IF
    END IF

    ! Assembly of the face terms:
    !----------------------------
    IF (GetLogical(GetSolverParams(),'Discontinuous Galerkin',Found)) THEN
      IF (GetLogical(GetSolverParams(),'Average Within Materials',Found)) THEN
        FORCE = 0.0d0
        CALL AddLocalFaceTerms( STIFF, FORCE(1,:) )
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
      !CALL ListAddConstReal( Model % Simulation,'res: cos mode fourier loss', ComponentLoss(1)) 
      !CALL ListAddConstReal( Model % Simulation,'res: sin mode fourier loss', ComponentLoss(2))       
    
      !---------------------------------------------------------------------------------
      ! Screen ouput for componentwise and bodywise losses 
      !--------------------------------------------------------------------------------
      WRITE( Message,'(A,ES12.3)') 'Loss for cos mode: ', ComponentLoss(1)
      CALL Info('BSolver', Message, Level=6 )
      WRITE( Message,'(A,ES12.3)') 'Loss for sin mode: ', ComponentLoss(2)
      CALL Info('BSolver', Message, Level=6 )
      WRITE( Message,'(A,ES12.3)') 'Total loss: ',TotalLoss
      CALL Info('BSolver',Message, Level=5 )

      CALL Info('FourierLosses','Losses by bodies',Level=6)
      DO j=1,Model % NumberOfBodies
         IF( BodyLoss(j) < TINY( TotalLoss ) ) CYCLE
         WRITE( Message,'(A,I0,A,ES12.3)') 'Body ',j,' : ',BodyLoss(j)
         WRITE (bodyNumber, "(I0)") j
         CALL ListAddConstReal( Model % Simulation,'res: Loss in Body '//TRIM(bodyNumber)//':', BodyLoss(j) )
         CALL Info('FourierLosses', Message, Level=6 )
      END DO

      DEALLOCATE( BodyLoss )
    END IF

    IF (LorentzForceCompute) THEN
       DO j=1,Model % NumberOfBodies
         DO i = 1, 2
           BodyLorentzForcesRe(i,j) = ParallelReduction(BodyLorentzForcesRe(i,j))
           BodyLorentzForcesIm(i,j) = ParallelReduction(BodyLorentzForcesIm(i,j))
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
                 in Component '//TRIM(i2s(j)), ComponentLorenzForcesRe(1,j) )
                         
           CALL ListAddConstReal( Model % Simulation,'res: Lorentz Force 2 re & 
                 in Component '//TRIM(i2s(j)), ComponentLorenzForcesRe(2,j) )

           CALL ListAddConstReal( Model % Simulation,'res: Lorentz Force 1 im & 
                 in Component '//TRIM(i2s(j)), ComponentLorenzForcesIm(1,j) )
                         
           CALL ListAddConstReal( Model % Simulation,'res: Lorentz Force 2 im & 
                 in Component '//TRIM(i2s(j)), ComponentLorenzForcesIm(2,j) )

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
                 in Component '//TRIM(i2s(j)), CirCompComplexPower(1,j) )
                         
           CALL ListAddConstReal( Model % Simulation,'res: Power im & 
                 in Component '//TRIM(i2s(j)), CirCompComplexPower(2,j) )
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
      


    DEALLOCATE( POT, STIFF, FORCE, Basis, dBasisdx, mu, Cond, sigma_33, sigmaim_33 )

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
      ProxNu(2) = AIMAG(imag_value) 
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
      INTEGER :: i, j, p, q, t, nFace, nParent
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

