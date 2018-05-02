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
SUBROUTINE MagnetoDynamics2D_Init( Model,Solver,dt,TransientSimulation ) ! {{{
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
END SUBROUTINE MagnetoDynamics2D_Init ! }}}
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solver the magnetic vector potential in cartesian 2D case.
!> The solver may take into account rotating boundary conditions.
!> Also optionally compute moments and inertia. 
!------------------------------------------------------------------------------
SUBROUTINE MagnetoDynamics2D( Model,Solver,dt,TransientSimulation ) ! {{{
!------------------------------------------------------------------------------
  USE DefUtils
  USE CircuitUtils
  USE ZirkaUtils
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
  INTEGER :: i,j,k,n, nb, nd, t, istat, Active, NonlinIter, iter, MinNonlinIter

  TYPE(ValueList_t), POINTER :: BC
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:), &
               NEWX(:), NEWY(:), POT(:)

  TYPE(Mesh_t),   POINTER :: Mesh

  LOGICAL :: NewtonRaphson = .FALSE., CSymmetry
  INTEGER :: CoupledIter
  TYPE(Variable_t), POINTER :: IterV, CoordVar

  TYPE(Matrix_t),POINTER::CM
  TYPE(GlobalHysteresisModel_t), POINTER :: ZirkaModel

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
  MinNonlinIter = GetInteger(GetSolverParams(), &
           'Nonlinear System Min Iterations',Found)
  IF(.NOT.Found) MinNonlinIter = 1


  CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )

  call info('MagnetoDynamics2D','Initializing Zirka hysteresis models', Level=10)
  call InitHysteresis(Model, Solver)

  DO iter = 1,NonlinIter
    !IF(Iter > 1) NewtonRaphson=.TRUE.
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
 
    IF( Solver % Variable % NonlinConverged > 0 .and. iter >= MinNonLinIter) EXIT
  END DO

  ! For cylindrical symmetry the model lumping has not been implemented
  IF( .NOT. CSymmetry ) THEN
    CALL CalculateLumped(Model % NumberOfBodyForces)
  END IF

  CoordVar => VariableGet(Mesh % Variables,'Coordinates')
  CALL DriveHysteresis(model, solver)

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
 SUBROUTINE CalculateLumped(nbf) ! {{{
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
 END SUBROUTINE CalculateLumped ! }}}
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE InertialMoment(U,A,Element,n,nd) ! {{{
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
  END SUBROUTINE InertialMoment ! }}}
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE Torque(U,Area,Element,n,nd) ! {{{
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
  END SUBROUTINE Torque !}}}
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE Potential( U, A, Element,n,nd) ! {{{
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
  END SUBROUTINE Potential ! }}}
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(Element, n, nd) ! {{{
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

    REAL(KIND=dp), POINTER :: Bval(:), Hval(:), Cval(:), &
      CubicCoeff(:), HB(:,:)
    TYPE(ValueListEntry_t), POINTER :: Lst
    TYPE(ValueList_t), POINTER :: Material, BodyForce

    TYPE(Nodes_t), SAVE :: Nodes
!$omp threadprivate(Nodes)
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilType
    LOGICAL :: CoilBody    
    TYPE(ValueList_t), POINTER :: CompParams

    REAL(KIND=dp) :: Bt(nd,2), Ht(nd,2)
    REAL(KIND=dp) :: nu_tensor(2,2)
    REAL(KIND=dp) :: B_ip(2), Alocal, H_ip(2)

    ! Zirka related
    LOGICAL :: Zirka
    TYPE(Variable_t), POINTER :: hystvar
    TYPE(GlobalHysteresisModel_t), pointer :: zirkamodel

!------------------------------------------------------------------------------
    CubicCoeff => NULL()
    HB => NULL()
    CALL GetElementNodes( Nodes,Element )
    STIFF = 0._dp
    JAC  = 0._dp
    FORCE = 0._dp
    IF(TransientSimulation) MASS = 0._dp

    Material => GetMaterial(Element)

    CALL GetConstRealArray( Material, HB, 'H-B curve', HBCurve )
    Zirka = ListGetLogical(Material, 'Zirka material', Zirka)

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

    if (zirka) then
      CALL GetLocalSolution(POT,UElement=Element,USolver=Solver)
      zirkamodel => GetZirkaPointer(Material)
      hystvar => GetZirkaVariable(Material)
    end if

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

      if(Zirka .or. HBCUrve) then
        Agrad = 0.0_dp
        Agrad = MATMUL( POT,dBasisdx )
        Alocal = SUM( POT(1:n) * Basis(1:n) )
        ! Sign? This convention: \vec A = A u_z
        ! -----
        B_ip(1) = Agrad(2) 
        B_ip(2) = -Agrad(1)
        IF( CSymmetry ) then
          B_ip = -B_ip
          B_ip(2) = B_ip(2) + Alocal/x
        end if
      end if

      IF (HBcurve ) THEN
        ! -----
        Babs = MAX( SQRT(SUM(B_ip**2)), 1.d-8 )
        mu = InterpolateCurve(Bval,Hval,Babs,CubicCoeff=Cval)/Babs
        muder = (DerivateCurve(Bval,Hval,Babs,CubicCoeff=Cval)-mu)/Babs
        nu_tensor(1,1) = mu ! Mu is really nu!!! too lazy to correct now...
        nu_tensor(2,2) = mu
      ELSEIF(Zirka) THEN
        call GetZirkaHBAtIP(t, solver, element, hystvar, zirkamodel, B_ip, H_ip, nu_tensor)
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
      Bt(:,1) =  dbasisdx(:,2)
      Bt(:,2) = -dbasisdx(:,1)
      IF ( CSymmetry ) then
        Bt(:,1:2) = -Bt(:,1:2)
        Bt(:,2) = Bt(:,2) + Basis(:)/x
      end if

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
      IF(zirka) then
        FORCE(1:nd) = FORCE(1:nd) - (H_ip(1)*Bt(1:nd,1) + H_ip(2)*Bt(1:nd,2)) * IP % s(t) * detJ
      END IF
    END DO

    IF (HBcurve .AND. NewtonRaphson) THEN
      STIFF = STIFF + JAC
      FORCE = FORCE + MATMUL(JAC,POT)
    END IF

    IF(Zirka) THEN
      FORCE = FORCE + MATMUL(STIFF, POT)
    END IF

    IF(TransientSimulation) THEN
      CALL Default1stOrderTime( MASS, STIFF, FORCE,UElement=Element, USolver=Solver )
    END IF
    CALL DefaultUpdateEquations( STIFF, FORCE,UElement=Element, USolver=Solver)

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix ! }}}
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
  REAL(KIND=dp), intent(INOUT) :: dHdB(2,2)
!-------------------------------------------------------------------------------
  INTEGER :: ipindex, n_dir, k,l
  REAL(KIND=dp) :: dH, B0(3)
!-------------------------------------------------------------------------------
  ipindex = getipindex(i_IP, usolver=solver, element=element, ipvar=hystvar)
  IF (ipindex /= 0 ) THEN
  H_ip = 0.0_dp
    Do n_dir = 1, ubound(zirkamodel % curves, 1)
      B0 = zirkamodel % curves(n_dir, ipindex) % B0
      associate(Bdir => sum(B_ip*B0(1:2)))
        ! H_ip(1:2) = H_ip(1:2) + zirkamodel % curves(n_dir,ipindex) % &
        !     eval(sum(B_ip*B0(1:2)), cached = .true., dhdb=dH) * &
        !     B0(1:2)
        H_ip(1:2) = H_ip(1:2) + zirkamodel % curves(n_dir,ipindex) % &
            eval(Bdir, cached = .true., dhdb=dH) * &
            B0(1:2)
      end associate
      DO k = 1,2
        DO l = 1,2
          dHdB(k,l) = dHdB(k,l) + dH*B0(k)*B0(l)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE ! }}}
!-------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC(Element, n, nd ) ! {{{
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    LOGICAL :: Stat, found_mu_2
    INTEGER :: i,p,q,t
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), R(2,2,n), R_ip, &
            Inf_ip,Coord(3),Normal(3),mu,u,v, mu_2

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
  END SUBROUTINE LocalMatrixBC ! }}}
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixAirGapBC(Element, BC, n, nd ) ! {{{
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
  END SUBROUTINE LocalMatrixAirGapBC ! }}}
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE SetMagneticFluxDensityBC() ! {{{
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
  END SUBROUTINE SetMagneticFluxDensityBC ! }}}
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE GetReluctivity(Material,Acoef,n,Element) ! {{{
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
  END SUBROUTINE GetReluctivity ! }}}
!------------------------------------------------------------------------------

END SUBROUTINE MagnetoDynamics2D ! }}}
!------------------------------------------------------------------------------




!> \}

