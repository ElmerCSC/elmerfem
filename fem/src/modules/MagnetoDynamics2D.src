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
! *  Authors: Juha Ruokolainen, Mika Malinen, Peter Råback
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

  Params => Solver % Values
  CALL ListAddInteger( Params, 'Variable Dofs',1 )
  IF( .NOT. ListCheckPresent(  Params,'Variable') ) THEN
    CALL ListAddString( Params,'Variable','Potential')
  END IF

  CALL ListAddLogical( Solver % Values, 'Apply Mortar BCs',.TRUE.)
!------------------------------------------------------------------------------
END SUBROUTINE MagnetoDynamics2D_Init
!------------------------------------------------------------------------------
!> \ingroup Solvers
!> \{
!> \}

!------------------------------------------------------------------------------
!> Solver the magnetic vector potential in cartesian 2D case.
!> The solver may take into account rotating boundary conditions.
!> Also optionally compute moments and inertia. 
!------------------------------------------------------------------------------
SUBROUTINE MagnetoDynamics2D( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
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

  TYPE(Matrix_t),POINTER::A,B,CM,S

!------------------------------------------------------------------------------

  CALL Info( 'MagnetoDynamics2D','------------------------------------------------', Level=4 )
  CALL Info( 'MagnetoDynamics2D', 'Solving equation for magnetic vector potential', Level=4 )
  CALL Info( 'MagnetoDynamics2D','------------------------------------------------', Level=4 )


  ! Allocate some permanent storage, this is done first time only:
  ! --------------------------------------------------------------
  NULLIFY(BC)
  Mesh => GetMesh()

  IF(GetCoupledIter()>1) NewtonRaphson=.TRUE.

  NonlinIter = GetInteger(GetSolverParams(), &
           'Nonlinear System Max Iterations',Found)
  IF(.NOT.Found) NonlinIter=1

  CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )

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
      END IF
    END DO
!$omp end parallel do

    CALL DefaultFinishAssembly()

    CALL DefaultDirichletBCs()
    Norm = DefaultSolve()
 
    IF( Solver % Variable % NonlinConverged == 1 ) EXIT
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

CONTAINS

!------------------------------------------------------------------------------
 SUBROUTINE CalculateLumped(nbf)
!------------------------------------------------------------------------------
   INTEGER::nbf
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: torq,a(nbf),u(nbf),IMoment,IA
   INTEGER :: i,bfid,n,nd
   LOGICAL :: Found
   TYPE(ValueList_t),POINTER::Params
!------------------------------------------------------------------------------

   U=0._dp; a=0._dp; torq=0._dp; IMoment=0._dp;IA=0
   DO i=1,GetNOFActive()
     Element => GetActiveElement(i)
     nd = GetElementNOFDOFs(Element)
     n  = GetElementNOFNodes(Element)

     CALL Torque(Torq,Element,n,nd)

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

   DO i=1,nbf
     IF(a(i)>0) THEN
       CALL ListAddConstReal(Model % Simulation,'res: Potential / bodyforce ' &
                     //TRIM(i2s(i)),u(i)/a(i))
       CALL ListAddConstReal(Model % Simulation,'res: area / bodyforce ' &
                     //TRIM(i2s(i)),a(i))
     END IF
   END DO
   CALL ListAddConstReal(Model % Simulation,'res: Air Gap Torque', Torq)
   CALL ListAddConstReal(Model % Simulation,'res: Inertial Volume', IA)
   CALL ListAddConstReal(Model % Simulation,'res: Inertial Moment', IMoment)

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
  SUBROUTINE Torque(U,Element,n,nd)
!------------------------------------------------------------------------------
    INTEGER :: n,nd
    REAL(KIND=dp)::U
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
    TYPE(Element_t) :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    INTEGER :: i,p,q,t,siz

    TYPE(GaussIntegrationPoints_t) :: IP

    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), &
      LOAD(nd),R(n),C(n), mu,muder,Babs,POT(nd), &
        JAC(nd,nd),Agrad(3),C_ip,M(2,n),M_ip(2),x

    LOGICAL :: Cubic, HBcurve, Found, Stat

    REAL(KIND=dp), POINTER :: Bval(:), Hval(:), Cval(:)
    TYPE(ValueList_t), POINTER :: Material, Lst, BodyForce

    TYPE(Nodes_t), SAVE :: Nodes
!$omp threadprivate(Nodes)
!------------------------------------------------------------------------------

    CALL GetElementNodes( Nodes,Element )
    STIFF = 0._dp
    JAC  = 0._dp
    FORCE = 0._dp
    IF(TransientSimulation) MASS = 0._dp

    Material => GetMaterial(Element)

    Lst => ListFind(Material,'H-B Curve',HBcurve)
    IF(HBcurve) THEN
      CALL GetLocalSolution(POT,UElement=Element,USolver=Solver)
      Cval => Lst % CubicCoeff
      Bval => Lst % TValues
      Hval => Lst % FValues(1,1,:)
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

      IF (HBcurve) THEN
        Agrad = MATMUL( POT,dBasisdx )
        Babs = MAX( SQRT(SUM(Agrad**2)), 1.d-8 )
        mu = InterpolateCurve(Bval,Hval,Babs,CubicCoeff=Cval)/Babs
        muder = (DerivateCurve(Bval,Hval,Babs,CubicCoeff=Cval)-mu)/Babs
      ELSE
        muder=0._dp
        mu = SUM( Basis(1:n) * R(1:n) )
      END IF

      C_ip = SUM( Basis(1:n) * C(1:n) )
      M_ip = MATMUL( M,Basis(1:n) )


      ! Finally, the elemental matrix & vector:
      !----------------------------------------
      IF(TransientSimulation .AND. C_ip/=0._dp) THEN
        DO p=1,nd
          DO q=1,nd
            MASS(p,q) = MASS(p,q) + IP % s(t) * detJ * C_ip * Basis(q)*Basis(p)
          END DO
        END DO
      END IF

      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + IP % s(t) * DetJ * &
             mu*MATMUL(dBasisdx, TRANSPOSE(dBasisdx))

      IF( Csymmetry ) THEN
        DO p = 1,nd
          DO q = 1,nd
            STIFF(p,q) = STIFF(p,q) + IP % s(t) * DetJ * &
                (mu/x) * ( Basis(p)*dBasisdx(q,1) + &
                Basis(q)*dBasisdx(p,1) + Basis(p)*Basis(q)/x )              
          END DO
        END DO
      END IF

      ! Csymmetry is not yet considered in the Newton linearization
      IF (HBcurve .AND. NewtonRaphson) THEN
        DO p=1,nd
          DO q=1,nd
            JAC(p,q) = JAC(p,q) + IP % s(t) * DetJ * &
              muder/babs*SUM(Agrad*dBasisdx(q,:))*SUM(Agrad*dBasisdx(p,:))
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
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), R(n), R_ip, &
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
    Material => GetMaterial(Parent)
    CALL GetReluctivity(Material,R,n,Element)

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                 IP % W(t), detJ, Basis )


      mu = SUM(Basis(1:n)*R(1:n))

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
 SUBROUTINE GetReluctivity(Material,Acoef,n,Element)
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Material
    INTEGER :: n
    REAL(KIND=dp) :: Acoef(:)
    TYPE(Element_t), OPTIONAL :: Element
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

    Acoef(1:n) = GetReal( Material, 'Relative Permeability', Found,Element )
    IF ( Found ) THEN
      Acoef(1:n) = Avacuum * Acoef(1:n)
    ELSE
      Acoef(1:n) = GetReal( Material, 'Permeability', Found,Element )
    END IF
    IF ( Found ) THEN
      Acoef(1:n) = 1._dp / Acoef(1:n)
    ELSE
      Acoef(1:n) = GetReal( Material, 'Reluctivity', Found,Element )
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

  Params => Solver % Values
  CALL ListAddInteger( Params, 'Variable Dofs',2 )
  IF( .NOT. ListCheckPresent(  Params,'Variable') ) THEN
    CALL ListAddString( Params,'Variable',&
        'Potential[Potential re:1 Potential im:1]')
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
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: TransientSimulation !< Steady state or transient simulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  LOGICAL :: AllocationsDone = .FALSE., Found,FoundMortar
  TYPE(Element_t),POINTER :: Element

  REAL(KIND=dp) :: Norm
  INTEGER :: i,j,k,ip,jp,n, nb, nd, t, istat, Active, iter, NonlinIter

  TYPE(ValueList_t), POINTER :: BC
  TYPE(Mesh_t),   POINTER :: Mesh
  COMPLEX(KIND=dp), PARAMETER :: im=(0._dp,1._dp)


  LOGICAL, SAVE :: NewtonRaphson = .FALSE., CSymmetry
  INTEGER :: CoupledIter
  TYPE(Variable_t), POINTER :: IterV, CoordVar

  TYPE(Matrix_t),POINTER::A,B,CM,S

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

  S=>Solver % Matrix
  S % COMPLEX = .TRUE.

  DO i=1,Model % NumberOFBCs
    j=GetInteger(Model % BCs(i) % Values,'Mortar BC',FoundMortar)
    IF(FoundMortar) THEN
      S % COMPLEX = .FALSE.

      A => PeriodicProjector(Model,Mesh,i,j,2,.TRUE.)
      FoundMortar = ASSOCIATED(A)
      EXIT
    END IF
  END DO

  IF(FoundMortar) THEN
    CM => AllocateMatrix()
    CM % FORMAT = MATRIX_LIST
    DO i=1,A % NumberOfRows
      DO j=A % Rows(i),A % Rows(i+1)-1
        jp = Solver % Variable % Perm(A % Cols(j))
        CALL AddToMatrixElement(CM,2*(i-1)+1,2*(jp-1)+1,A % Values(j))
        CALL AddToMatrixElement(CM,2*(i-1)+1,2*(jp-1)+2,A % Values(j))
      END DO
    END DO
    CALL FreeMatrix(A)

    CALL List_toCRSMatrix(CM)
    ALLOCATE(CM % RHS(CM % NumberOfRows))
    CM % RHS = 0._dp
    S % ConstraintMatrix => CM
  END IF

  IF(GetCoupledIter()>1) NewtonRaphson=.TRUE.

  NonlinIter = GetInteger(GetSolverParams(), &
           'Nonlinear system max iterations',Found)
  IF(.NOT.Found) NonlinIter=1

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
      END IF
    END DO
!$omp end parallel do

    CALL DefaultFinishAssembly()

    CALL DefaultDirichletBCs()
    Norm = DefaultSolve()
 
    IF( Solver % Variable % NonlinConverged == 1 ) EXIT
  END DO

   IF(FoundMortar) THEN
     CALL FreeMatrix(CM)
     S % ConstraintMatrix => NULL()
   END IF

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

CONTAINS

!------------------------------------------------------------------------------
 SUBROUTINE CalculateLumped(nbf)
!------------------------------------------------------------------------------
   INTEGER::nbf
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: a(nbf),IMoment,IA
   COMPLEX(KIND=dp)::u(nbf),torq
   INTEGER :: i,bfid,n,nd
   TYPE(ValueList_t),POINTER::Params
!------------------------------------------------------------------------------

   U=0._dp; a=0._dp; torq=0._dp; IMoment=0._dp;IA=0
   DO i=1,GetNOFActive()
     Element => GetActiveElement(i)
     nd = GetElementNOFDOFs(Element)
     n  = GetElementNOFNodes(Element)

     CALL Torque(Torq,Element,n,nd)

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
     IF(a(i)>0) THEN
       CALL ListAddConstReal(Model % Simulation,'res: Potential re / bodyforce ' &
                     //TRIM(i2s(i)),REAL(u(i))/a(i))
       CALL ListAddConstReal(Model % Simulation,'res: Potential im / bodyforce ' &
                     //TRIM(i2s(i)),AIMAG(u(i))/a(i))
     END IF
   END DO
   CALL ListAddConstReal(Model % Simulation,'res: Air Gap Torque re', REAL(Torq))
   CALL ListAddConstReal(Model % Simulation,'res: Air Gap Torque im', AIMAG(Torq))
   CALL ListAddConstReal(Model % Simulation,'res: Inertial Volume', IA)
   CALL ListAddConstReal(Model % Simulation,'res: Inertial Moment', IMoment)
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
  SUBROUTINE Torque(U,Element,n,nd)
!------------------------------------------------------------------------------
    INTEGER :: n,nd
    COMPLEX(KIND=dp)::U
    TYPE(Element_t)::Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: dBasisdx(nd,3),Basis(nd), DetJ, &
             POT(2,nd),x,y,r,r0,r1
    COMPLEX(KIND=dp)::POTC(nd),Br,Bp,Bx,By
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
      U = U + IP % s(t)*detJ*r*Br*Bp/(PI*4.0d-7*(r1-r0))
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

    REAL(KIND=dp) :: POT(2,nd),R(n),Babs,mu,muder,Omega

    LOGICAL :: Cubic, HBcurve, Found, Stat

    REAL(KIND=dp), POINTER :: Bval(:), Hval(:), Cval(:)
    TYPE(ValueList_t), POINTER :: Material, Lst, BodyForce

    TYPE(Nodes_t), SAVE :: Nodes
!$omp threadprivate(Nodes)
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes,Element )
    STIFF = 0._dp
    JAC  = 0._dp
    FORCE = 0._dp
    IF(TransientSimulation) MASS = 0._dp

    Material => GetMaterial(Element)

    Omega = GetAngularFrequency(Found=Found)

    Lst => ListFind(Material,'H-B Curve',HBcurve)
    IF(HBcurve) THEN
      CALL GetLocalSolution(POT,UElement=Element)
      POTC=POT(1,:)+im*POT(2,:)
      Cval => Lst % CubicCoeff
      Bval => Lst % TValues
      Hval => Lst % FValues(1,1,:)
    ELSE
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

      IF (HBcurve) THEN
        Agrad = MATMUL( POTC,dBasisdx )
        Babs = MAX( SQRT(SUM(ABS(Agrad)**2)), 1.d-8 )
        mu = InterpolateCurve(Bval,Hval,Babs,CubicCoeff=Cval)/Babs
        muder = (DerivateCurve(Bval,Hval,Babs,CubicCoeff=Cval)-mu)/Babs
      ELSE
        muder=0._dp
        mu = SUM( Basis(1:n) * R(1:n) )
      END IF

      C_ip = SUM( Basis(1:n) * C(1:n) )
      M_ip = MATMUL( M,Basis(1:n) )

      DO p=1,nd
        DO q=1,nd
          STIFF(p,q) = STIFF(p,q) + IP % s(t) * detJ * im * omega * C_ip * Basis(q)*Basis(p)
        END DO
      END DO

      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + IP % s(t) * DetJ * &
             mu*MATMUL(dBasisdx, TRANSPOSE(dBasisdx))

      IF( Csymmetry ) THEN
        DO p = 1,nd
          DO q = 1,nd
            STIFF(p,q) = STIFF(p,q) + IP % s(t) * DetJ * &
                 (mu/x) * ( Basis(p)*dBasisdx(q,1) + &
                 Basis(q)*dBasisdx(p,1) + Basis(p)*Basis(q)/x )              
          END DO
        END DO
      END IF

      IF (HBcurve.AND.NewtonRaphson) THEN
        DO p=1,nd
          DO q=1,nd
            JAC(p,q) = JAC(p,q) + IP % s(t) * DetJ * &
              muder/babs*SUM(Agrad*dBasisdx(q,:))*SUM(CONJG(Agrad)*dBasisdx(p,:))
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

    REAL(KIND=dp) :: R(n), R_ip, &
            Inf_ip,Coord(3),Normal(3),mu,u,v

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
    Material => GetMaterial(Parent)
    CALL GetReluctivity(Material,R,n,Element)

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                 IP % W(t), detJ, Basis )

      mu = SUM(Basis(1:n)*R(1:n))

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
 SUBROUTINE GetReluctivity(Material,Acoef,n,Element)
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Material
    INTEGER :: n
    REAL(KIND=dp) :: Acoef(:)
    TYPE(Element_t), OPTIONAL :: Element
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

    Acoef(1:n) = GetReal( Material, 'Relative Permeability', Found, Element )
    IF ( Found ) THEN
      Acoef(1:n) = Avacuum * Acoef(1:n)
    ELSE
      Acoef(1:n) = GetReal( Material, 'Permeability', Found, Element )
    END IF
    IF ( Found ) THEN
      Acoef(1:n) = 1._dp / Acoef(1:n)
    ELSE
      Acoef(1:n) = GetReal( Material, 'Reluctivity', Found, Element )
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
  END IF


END SUBROUTINE Bsolver_init


!------------------------------------------------------------------------------
!> Given the vector potential computes its gradient i.e. the magnetic
!> field intensity.  
!------------------------------------------------------------------------------
SUBROUTINE Bsolver( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------

  USE DefUtils

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
  REAL(KIND=dp) :: at0,at1,at2,CPUTime,RealTime
  
  TYPE(Variable_t), POINTER :: FluxSol, HeatingSol, JouleSol, AzSol
  LOGICAL ::  FoundMortar, CSymmetry, LossEstimation, JouleHeating
  TYPE(Matrix_t),POINTER::A,B,CM,S
  REAL(KIND=dp) :: Omega
  
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
    END IF
  END IF

  !------------------------------------------------------------------------------
  ! In the case of time-harmonic analysis losses may be estimated in terms of B
  !------------------------------------------------------------------------------ 
  LossEstimation = GetLogical(SolverParams,'Loss Estimation',GotIt)
  IF( LossEstimation .AND. FluxDofs /= 4) THEN
    CALL Fatal( 'BSolver', 'Real solution, loss estimation omitted' )
  END IF


  ALLOCATE(ForceVector(SIZE(Solver % Matrix % RHS),TotDOFs))  
  ForceVector = 0.0_dp
  SaveRHS => Solver % Matrix % RHS

  at0 = RealTime()
  CALL BulkAssembly()
  CALL DefaultFinishAssembly()
  at1 = RealTime()
  WRITE(Message,* ) 'Assembly Time: ',at1-at0
  CALL Info( 'BSolver', Message, Level=5 )
!        
!------------------------------------------------------------------------------     

  S=>Solver % Matrix

!  DO i=1,Model % NumberOFBCs
!    j=GetInteger(Model % BCs(i) % Values,'Mortar BC',FoundMortar)
!    IF(FoundMortar) THEN
!      A => PeriodicProjector(Model,Solver % Mesh,i,j,dim,.TRUE.)
!      FoundMortar = ASSOCIATED(A)
!      EXIT
!    END IF
!  END DO

!  IF(FoundMortar) THEN
!    CM=>A
!    ALLOCATE(CM % RHS(CM % NumberOfRows))
!    CM % RHS = 0._dp
!    S % ConstraintMatrix => CM
!    CM % Cols = Solver % Variable % Perm(CM % Cols)
!    CM % Ordered = .FALSE.
!    CALL CRS_SortMatrix(CM,.TRUE.)
!  END IF

   CALL DefaultDirichletBCs()
      
   TotNorm = 0.0_dp
   DO i=1,TotDofs
     Solver % Matrix % RHS => ForceVector(:,i)
     UNorm = DefaultSolve()
     TotNorm = TotNorm + SUM(Solver % Variable % Values**2)
     IF( i <= FluxDofs ) THEN
       FluxSol % Values(i::FluxDofs) = Solver % Variable % Values
     ELSE IF( i == FluxDofs + 1 ) THEN
       JouleSol % Values = Solver % Variable % Values
     ELSE
       HeatingSol % Values = Solver % Variable % Values
     END IF
   END DO
   DEALLOCATE( ForceVector )  

   Solver % Matrix % RHS => SaveRHS
   TotNorm = SQRT(TotNorm)
   Solver % Variable % Norm = Totnorm

!  IF(FoundMortar) THEN
!    CALL FreeMatrix(CM)
!    S % ConstraintMatrix=>NULL()
!  END IF

!------------------------------------------------------------------------------     

  at2 = RealTime()
  WRITE(Message,* ) 'Solution Time: ',at2-at1
  CALL Info( 'BSolver', Message, Level=5 )
  
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
    REAL(KIND=dp) :: weight,coeff,detJ,BAtIp(6),PotAtIp(2),CondAtIp,&
        Omega,TotalHeating, DesiredHeating, HeatingCoeff
    REAL(KIND=dp) :: Freq, FreqPower, FieldPower, ComponentLoss(2), LossCoeff, &
        ValAtIp, TotalLoss, x
    LOGICAL :: Found, SetHeating
    TYPE(ValueList_t), POINTER :: Material

    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:,:)
    REAL(KIND=dp), ALLOCATABLE :: POT(:,:)
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
    REAL(KIND=dp), ALLOCATABLE :: Cond(:)
    REAL(KIND=dp), ALLOCATABLE :: BodyLoss(:)

    SAVE Nodes

    n = 2*MAX(Solver % Mesh % MaxElementDOFs,Solver % Mesh % MaxElementNodes)
    ALLOCATE( STIFF(n,n), FORCE(Totdofs,n) )
    ALLOCATE( POT(2,n), Basis(n), dBasisdx(n,3) )
    
    IF( JouleHeating ) THEN
      ALLOCATE( Cond(n) ) 
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


    DO elem = 1,GetNOFActive()
         
      ! Element information
      ! ---------------------
      Element => GetActiveElement(elem)
      CALL GetElementNodes( Nodes )
      nd = GetElementNOFDOFs()
      n  = GetElementNOFNodes()
      
      ! Integrate local stresses:
      ! -------------------------
      IntegStuff = GaussPoints( Element )
      STIFF  = 0.0_dp
      FORCE  = 0.0_dp

      CALL GetLocalSolution( POT, VarName )

      IF( JouleHeating ) THEN
        Material => GetMaterial()
        Cond(1:n) = GetReal( Material, 'Electric Conductivity', Found, Element)
      END IF

      IF( LossEstimation ) THEN
        BodyId = GetBody() 
        LossCoeff = ListGetFun( Material,'Fourier Loss Coefficient',Freq,Found ) 
      END IF

      DO t=1,IntegStuff % n
        Found = ElementInfo( Element, Nodes, IntegStuff % u(t), &
            IntegStuff % v(t), IntegStuff % w(t), detJ, Basis, dBasisdx )
        
        Weight = IntegStuff % s(t) * detJ
        IF( CSymmetry ) THEN
          x = SUM( Basis(1:n) * Nodes % x(1:n) )
          Weight = Weight * x
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
          CondAtIp = SUM( Basis(1:n) * Cond(1:n) )
          PotAtIp(1) = SUM( POT(1,1:nd) * Basis(1:nd ) )
          PotAtIp(2) = SUM( POT(2,1:nd) * Basis(1:nd ) )
          BAtIp(5) = 0.5_dp * Omega**2 * ( PotAtIp(1)**2 + PotAtIp(2)**2 )  
          BAtIp(6) = CondAtIp * BAtIp(5) 
          TotalHeating = TotalHeating + Weight * BAtIp(6)
        END IF

        IF( LossEstimation ) THEN
          DO i=1,2
            ValAtIP = SUM( BAtIP(2*i-1:2*i) ** 2 )
            Coeff = Weight * LossCoeff * ( Freq ** FreqPower ) * ( ValAtIp ** FieldPower )
            ComponentLoss(i) = ComponentLoss(i) + Coeff
            BodyLoss(BodyId) = BodyLoss(BodyId) + Coeff
          END DO          
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
         CALL Info('FourierLosses', Message, Level=6 )
      END DO

      DEALLOCATE( BodyLoss )
    END IF

    DEALLOCATE( POT, STIFF, FORCE, Basis, dBasisdx )

    IF( JouleHeating ) THEN
      DEALLOCATE( Cond ) 
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE BulkAssembly
!------------------------------------------------------------------------------

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
