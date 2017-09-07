!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
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
! ******************************************************************************
! *
! *  Authors: Olivier Gagliardini             
! *  Email:   gagliar@lgge.obs.ujf-grenoble.fr
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 30. April 2010
! * 
! *****************************************************************************
!> SSolver to inquire the velocity from the SSA solution            
SUBROUTINE SSABasalSolver( Model,Solver,dt,TransientSimulation )
  !------------------------------------------------------------------------------
  !******************************************************************************
  !
  !  Solve the in-plane basal velocity with the SSA solution !
  !  To be computed only at the base. Use then the SSASolver to export verticaly 
  !  the basal velocity and compute the vertical velocity and pressure (if needed)
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
  TYPE(Nodes_t)   :: ElementNodes
  TYPE(Element_t),POINTER :: CurrentElement, Element, ParentElement, BoundaryElement
  TYPE(Matrix_t),POINTER  :: StiffMatrix
  TYPE(ValueList_t), POINTER :: SolverParams, BodyForce, Material, BC
  TYPE(Variable_t), POINTER :: PointerToVariable, ZsSol, ZbSol, VeloSol, Nsol

  LOGICAL :: AllocationsDone = .FALSE., Found, GotIt, CalvingFront, UnFoundFatal=.TRUE.
  LOGICAL :: Newton

  INTEGER :: i, n, m, t, istat, DIM, p, STDOFs, iFriction
  INTEGER :: NonlinearIter, NewtonIter, iter, other_body_id

  INTEGER, POINTER :: Permutation(:), ZsPerm(:), ZbPerm(:), &
       NodeIndexes(:), NPerm(:)

  REAL(KIND=dp), POINTER :: ForceVector(:)
  REAL(KIND=dp), POINTER :: VariableValues(:), Zs(:), Zb(:), Nval(:)

  REAL(KIND=dp) :: UNorm, cn, dd, NonlinearTol, NewtonTol, MinSRInv, MinH, rhow, sealevel, &
       PrevUNorm, relativeChange, minv, fm, PostPeak, MinN

  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:), &
       NodalGravity(:), NodalViscosity(:), NodalDensity(:), &
       NodalZs(:), NodalZb(:), NodalU(:), NodalV(:),  &
       NodalBeta(:), NodalLinVelo(:), NodalC(:), NodalN(:)


  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, Friction, ZsName, ZbName
#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: at, at0
#else
  REAL(KIND=dp) :: at, at0, CPUTime, RealTime
#endif     

  SAVE rhow,sealevel
  SAVE STIFF, LOAD, FORCE, AllocationsDone, DIM, SolverName, ElementNodes
  SAVE NodalGravity, NodalViscosity, NodalDensity, &
       NodalZs, NodalZb,   &
       NodalU, NodalV, NodeIndexes, &
       NodalBeta, NodalLinVelo, NodalC, NodalN

  !------------------------------------------------------------------------------
  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values
  STDOFs = PointerToVariable % DOFs 
  WRITE(SolverName, '(A)') 'SSASolver-SSABasalSolver'

  !------------------------------------------------------------------------------
  !    Get variables needed for solution
  !------------------------------------------------------------------------------
  DIM = CoordinateSystemDimension()

  ! IF DIM = STDOFs+1  Normal-Tangential can not be used => trick temporary set Model Dimension to STDOFs
  IF (DIM.eq.(STDOFs+1)) CurrentModel % Dimension = STDOFs

  SolverParams => GetSolverParams()

  ZbName = GetString(SolverParams, 'Bottom Surface Name', GotIt)
  IF (GotIt) THEN
    CALL INFO(SolverName, 'Bottom Surface Name found', level=4)            
  ELSE
    CALL INFO(SolverName, 'Bottom Surface Name not found - using default Zb', level=1) 
    WRITE(ZbName,'(A)') 'Zb'
  END IF
  ZbSol => VariableGet( Solver % Mesh % Variables, ZbName,UnFoundFatal=UnFoundFatal)
  Zb => ZbSol % Values
  ZbPerm => ZbSol % Perm

  ZsName = GetString(SolverParams, 'Top Surface Name', GotIt)
  IF (GotIt) THEN
    CALL INFO(SolverName, 'Top Surface Name found', level=4)            
  ELSE
    CALL INFO(SolverName, 'Top Surface Name not found - using default Zs', level=1) 
    WRITE(ZsName,'(A)') 'Zs'
  END IF
  ZsSol => VariableGet( Solver % Mesh % Variables, ZsName,UnFoundFatal=UnFoundFatal)
  Zs => ZsSol % Values
  ZsPerm => ZsSol % Perm
  NSol => VariableGet( Solver % Mesh % Variables, 'Effective Pressure' )
  IF (ASSOCIATED(NSol)) THEN
    Nval => NSol % Values
    NPerm => NSol % Perm
  END IF
  !--------------------------------------------------------------
  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed  ) THEN

    ! Get some constants
    rhow = GetConstReal( Model % Constants, 'Water Density', Found )
    If (.NOT.Found) Then
      WRITE(Message,'(A)') 'Constant Water Density not found. &
           &Setting to 1.03225e-18'
      CALL INFO(SolverName, Message, level=20)
      rhow = 1.03225e-18_dp
    End if

    sealevel = GetCReal( Model % Constants, 'Sea Level', Found )
    If (.NOT.Found) Then
      WRITE(Message,'(A)') 'Constant >Sea Level< not found. &
           &Setting to 0.0'
      CALL INFO(SolverName, Message, level=20)
      sealevel=0.0_dp
    End if

    ! Allocate

    N = Model % MaxElementNodes
    M = Model % Mesh % NumberOfNodes
    IF (AllocationsDone) DEALLOCATE(FORCE, LOAD, STIFF, NodalGravity, &
         NodalViscosity, NodalDensity,  &
         NodalZb, NodalZs,  NodalU, NodalV, &
         NodalBeta, NodalLinVelo, NodalC, NodalN, &
         ElementNodes % x, &
         ElementNodes % y, ElementNodes % z )

    ALLOCATE( FORCE(STDOFs*N), LOAD(N), STIFF(STDOFs*N,STDOFs*N), &
         NodalGravity(N), NodalDensity(N), NodalViscosity(N), &
         NodalZb(N), NodalZs(N), NodalU(N), NodalV(N), &
         NodalBeta(N), NodalLinVelo(N), NodalC(N), NodalN(N), &
         ElementNodes % x(N), ElementNodes % y(N), ElementNodes % z(N), &
         STAT=istat )
    IF ( istat /= 0 ) THEN
      CALL Fatal( SolverName, 'Memory allocation error.' )
    END IF

    AllocationsDone = .TRUE.
    CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
  END IF

  StiffMatrix => Solver % Matrix
  ForceVector => StiffMatrix % RHS

  !------------------------------------------------------------------------------
  !    Do some additional initialization, and go for it
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  NonlinearTol = GetConstReal( Solver % Values, &
       'Nonlinear System Convergence Tolerance' )

  NonlinearIter = GetInteger( Solver % Values, &
       'Nonlinear System Max Iterations',GotIt )

  IF ( .NOT.GotIt ) NonlinearIter = 1

  NewtonTol = ListGetConstReal( Solver % Values, &
       'Nonlinear System Newton After Tolerance', minv=0.0d0 )

  NewtonIter = ListGetInteger( Solver % Values, &
       'Nonlinear System Newton After Iterations', GotIt )
  if (.NOT.Gotit) NewtonIter = NonlinearIter + 1

  Newton=.False.
  !------------------------------------------------------------------------------
  DO iter=1,NonlinearIter

    at  = CPUTime()
    at0 = RealTime()

    CALL Info( SolverName, ' ', Level=4 )
    CALL Info( SolverName, ' ', Level=4 )
    CALL Info( SolverName, &
         '-------------------------------------',Level=4 )
    WRITE( Message, * ) 'SSA BASAL VELOCITY NON-LINEAR ITERATION', iter
    CALL Info( SolverName, Message, Level=4 )
    If (Newton) Then
      WRITE( Message, * ) 'Newton linearisation is used'
      CALL Info( SolverName, Message, Level=4 )
    Endif
    CALL Info( SolverName, ' ', Level=4 )
    CALL Info( SolverName, &
         '-------------------------------------',Level=4 )
    CALL Info( SolverName, ' ', Level=4 )


    !Initialize the system and do the assembly:
    !------------------------------------------
    CALL DefaultInitialize()

    ! bulk assembly
    DO t=1,Solver % NumberOfActiveElements
      Element => GetActiveElement(t)
      IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
      n = GetElementNOFNodes()

      NodeIndexes => Element % NodeIndexes

      ! set coords of highest occurring dimension to zero (to get correct path element)
      !-------------------------------------------------------------------------------
      ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
      IF (STDOFs == 1) THEN !1D SSA
        ElementNodes % y(1:n) = 0.0_dp
        ElementNodes % z(1:n) = 0.0_dp
      ELSE IF (STDOFs == 2) THEN !2D SSA
        ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
        ElementNodes % z(1:n) = 0.0_dp
      ELSE
        WRITE(Message,'(a,i1,a)')&
             'It is not possible to compute SSA problems with DOFs=',&
             STDOFs, ' . Aborting'
        CALL Fatal( SolverName, Message)
        STOP
      END IF

      ! Read the gravity in the Body Force Section 
      BodyForce => GetBodyForce()
      NodalGravity = 0.0_dp
      IF ( ASSOCIATED( BodyForce ) ) THEN
        IF (STDOFs==1) THEN 
          NodalGravity(1:n) = ListGetReal( &
               BodyForce, 'Flow BodyForce 2', n, NodeIndexes, Found)
        ELSE 
          NodalGravity(1:n) = ListGetReal( &
               BodyForce, 'Flow BodyForce 3', n, NodeIndexes, Found)
        END IF
      END IF

      ! Read the Viscosity eta, density, and exponent m in MMaterial Section
      ! Same definition as NS Solver in Elmer - n=1/m , A = 1/ (2 eta^n) 
      Material => GetMaterial(Element)

      cn = ListGetConstReal( Material, 'Viscosity Exponent',Found)
      MinSRInv = ListGetConstReal( Material, 'Critical Shear Rate',Found)
      MinH = ListGetConstReal( Material, 'SSA Critical Thickness',Found)
      If (.NOT.Found) MinH=EPSILON(MinH)


      NodalDensity=0.0_dp
      NodalDensity(1:n) = ListGetReal( Material, 'SSA Mean Density',n,NodeIndexes,Found,&
           UnFoundFatal=UnFoundFatal)

      NodalViscosity=0.0_dp
      NodalViscosity(1:n) = ListGetReal( Material, 'SSA Mean Viscosity',n, NodeIndexes,Found,&
           UnFoundFatal=UnFoundFatal)

      Friction = GetString(Material, 'SSA Friction Law', Found) 
      IF (.NOT.Found) &
           CALL FATAL(SolverName,'Could not find Material keyword >SSA Friction Law<')

      SELECT CASE(Friction) 
      CASE('linear')
        iFriction = 1
        fm = 1.0_dp
      CASE('weertman')
        iFriction = 2
      CASE('coulomb')
        iFriction = 3
      CASE DEFAULT
        CALL FATAL(SolverName,'Friction should be linear, Weertman or Coulomb')
      END SELECT

      ! for all friction law
      NodalBeta = 0.0_dp
      NodalBeta(1:n) = ListGetReal( &
           Material, 'SSA Friction Parameter', n, NodeIndexes(1:n), Found,&
           UnFoundFatal=UnFoundFatal)

      ! for Weertman and Coulomb friction
      IF (iFriction > 1) THEN
        fm = ListGetConstReal( Material, 'SSA Friction Exponent', Found , UnFoundFatal=UnFoundFatal)

        MinN = ListGetConstReal( Material, 'SSA Min Effective Pressure', Found, UnFoundFatal=UnFoundFatal)
        !Previous default value: MinN = 1.0e-6_dp

        NodalLinVelo = 0.0_dp
        NodalLinVelo(1:n) = ListGetReal( &
             Material, 'SSA Friction Linear Velocity', n, NodeIndexes(1:n), Found,&
             UnFoundFatal=UnFoundFatal)
      END IF

      ! only for Coulomb friction
      IF (iFriction > 2) THEN
        PostPeak = ListGetConstReal( Material, 'SSA Friction Post-Peak', Found, UnFoundFatal=UnFoundFatal )

        NodalC = 0.0_dp
        NodalC(1:n) = ListGetReal( &
             Material, 'SSA Friction Maximum Value', n, NodeIndexes(1:n), Found,&
             UnFoundFatal=UnFoundFatal)

        ! Get the effective pressure
        IF (ASSOCIATED(NSol)) THEN
          NodalN(1:n) = Nval(NPerm(NodeIndexes(1:n)))
        ELSE
          CALL FATAL(SolverName,'Could not find variable >Effective Pressure<')
        END IF
      END IF


      ! Get the Nodal value of Zb and Zs
      NodalZb(1:n) = Zb(ZbPerm(NodeIndexes(1:n)))
      NodalZs(1:n) = Zs(ZsPerm(NodeIndexes(1:n)))

      ! Previous Velocity 
      NodalU(1:n) = VariableValues(STDOFs*(Permutation(NodeIndexes(1:n))-1)+1)
      NodalV = 0.0_dp
      IF (STDOFs.EQ.2) NodalV(1:n) = VariableValues(STDOFs*(Permutation(NodeIndexes(1:n))-1)+2)


      CALL LocalMatrixUVSSA (  STIFF, FORCE, Element, n, ElementNodes, NodalGravity, &
           NodalDensity, NodalViscosity, NodalZb, NodalZs, NodalU, NodalV, &
           iFriction, NodalBeta, fm, NodalLinVelo, PostPeak, NodalC, NodalN, &
           cn, MinSRInv, MinH , STDOFs, Newton)

      CALL DefaultUpdateEquations( STIFF, FORCE )

    END DO
    CALL DefaultFinishBulkAssembly()

    !  
    ! Neumann condition
    !
    DO t=1,GetNOFBoundaryElements()


      BoundaryElement => GetBoundaryElement(t)
      IF (STDOFS.NE.1) then
        IF ( .NOT. ActiveBoundaryElement() ) CYCLE
      END IF
      IF ( GetElementFamily() == 1 ) CYCLE

      NodeIndexes => BoundaryElement % NodeIndexes
      IF (ParEnv % myPe .NE. BoundaryElement % partIndex) CYCLE

      n = GetElementNOFNodes()
      FORCE = 0.0_dp
      STIFF = 0.0_dp

      ! set coords of highest occurring dimension to zero (to get correct path element)
      !-------------------------------------------------------------------------------
      ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
      IF (STDOFs == 1) THEN
        ElementNodes % y(1:n) = 0.0_dp
        ElementNodes % z(1:n) = 0.0_dp
      ELSE IF (STDOFs == 2) THEN
        ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
        ElementNodes % z(1:n) = 0.0_dp
      ELSE
        WRITE(Message,'(a,i1,a)')&
             'It is not possible to compute SSA with SSA var DOFs=',&
             STDOFs, '. Aborting'
        CALL Fatal( SolverName, Message)
        STOP
      END IF


      BC => GetBC()
      IF (.NOT.ASSOCIATED( BC ) ) CYCLE

      ! Find the nodes for which 'Calving Front' = True             
      CalvingFront=.False. 
      CalvingFront = ListGetLogical( BC, 'Calving Front', GotIt )
      IF (CalvingFront) THEN
        NodalZs(1:n) = Zs(ZsPerm(NodeIndexes(1:n)))
        NodalZb(1:n) = Zb(ZbPerm(NodeIndexes(1:n)))
        ! Need to access Parent Element to get Material properties
        other_body_id = BoundaryElement % BoundaryInfo % outbody
        IF (other_body_id < 1) THEN ! only one body in calculation
          ParentElement => BoundaryElement % BoundaryInfo % Right
          IF ( .NOT. ASSOCIATED(ParentElement) ) ParentElement => BoundaryElement % BoundaryInfo % Left
        ELSE ! we are dealing with a body-body boundary and assume that the normal is pointing outwards
          ParentElement => BoundaryElement %  BoundaryInfo % Right
          IF (ParentElement % BodyId == other_body_id) ParentElement =>  BoundaryElement % BoundaryInfo % Left
        END IF

        ! Read Density in the Material Section
        Material => GetMaterial(ParentElement)

        NodalDensity=0.0_dp
        NodalDensity(1:n) = ListGetReal( Material, 'SSA Mean Density',n, NodeIndexes,Found,&
             UnFoundFatal=UnFoundFatal)

        MinH = ListGetConstReal( Material, 'SSA Critical Thickness',Found)
        If (.NOT.Found) MinH=EPSILON(MinH)

        ! Read the gravity in the Body Force Section 
        BodyForce => GetBodyForce(ParentElement)
        NodalGravity = 0.0_dp
        IF ( ASSOCIATED( BodyForce ) ) THEN
          IF (STDOFs==1) THEN 
            NodalGravity(1:n) = ListGetReal( &
                 BodyForce, 'Flow BodyForce 2', n, NodeIndexes, Found)
          ELSE 
            NodalGravity(1:n) = ListGetReal( &
                 BodyForce, 'Flow BodyForce 3', n, NodeIndexes, Found)
          END IF
        END IF

        CALL LocalMatrixBCSSA(  STIFF, FORCE, BoundaryElement, n, ElementNodes,&
             NodalDensity, NodalGravity, NodalZb, NodalZs, rhow, sealevel,MinH )
        CALL DefaultUpdateEquations( STIFF, FORCE )
      END IF
    END DO

    CALL DefaultFinishAssembly()

    ! Dirichlet 
    CALL DefaultDirichletBCs()


    !------------------------------------------------------------------------------
    !     Solve the system and check for convergence
    !------------------------------------------------------------------------------
    PrevUNorm = UNorm

    UNorm = DefaultSolve()


    RelativeChange = Solver % Variable % NonlinChange
    !IF ( PrevUNorm + UNorm /= 0.0d0 ) THEN
    !   RelativeChange = 2.0d0 * ABS( PrevUNorm - UNorm) / ( PrevUnorm + UNorm)
    !ELSE
    !   RelativeChange = 0.0d0
    !END IF

    WRITE( Message, * ) 'Result Norm   : ', UNorm, PrevUNorm
    CALL Info(SolverName, Message, Level=4 )
    WRITE( Message, * ) 'Relative Change : ', RelativeChange
    CALL Info(SolverName, Message, Level=4 )


    IF ( RelativeChange < NewtonTol .OR. &
         iter > NewtonIter ) Newton = .TRUE.

    !------------------------------------------------------------------------------
    IF ( RelativeChange < NonLinearTol ) EXIT
    !------------------------------------------------------------------------------

  END DO ! Loop Non-Linear Iterations

!!! reset Model Dimension to dim
  IF (DIM.eq.(STDOFs+1)) CurrentModel % Dimension = DIM

CONTAINS

  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixUVSSA(  STIFF, FORCE, Element, n, Nodes, gravity, &
       Density, Viscosity, LocalZb, LocalZs, LocalU, LocalV, &
       Friction, LocalBeta, fm, LocalLinVelo, fq, LocalC, LocalN, &
       cm, MinSRInv, MinH, STDOFs , Newton )
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), gravity(:), Density(:), &
         Viscosity(:), LocalZb(:), LocalZs(:), &
         LocalU(:), LocalV(:) , LocalBeta(:), &
         LocalLinVelo(:), LocalC(:), LocalN(:)
    INTEGER :: n, cp , STDOFs, Friction
    REAL(KIND=dp) :: cm, fm, fq
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: Newton
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), detJ 
    REAL(KIND=dp) :: g, rho, eta, h, dhdx, dhdy , muder
    REAL(KIND=dp) :: beta, LinVelo, fC, fN, Velo(2), ub, alpha, fB
    REAL(KIND=dp) :: gradS(2), A(2,2), StrainA(2,2), StrainB(2,2), Exx, Eyy, Exy, Ezz, Ee, MinSRInv ,MinH                           
    REAL(KIND=dp) :: Jac(2*n,2*n), SOL(2*n), Slip, Slip2
    LOGICAL :: Stat, NewtonLin, fNewtonLIn
    INTEGER :: i, j, t, p, q , dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    !------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()

    STIFF = 0.0_dp
    FORCE = 0.0_dp
    Jac=0.0_dp

    ! Use Newton Linearisation
    NewtonLin = (Newton.AND.(cm.NE.1.0_dp))
    fNewtonLin = (Newton.AND.(fm.NE.1.0_dp))

    IP = GaussPoints( Element )
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

      ! Needed Integration Point value

      g = ABS(SUM( Gravity(1:n) * Basis(1:n) ))
      rho = SUM( Density(1:n) * Basis(1:n) )
      eta = SUM( Viscosity(1:n) * Basis(1:n) )
      gradS = 0.0_dp
      gradS(1) = SUM( LocalZs(1:n) * dBasisdx(1:n,1) )
      IF (STDOFs == 2) gradS(2) = SUM( LocalZs(1:n) * dBasisdx(1:n,2) )
      h = SUM( (LocalZs(1:n)-LocalZb(1:n)) * Basis(1:n) )
      h=max(h,MinH)

      beta = SUM( LocalBeta(1:n) * Basis(1:n) )
      IF (iFriction > 1) THEN
        LinVelo = SUM( LocalLinVelo(1:n) * Basis(1:n) )
        IF ((iFriction == 2).AND.(fm==1.0_dp)) iFriction=1
        Velo = 0.0_dp
        Velo(1) = SUM(LocalU(1:n) * Basis(1:n))
        IF (STDOFs == 2) Velo(2) = SUM(LocalV(1:n) * Basis(1:n))
        ub = SQRT(Velo(1)*Velo(1)+Velo(2)*Velo(2))
        Slip2=1.0_dp
        IF (ub < LinVelo) then 
          ub = LinVelo
          Slip2=0.0_dp
        ENDIF
      END IF

      IF (iFriction==3) THEN
        fC = SUM( LocalC(1:n) * Basis(1:n) )
        fN = SUM( LocalN(1:n) * Basis(1:n) )
        ! Effective pressure should be >0 (for the friction law)
        fN = MAX(fN, MinN)
      END IF

      IF (iFriction==1) THEN
        Slip = beta
        fNewtonLin = .FALSE.
      ELSE IF (iFriction==2) THEN
        Slip = beta * ub**(fm-1.0_dp) 
        Slip2 = Slip2*Slip*(fm-1.0_dp)/(ub*ub)
      ELSE IF (iFriction==3) THEN
        IF (PostPeak.NE.1.0_dp) THEN
          alpha = (PostPeak-1.0_dp)**(PostPeak-1.0_dp) / PostPeak**PostPeak
        ELSE
          alpha = 1.0_dp
        END IF
        fB = alpha * (beta / (fC*fN))**(PostPeak/fm)
        Slip = beta * ub**(fm-1.0_dp) / (1.0_dp + fB * ub**PostPeak)**fm
        Slip2 = Slip2 * Slip * ((fm-1.0_dp) / (ub*ub) - &
             fm*PostPeak*fB*ub**(PostPeak-2.0_dp)/(1.0_dp+fB*ub**PostPeak))
      END IF

      !------------------------------------------------------------------------------
      ! In the non-linear case, effective viscosity       
      IF (cm.NE.1.0_dp) THEN
        Exx = SUM(LocalU(1:n)*dBasisdx(1:n,1))
        Eyy = 0.0_dp
        Exy = 0.0_dp
        IF (STDOFs.EQ.2) THEN
          Eyy = SUM(LocalV(1:n)*dBasisdx(1:n,2))
          Ezz = -Exx - Eyy
          Exy = SUM(LocalU(1:n)*dBasisdx(1:n,2))
          Exy = 0.5_dp*(Exy + SUM(LocalV(1:n)*dBasisdx(1:n,1)))
          Ee = 0.5_dp*(Exx*Exx + Eyy*Eyy + Ezz*Ezz) + Exy*Exy
        ELSE
          Ee = Exx*Exx
        END IF
        muder = eta * 0.5_dp * (2.0_dp**cm) * ((cm-1.0_dp)/2.0_dp) *  Ee**((cm-1.0_dp)/2.0_dp - 1.0_dp)
        IF (sqrt(Ee) < MinSRInv) THEN
          Ee = MinSRInv*MinSRInv
          muder = 0.0_dp
        END IF
        eta = eta * 0.5_dp * (2.0_dp**cm) * Ee**((cm-1.0_dp)/2.0_dp)
      END IF

      StrainA=0.0_dp
      StrainB=0.0_dp
      IF (NewtonLin) THEN
        StrainA(1,1)=SUM(2.0_dp*dBasisdx(1:n,1)*LocalU(1:n))

        IF (STDOFs.EQ.2) THEN
          StrainB(1,1)=SUM(0.5_dp*dBasisdx(1:n,2)*LocalU(1:n))

          StrainA(1,2)=SUM(dBasisdx(1:n,2)*LocalV(1:n))
          StrainB(1,2)=SUM(0.5_dp*dBasisdx(1:n,1)*LocalV(1:n))

          StrainA(2,1)=SUM(dBasisdx(1:n,1)*LocalU(1:n))
          StrainB(2,1)=SUM(0.5_dp*dBasisdx(1:n,2)*LocalU(1:n))

          StrainA(2,2)=SUM(2.0_dp*dBasisdx(1:n,2)*LocalV(1:n))
          StrainB(2,2)=SUM(0.5_dp*dBasisdx(1:n,1)*LocalV(1:n))
        END IF
      END IF

      A = 0.0_dp
      DO p=1,n
        DO q=1,n
          A(1,1) = 2.0_dp*dBasisdx(q,1)*dBasisdx(p,1)  
          IF (STDOFs.EQ.2) THEN
            A(1,1) = A(1,1) + 0.5_dp*dBasisdx(q,2)*dBasisdx(p,2)
            A(1,2) = dBasisdx(q,2)*dBasisdx(p,1) + &
                 0.5_dp*dBasisdx(q,1)*dBasisdx(p,2)
            A(2,1) = dBasisdx(q,1)*dBasisdx(p,2) + &
                 0.5_dp*dBasisdx(q,2)*dBasisdx(p,1)
            A(2,2) = 2.0*dBasisdx(q,2)*dBasisdx(p,2) +&
                 0.5_dp*dBasisdx(q,1)*dBasisdx(p,1)  
          END IF
          A = 2.0_dp * h * eta * A
          DO i=1,STDOFs
            STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+i) = STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+i) +&
                 Slip * Basis(q) * Basis(p) * IP % S(t) * detJ
            DO j=1,STDOFs
              STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+j) = STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+j) +& 
                   A(i,j) * IP % S(t) * detJ 
            END DO
          END DO

          IF ((fNewtonLin).AND.(iFriction > 1)) THEN
            DO i=1,STDOFs
              Do j=1,STDOFs
                STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+j) = STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+j) +&
                     Slip2 * Velo(i) * Velo(j) * Basis(q) * Basis(p) * IP % S(t) * detJ
              End do
            END DO
          END IF

          IF (NewtonLin) then
            IF (STDOFs.EQ.1) THEN
              Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+1) = Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+1) +&
                   IP % S(t) * detJ * 2.0_dp * h * StrainA(1,1)*dBasisdx(p,1) * &
                   muder * 2.0_dp * Exx*dBasisdx(q,1) 
            ELSE IF (STDOFs.EQ.2) THEN
              Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+1) = Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+1) +&
                   IP % S(t) * detJ * 2.0_dp * h * ((StrainA(1,1)+StrainA(1,2))*dBasisdx(p,1)+ &
                   (StrainB(1,1)+StrainB(1,2))*dBasisdx(p,2)) * muder *((2.0_dp*Exx+Eyy)*dBasisdx(q,1)+Exy*dBasisdx(q,2)) 

              Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+2) = Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+2) +&
                   IP % S(t) * detJ * 2.0_dp * h * ((StrainA(1,1)+StrainA(1,2))*dBasisdx(p,1)+ & 
                   (StrainB(1,1)+StrainB(1,2))*dBasisdx(p,2)) * muder *((2.0_dp*Eyy+Exx)*dBasisdx(q,2)+Exy*dBasisdx(q,1)) 

              Jac((STDOFs)*(p-1)+2,(STDOFs)*(q-1)+1) = Jac((STDOFs)*(p-1)+2,(STDOFs)*(q-1)+1) +&
                   IP % S(t) * detJ * 2.0_dp * h * ((StrainA(2,1)+StrainA(2,2))*dBasisdx(p,2)+ & 
                   (StrainB(2,1)+StrainB(2,2))*dBasisdx(p,1)) * muder *((2.0_dp*Exx+Eyy)*dBasisdx(q,1)+Exy*dBasisdx(q,2)) 

              Jac((STDOFs)*(p-1)+2,(STDOFs)*(q-1)+2) = Jac((STDOFs)*(p-1)+2,(STDOFs)*(q-1)+2) +&
                   IP % S(t) * detJ * 2.0_dp * h * ((StrainA(2,1)+StrainA(2,2))*dBasisdx(p,2)+ &
                   (StrainB(2,1)+StrainB(2,2))*dBasisdx(p,1)) * muder *((2.0_dp*Eyy+Exx)*dBasisdx(q,2)+Exy*dBasisdx(q,1)) 
            END IF
          END IF

        END DO

        DO i=1,STDOFs
          FORCE((STDOFs)*(p-1)+i) =   FORCE((STDOFs)*(p-1)+i) - &   
               rho*g*h*gradS(i) * IP % s(t) * detJ * Basis(p) 
        END DO

        IF ((fNewtonLin).AND.(iFriction>1)) THEN
          DO i=1,STDOFs
            FORCE((STDOFs)*(p-1)+i) =   FORCE((STDOFs)*(p-1)+i) + &   
                 Slip2 * Velo(i) * ub * ub * IP % s(t) * detJ * Basis(p) 
          END DO
        END IF

      END DO
    END DO

    IF (NewtonLin) THEN
      SOL(1:STDOFs*n:STDOFs)=LocalU(1:n)
      IF (STDOFs.EQ.2) SOL(2:STDOFs*n:STDOFs)=LocalV(1:n)

      STIFF(1:STDOFs*n,1:STDOFs*n) = STIFF(1:STDOFs*n,1:STDOFs*n) + &
           Jac(1:STDOFs*n,1:STDOFs*n)
      FORCE(1:STDOFs*n) = FORCE(1:STDOFs*n) + &
           MATMUL(Jac(1:STDOFs*n,1:STDOFs*n),SOL(1:STDOFs*n))
    END IF
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixUVSSA
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBCSSA(  STIFF, FORCE, Element, n, ENodes, Density, & 
       Gravity, LocalZb, LocalZs, rhow, sealevel, MinH)
    !------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) ::  ENodes
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:),  density(:), Gravity(:), LocalZb(:),&
         LocalZs(:),rhow, sealevel, MinH
    INTEGER :: n
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3), &
         DetJ,Normal(3), rhoi, g, alpha, h, h_im,norm
    LOGICAL :: Stat
    INTEGER :: t, i
    TYPE(GaussIntegrationPoints_t) :: IP

    !------------------------------------------------------------------------------
    STIFF = 0.0_dp
    FORCE = 0.0_dp

    ! The front force is a concentrated nodal force in 1D-SSA and
    ! a force distributed along a line in 2D-SSA    

    ! 1D-SSA Case : concentrated force at each nodes
    IF (STDOFs==1) THEN  !1D SSA but should be 2D problem (does elmer work in 1D?)
      DO i = 1, n
        g = ABS( Gravity(i) )
        rhoi = Density(i)
        h = LocalZs(i)-LocalZb(i) 
        h = max(h,MinH)
        h_im = max(0.0_dp,sealevel-LocalZb(i))
        alpha=0.5_dp * g * (rhoi * h*h - rhow * h_im*h_im)
        FORCE(i) = FORCE(i) + alpha
      END DO

      ! 2D-SSA Case : force distributed along the line       
      ! This will work in DIM=3D only if working with Extruded Mesh and Preserve
      ! Baseline as been set to True to keep the 1D-BC 
    ELSE IF (STDOFs==2) THEN

      IP = GaussPoints( Element )
      DO t=1,IP % n
        stat = ElementInfo( Element, ENodes, IP % U(t), IP % V(t), &
             IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

        g = ABS(SUM( Gravity(1:n) * Basis(1:n) ))
        rhoi = SUM( Density(1:n) * Basis(1:n) )
        h = SUM( (LocalZs(1:n)-LocalZb(1:n)) * Basis(1:n))
        h_im = max(0.0_dp , SUM( (sealevel-LocalZb(1:n)) * Basis(1:n)) )
        alpha=0.5_dp * g * (rhoi * h*h - rhow * h_im*h_im)

        ! Normal in the (x,y) plane
        Normal = NormalVector( Element, ENodes, IP % U(t), IP % V(t), .TRUE.)
        norm=SQRT(normal(1)*normal(1) +normal(2)*normal(2))
        Normal(1) = Normal(1)/norm
        Normal(2) = Normal(2)/norm

        DO p=1,n
          DO i=1,STDOFs
            FORCE(STDOFs*(p-1)+i) =   FORCE(STDOFs*(p-1)+i) +&   
                 alpha * Normal(i) * IP % s(t) * detJ * Basis(p) 
          END DO
        END DO
      END DO

    ELSE   

      CALL FATAL('SSASolver-SSABasalSolver','Do not work for STDOFs <> 1 or 2')

    END IF
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBCSSA


  !------------------------------------------------------------------------------
END SUBROUTINE SSABasalSolver
!------------------------------------------------------------------------------

! *****************************************************************************
!>   Compute the depth integrated viscosity = sum_zb^zs eta dz
!>     and the depth integrated density = sum_zb^zs rho dz
SUBROUTINE GetMeanValueSolver( Model,Solver,dt,TransientSimulation )
  !------------------------------------------------------------------------------
  !******************************************************************************
  !
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
  TYPE(Element_t),POINTER :: CurrentElement, Element, ParentElement, &
       BoundaryElement
  TYPE(Matrix_t),POINTER  :: StiffMatrix
  TYPE(ValueList_t), POINTER :: SolverParams, BodyForce, Material
  TYPE(Variable_t), POINTER :: PointerToVariable, IntViscoSol, IntDensSol,&
       DepthSol  

  LOGICAL :: AllocationsDone = .FALSE., Found,UnFoundFatal=.TRUE.

  INTEGER :: i, n, m, t, istat, DIM, COMP, other_body_id   
  INTEGER, POINTER :: Permutation(:), NodeIndexes(:), IntViscoPerm(:),&
       IntDensPerm(:), DepthPerm(:) 

  REAL(KIND=dp), POINTER :: ForceVector(:)
  REAL(KIND=dp), POINTER :: VariableValues(:), IntVisco(:), IntDens(:), Depth(:)
  REAL(KIND=dp) :: Norm, cn, dd 

  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:), &
       NodalVar(:) 

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName

  SAVE STIFF, LOAD, FORCE, AllocationsDone, DIM, SolverName
  SAVE NodalVar 
  !------------------------------------------------------------------------------
  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values
  WRITE(SolverName, '(A)') 'SSASolver-IntValue'

  IntViscoSol => VariableGet( Solver % Mesh % Variables, 'Mean Viscosity',UnFoundFatal=UnFoundFatal)
  IntVisco => IntViscoSol % Values
  IntViscoPerm => IntViscoSol % Perm

  IntDensSol => VariableGet( Solver % Mesh % Variables, 'Mean Density',UnFoundFatal=UnFoundFatal)
  IntDens => IntDensSol % Values
  IntDensPerm => IntDensSol % Perm

  DepthSol => VariableGet( Solver % Mesh % Variables, 'Depth',UnFoundFatal=UnFoundFatal)
  Depth => DepthSol % Values
  DepthPerm => DepthSol % Perm

  !--------------------------------------------------------------
  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed  ) THEN
    N = Solver % Mesh % MaxElementNodes ! just big enough for elemental arrays
    M = Model % Mesh % NumberOfNodes
    IF (AllocationsDone) DEALLOCATE(FORCE, LOAD, STIFF, NodalVar) 

    ALLOCATE( FORCE(N), LOAD(N), STIFF(N,N), NodalVar(N), &
         STAT=istat )
    IF ( istat /= 0 ) THEN
      CALL Fatal( SolverName, 'Memory allocation error.' )
    END IF
    AllocationsDone = .TRUE.
    CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
  END IF

  StiffMatrix => Solver % Matrix
  ForceVector => StiffMatrix % RHS

  ! Loop for viscosity and density
  DO COMP=1, 2
    ! No non-linear iteration, no time dependency  
    VariableValues = 0.0d0
    Norm = Solver % Variable % Norm

    !Initialize the system and do the assembly:
    !------------------------------------------
    CALL DefaultInitialize()
    ! bulk assembly
    DO t=1,Solver % NumberOfActiveElements
      Element => GetActiveElement(t)
      IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
      n = GetElementNOFNodes()

      NodeIndexes => Element % NodeIndexes
      Material => GetMaterial(Element)

      IF (COMP==1) THEN
        ! Read the Viscosity eta, 
        ! Same definition as NS Solver in Elmer - n=1/m , A = 1/ (2 eta^n) 
        NodalVar = 0.0D0
        NodalVar(1:n) = ListGetReal( &
             Material, 'Viscosity', n, NodeIndexes, Found )
      ELSE IF (COMP==2) THEN
        NodalVar = 0.0D0
        NodalVar(1:n) = ListGetReal( &
             Material, 'Density', n, NodeIndexes, Found )
      END IF

      CALL LocalMatrix (  STIFF, FORCE, Element, n, NodalVar )
      CALL DefaultUpdateEquations( STIFF, FORCE )
    END DO

    ! Neumann conditions 
    DO t=1,Solver % Mesh % NUmberOfBoundaryElements
      BoundaryElement => GetBoundaryElement(t)
      IF ( GetElementFamily() == 1 ) CYCLE
      NodeIndexes => BoundaryElement % NodeIndexes
      IF (ParEnv % myPe .NE. BoundaryElement % partIndex) CYCLE
      n = GetElementNOFNodes()

      ! Find the Parent element     
      other_body_id = BoundaryElement % BoundaryInfo % outbody
      IF (other_body_id < 1) THEN ! only one body in calculation
        ParentElement => BoundaryElement % BoundaryInfo % Right
        IF ( .NOT. ASSOCIATED(ParentElement) ) ParentElement => BoundaryElement % BoundaryInfo % Left
      ELSE ! we are dealing with a body-body boundary and assume that the normal is pointing outwards
        ParentElement => BoundaryElement % BoundaryInfo % Right
        IF (ParentElement % BodyId == other_body_id) ParentElement => BoundaryElement % BoundaryInfo % Left
      END IF

      Material => GetMaterial(ParentElement)

      IF (COMP==1) THEN
        ! Read the Viscosity eta, 
        ! Same definition as NS Solver in Elmer - n=1/m , A = 1/ (2 eta^n) 
        NodalVar = 0.0D0
        NodalVar(1:n) = ListGetReal( &
             Material, 'Viscosity', n, NodeIndexes, Found )
      ELSE IF (COMP==2) THEN
        NodalVar = 0.0D0
        NodalVar(1:n) = ListGetReal( &
             Material, 'Density', n, NodeIndexes, Found )
      END IF
      CALL LocalMatrixBC(  STIFF, FORCE, BoundaryElement, n, NodalVar)
      CALL DefaultUpdateEquations( STIFF, FORCE )
    END DO

    CALL DefaultFinishAssembly()
    ! Dirichlet 
    IF (COMP==1) THEN
      CALL SetDirichletBoundaries( Model, StiffMatrix, ForceVector, &
           'Mean Viscosity', 1,1, Permutation )
    ELSE
      CALL SetDirichletBoundaries( Model, StiffMatrix, ForceVector, &
           'Mean Density', 1,1, Permutation )
    END IF
    Norm = DefaultSolve()

    ! Save the solution on the right variable
    IF (COMP==1) THEN
      DO i = 1, Model % Mesh % NumberOfNodes
        IF (IntViscoPerm(i)>0) THEN
          IntVisco(IntViscoPerm(i)) = VariableValues(Permutation(i)) 
          IF (Depth(DepthPerm(i))>0.0_dp) IntVisco(IntViscoPerm(i)) = &
               IntVisco(IntViscoPerm(i)) / Depth(DepthPerm(i))
        END IF
      END DO
    ELSE IF (COMP==2) THEN
      DO i = 1, Model % Mesh % NumberOfNodes
        IF (IntDensPerm(i)>0) THEN
          IntDens(IntDensPerm(i)) = VariableValues(Permutation(i)) 
          IF (Depth(DepthPerm(i))>0.0_dp) IntDens(IntDensPerm(i)) = &
               IntDens(IntDensPerm(i)) / Depth(DepthPerm(i))


        END IF
      END DO
    END IF

  END DO !COMP


CONTAINS

  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, Element, n, var)
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), var(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), detJ, grad
    LOGICAL :: Stat
    INTEGER :: t, p,q ,dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    dim = CoordinateSystemDimension()

    IP = GaussPoints( Element )
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

      grad  = SUM( var(1:n) * dBasisdx(1:n,dim) )
      FORCE(1:n) = FORCE(1:n) + grad * IP % s(t) * DetJ  * Basis(1:n)

      DO p=1,n
        DO q=1,n
          STIFF(p,q) = STIFF(p,q) + IP % S(t) * detJ * dBasisdx(q,dim)*dBasisdx(p,dim)
        END DO
      END DO
    END DO

    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC(  STIFF, FORCE, Element, n, var ) 
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), var(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3), &
         DetJ,Normal(3), eta, grad 
    LOGICAL :: Stat
    INTEGER :: t, dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    dim = CoordinateSystemDimension()

    IP = GaussPoints( Element )
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

      grad  = SUM( var(1:n) * Basis(1:n) )

      Normal = NormalVector( Element, Nodes, IP % U(t), IP % V(t), .TRUE.)
      FORCE(1:n) = FORCE(1:n) - grad * IP % s(t) * DetJ * Normal(dim) * Basis(1:n)
    END DO
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
  !------------------------------------------------------------------------------
END SUBROUTINE GetMeanValueSolver
!------------------------------------------------------------------------------




! *****************************************************************************
SUBROUTINE SSASolver( Model,Solver,dt,TransientSimulation )
!DEC$ATTRIBUTES DLLEXPORT :: SSASolver
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Export vertically the SSABasal Velocity (given as a Dirichlet Boundary condition)
!  Compute also the vertical SSAvelocity and the cryostatic pressure
!  If SSA Coupling = Logical True, then coupled the found SSA solution
!   (SSAmask = -1) to the FS solution (SSAmask = 1)
!   and FS solution also where SSAmask = 0 (interface),
!   using w_FS as Dirichlet BC for w_SSA.
!   NOTE: may be discontinuous in coupled vertical velocity since incompressibility only gives dw/dz
!   for w_SSA, so only known up to a constant. Therefore w_SSA cannot be used to update
!   the surface (with FreeSurfSolver), use ThicknessSolver and impose floatation criterium
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
!> Updates vertical velocity from horizontal SSA solution
  
  USE DefUtils
  USE SolverUtils !necessary for SetSinglePoint

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: CurrentElement, Element
  TYPE(Matrix_t),POINTER  :: StiffMatrix
  TYPE(ValueList_t), POINTER :: SolverParams, BodyForce, Material
  TYPE(Variable_t), POINTER :: PointerToVariable, Grad1Sol, Grad2Sol, &
       DepthSol, VeloSol, MaskSol, CoupledSol

  LOGICAL :: AllocationsDone = .FALSE., Found,UnFoundFatal=.TRUE., GotIt, &
       SSACoupling

  INTEGER :: i, j, n, m, t, istat, DIM, p, Indexes(128), COMP
  INTEGER, POINTER :: Permutation(:), VeloPerm(:), &
       DepthPerm(:), GradSurface1Perm(:), GradSurface2Perm(:), &
       NodeIndexes(:), MaskPerm(:), CoupledPerm(:)

  REAL(KIND=dp), POINTER :: ForceVector(:)
  REAL(KIND=dp), POINTER :: VariableValues(:), Depth(:), GradSurface1(:), &
       GradSurface2(:), Velocity(:), PrevVelo(:,:), Mask(:), Coupled(:)
  REAL(KIND=dp) :: Norm, cn, dd

  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:), &
       NodalGravity(:), NodalDensity(:), &
       NodalDepth(:), NodalSurfGrad1(:), NodalSurfGrad2(:), &
       NodalU(:), NodalV(:)


  !added for glueing, remove all timing if complains, not interesting now
  !REAL(KIND=dp) :: gluetime, CPUTime
  INTEGER :: NSDOFs
  INTEGER, POINTER :: FlowPerm(:)
  TYPE(Variable_t), POINTER :: FlowSol
  REAL(KIND=dp), POINTER :: FlowSolution(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, FreeSurfGradName, MaskName, FSName


  SAVE STIFF, LOAD, FORCE, AllocationsDone, DIM, SolverName
  SAVE NodalGravity, NodalDensity, &
       NodalDepth, NodalSurfGrad1, NodalSurfGrad2, &
       NodalU, NodalV
  !------------------------------------------------------------------------------
  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values
  WRITE(SolverName, '(A)') 'SSASolver'
  !------------------------------------------------------------------------------
  !    Get variables needed for solution
  !------------------------------------------------------------------------------
  DIM = CoordinateSystemDimension()

  SolverParams => GetSolverParams()

  SSACoupling = GetLogical(SolverParams, 'SSA Coupling', GotIt)
  IF (.NOT.GotIt .OR. .NOT.SSACoupling) THEN
    SSACoupling = .FALSE.
    CALL INFO(SolverName,'No FS - SSA coupling',Level=4)
  ELSE
    CALL INFO(SolverName,'FS - SSA coupling active',Level=1)
  END IF


  VeloSol => VariableGet( Solver % Mesh % Variables, 'SSAFlow',UnFoundFatal=UnFoundFatal)
  Velocity => VeloSol % Values
  VeloPerm => VeloSol % Perm
  PrevVelo => veloSol % PrevValues

  IF (SSACoupling) THEN
    CoupledSol => VariableGet( Solver % Mesh % Variables, 'CoupledFlow',UnFoundFatal=UnFoundFatal)
    Coupled => CoupledSol % Values
    CoupledPerm => CoupledSol % Perm
  END IF

  DepthSol => VariableGet( Solver % Mesh % Variables, 'Depth',UnFoundFatal=UnFoundFatal)
  Depth => DepthSol % Values
  DepthPerm => DepthSol % Perm


  FreeSurfGradName = GetString(SolverParams, 'Free Surface Grad 1 name', GotIt)
  IF (GotIt) THEN
    CALL INFO(SolverName, 'Free Surface Grad 1 name found', level=8)
    Grad1Sol => VariableGet(Model % Mesh % Variables, FreeSurfGradName,UnFoundFatal=UnFoundFatal)
  ELSE
    !original implementation
    Grad1Sol => VariableGet( Solver % Mesh % Variables, 'FreeSurfGrad1',UnFoundFatal=UnFoundFatal)
  END IF
  GradSurface1 => Grad1Sol % Values
  GradSurface1Perm => Grad1Sol % Perm

  IF (dim > 2) THEN
    FreeSurfGradName = GetString(SolverParams, 'Free Surface Grad 2 name', GotIt)
    IF (GotIt) THEN
      CALL INFO(SolverName, 'Free Surface Grad 2 name found', level=8)
      Grad2Sol => VariableGet(Model % Mesh % Variables, FreeSurfGradName,UnFoundFatal=UnFoundFatal)
    ELSE
      !original implementation
      Grad2Sol => VariableGet( Solver % Mesh % Variables, 'FreeSurfGrad2',UnFoundFatal=UnFoundFatal)
    END IF
    GradSurface2 => Grad2Sol % Values
    GradSurface2Perm => Grad2Sol % Perm
  END IF

  !>>>> start coupling added by eef
  ! mask value geq 0: use FS, mask value -1; compute w_ssa, p_ssa and glue to fs
  IF (SSACoupling) THEN
    MaskName = GetString(SolverParams, 'Mask name', GotIt)
    IF (GotIt) THEN
      CALL info(SolverName, 'Mask found', level=8)
      MaskSol => VariableGet(Model % Mesh % Variables, MaskName,UnFoundFatal=UnFoundFatal)
      Mask => MaskSol % Values
      MaskPerm => MaskSol % Perm 
    ELSE
      print *, 'Warning: mask not found'
    END IF

    FSName = GetString(SolverParams, 'Full Stokes name', GotIt)
    IF (GotIt) THEN
      CALL info(SolverName, 'FS name found', level=8)
      FlowSol => VariableGet(Model % Mesh % Variables, FSName,UnFoundFatal=.FALSE.)
    ELSE
      CALL info(SolverName, 'Warning: FS name not found',level=3)
      FlowSol => VariableGet(Model % Mesh % Variables, 'Flow Solution',UnFoundFatal=UnFoundFatal)
      !take Flow Solution as default!
    END IF
    NSDOFs = FlowSol % DOFs
    FlowPerm => FlowSol % Perm
    FlowSolution  => FlowSol % Values
  END IF

  !--------------------------------------------------------------
  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed  ) THEN
    N = Solver % Mesh % MaxElementNodes ! just big enough for elemental arrays
    M = Model % Mesh % NumberOfNodes
    IF (AllocationsDone) DEALLOCATE(FORCE, LOAD, STIFF, NodalGravity, &
         NodalDensity, NodalDepth, &
         NodalSurfGrad1, NodalSurfGrad2, NodalU, NodalV)

    ALLOCATE( FORCE(N), LOAD(N), STIFF(N,N), &
         NodalGravity(N), NodalDensity(N), &
         NodalDepth(N), NodalSurfGrad1(N), NodalSurfGrad2(N), &
         NodalU(N), NodalV(N), STAT=istat )
    IF ( istat /= 0 ) THEN
      CALL Fatal( SolverName, 'Memory allocation error.' )
    END IF


    AllocationsDone = .TRUE.
    CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
  END IF



  StiffMatrix => Solver % Matrix
  ForceVector => StiffMatrix % RHS

  ! Loop over the velocity components and pressure
  ! If DIM = 2 u, w, p
  ! If DIM = 3 u, v, w, p
  !-----------------------------------------------
  DO  COMP = 1, DIM +1 !temporary hack since problem with pressureDIM+1

    ! No non-linear iteration, no time dependency
    VariableValues = 0.0d0
    Norm = Solver % Variable % Norm

    !Initialize the system and do the assembly:
    !------------------------------------------
    CALL DefaultInitialize()
    ! bulk assembly
    DO t=1,Solver % NumberOfActiveElements
      Element => GetActiveElement(t)
      !below for testing purpose, having problems with setting passive!
      IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
      n = GetElementNOFNodes()

      NodeIndexes => Element % NodeIndexes

      ! Read the gravity in the Body Force Section
      BodyForce => GetBodyForce()
      NodalGravity = 0.0_dp
      IF ( ASSOCIATED( BodyForce ) ) THEN
        IF (DIM==2) THEN
          NodalGravity(1:n) = ListGetReal( &
               BodyForce, 'Flow BodyForce 2', n, NodeIndexes, Found)
        ELSE
          NodalGravity(1:n) = ListGetReal( &
               BodyForce, 'Flow BodyForce 3', n, NodeIndexes, Found)
        END IF
      END IF

      ! Read the Viscosity eta, density, and exponent m in Material Section
      ! Same definition as NS Solver in Elmer - n=1/m , A = 1/ (2 eta^n)
      Material => GetMaterial()

      NodalDensity = 0.0D0
      NodalDensity(1:n) = ListGetReal( &
           Material, 'Density', n, NodeIndexes, Found )

      ! Get the Nodal value of Depth, FreeSurfGrad1 and FreeSurfGrad2
      NodalDepth(1:n) = Depth(DepthPerm(NodeIndexes(1:n)))
      NodalSurfGrad1(1:n) = GradSurface1(GradSurface1Perm(NodeIndexes(1:n)))
      NodalSurfGrad2 = 0.0D0
      IF (DIM==3) NodalSurfGrad2(1:n) = GradSurface2(GradSurface2Perm(NodeIndexes(1:n)))

      IF (COMP==1) THEN     ! u
        CALL LocalMatrixUV (  STIFF, FORCE, Element, n )

      ELSE IF (COMP==DIM) THEN  ! w
        NodalU(1:n) = Velocity((DIM+1)*(VeloPerm(NodeIndexes(1:n))-1)+1)
        NodalV = 0.0D0
        IF (DIM==3) NodalV(1:n) = Velocity((DIM+1)*(VeloPerm(NodeIndexes(1:n))-1)+2)
        CALL LocalMatrixW (  STIFF, FORCE, Element, n, NodalU, NodalV )

      ELSE IF (COMP==DIM+1) THEN ! p
        CALL LocalMatrixP (  STIFF, FORCE, Element, n )

      ELSE               ! v if dim=3
        CALL LocalMatrixUV (  STIFF, FORCE, Element, n )

      END IF

      CALL DefaultUpdateEquations( STIFF, FORCE )
    END DO

    ! Neumann conditions only for w and p
    IF (COMP .GE. DIM) THEN
      DO t=1,Solver % Mesh % NUmberOfBoundaryElements
        Element => GetBoundaryElement(t)
        IF ( GetElementFamily() == 1 ) CYCLE
        NodeIndexes => Element % NodeIndexes
        IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
        n = GetElementNOFNodes()
        STIFF = 0.0D00
        FORCE = 0.0D00

        IF (COMP==DIM) THEN
          ! only for the surface nodes
          dd = SUM(ABS(Depth(Depthperm(NodeIndexes(1:n)))))
          IF (dd < 1.0e-6) THEN
            NodalU(1:n) = Velocity((DIM+1)*(VeloPerm(NodeIndexes(1:n))-1)+1)
            NodalV = 0.0D0
            IF (DIM==3) NodalV(1:n) = Velocity((DIM+1)*(VeloPerm(NodeIndexes(1:n))-1)+2)
            CALL LocalMatrixBCW (  STIFF, FORCE, Element, n, NodalU, NodalV )
          END IF
        ELSE IF (COMP==DIM+1) THEN
          CALL LocalMatrixBCP(  STIFF, FORCE, Element, n, NodalDensity, &
               NodalGravity )
        END IF
        CALL DefaultUpdateEquations( STIFF, FORCE )
      END DO
    END IF

    CALL DefaultFinishAssembly()

    ! Dirichlet
    CALL SetDirichletBoundaries( Model, StiffMatrix, ForceVector, &
         ComponentName('SSAFlow',COMP), 1,1, Permutation )

    IF (SSACoupling .AND. (COMP==DIM)) THEN
      DO i = 1, Model % Mesh % NumberOfNodes
        IF (Mask(MaskPerm(i))==0) THEN 
          CALL SetDirichletPoint(StiffMatrix,ForceVector,1,1,Permutation,i,FlowSolution(NSDOFs*(FlowPerm(i)-1)+DIM))
          !setting w_FS as dirichlet for vertical component of SSAFlow
        END IF
      END DO
    END IF

    !Solve the system
    Norm = DefaultSolve()
    ! Save the solution on the right variable
    DO i = 1, Model % Mesh % NumberOfNodes 
      IF (VeloPerm(i)>0) THEN
        Velocity ((DIM+1)*(VeloPerm(i)-1) + COMP) = VariableValues(Permutation(i))
      END IF

    END DO

  END DO ! Loop p

  IF (SSACoupling) THEN
    !Glue solution together with fs-solution
    WRITE( Message, * ) 'Glue FS-solution and SSA-solution together'
    CALL Info( 'SSACoupler',Message, Level=4 )
    j=0

    DO i = 1, SIZE(Coupled)/NSDOFs !goes through rows in flowsolution
      IF (Mask(MaskPerm(i)) .EQ. -1) THEN !SSA
        DO j=1,NSDOFs
          Coupled(NSDOFs*(CoupledPerm(i)-1)+j)= Velocity(NSDOFs*(VeloPerm(i)-1)+j)
        END DO
      ELSE ! Full Stokes
        DO j=1,NSDOFs
          Coupled(NSDOFs*(CoupledPerm(i)-1)+j) = FlowSolution(NSDOFs*(FlowPerm(i)-1)+j)
        END DO
      END IF
    END DO

  END IF
CONTAINS

  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixUV(  STIFF, FORCE, Element, n )
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), detJ
    LOGICAL :: Stat
    INTEGER :: t, p, q , dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    dim = CoordinateSystemDimension()


    IP = GaussPoints( Element )
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

      DO p=1,n
        DO q=1,n
          STIFF(p,q) = STIFF(p,q) + IP % S(t) * detJ * dBasisdx(q,dim)*dBasisdx(p,dim)
        END DO
      END DO

    END DO
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixUV
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixW(  STIFF, FORCE, Element, n, VeloU, VeloV)
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), VeloU(:), VeloV(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), detJ, &
         dU2dxz, dV2dyz
    LOGICAL :: Stat
    INTEGER :: t, p,q , DIM
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    DIM = CoordinateSystemDimension()

    IP = GaussPoints( Element )
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .TRUE. )

      DO p=1,n
        DO q=1,n
          STIFF(p,q) = STIFF(p,q) + IP % S(t) * detJ * dBasisdx(q,dim)*dBasisdx(p,dim)
        END DO
      END DO

      dU2dxz = SUM(VeloU(1:n)*ddBasisddx(1:n,1,dim))
      dV2dyz = 0.0d0
      IF (DIM==3) dV2dyz = SUM(VeloV(1:n)*ddBasisddx(1:n,2,3))


      FORCE(1:n) = FORCE(1:n) + (dU2dxz + dV2dyz) * IP % s(t) * detJ * Basis(1:n)

    END DO

    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixW

  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixP(  STIFF, FORCE, Element, n)
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), detJ
    LOGICAL :: Stat
    INTEGER :: t, p,q ,dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    dim = CoordinateSystemDimension()

    IP = GaussPoints( Element )
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

      DO p=1,n
        DO q=1,n
          STIFF(p,q) = STIFF(p,q) + IP % S(t) * detJ * dBasisdx(q,dim)*dBasisdx(p,dim)
        END DO
      END DO
    END DO

    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixP
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBCW(  STIFF, FORCE, Element, n, VeloU, VeloV )
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), veloU(:), veloV(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3), &
         DetJ, Normal(3), grad, dUdx, dVdy
    LOGICAL :: Stat
    INTEGER :: t, DIM
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    DIM = CoordinateSystemDimension()

    IP = GaussPoints( Element )
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

      dUdx = SUM( VeloU(1:n) * dBasisdx(1:n,1) )
      dVdy = 0.0e0
      IF (DIM==3) dVdy = SUM( VeloV(1:n) * dBasisdx(1:n,2) )

      grad = - (dUdx + dVdy)

      Normal = NormalVector( Element, Nodes, IP % U(t), IP % V(t), .TRUE.)
      FORCE(1:n) = FORCE(1:n) + grad * IP % s(t) * DetJ * Normal(dim) * Basis(1:n)
    END DO
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBCW
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBCP(  STIFF, FORCE, Element, n, Density, &
       Gravity)
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), density(:), Gravity(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3), &
         DetJ,Normal(3), rho, g, grad
    LOGICAL :: Stat
    INTEGER :: t, dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    dim = CoordinateSystemDimension()

    IP = GaussPoints( Element )
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

      g = ABS(SUM( Gravity(1:n) * Basis(1:n) ))
      rho = SUM( Density(1:n) * Basis(1:n) )

      grad = - rho * g

      Normal = NormalVector( Element, Nodes, IP % U(t), IP % V(t), .TRUE.)
      FORCE(1:n) = FORCE(1:n) + grad * IP % s(t) * DetJ * Normal(dim) * Basis(1:n)
    END DO
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBCP
  !------------------------------------------------------------------------------
END SUBROUTINE SSASolver
!------------------------------------------------------------------------------

