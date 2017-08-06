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
SUBROUTINE AdjointSSA_SSASolver( Model,Solver,dt,TransientSimulation )
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
  TYPE(Variable_t), POINTER :: PointerToVariable, ZsSol, ZbSol, &
                               VeloSol

  LOGICAL :: AllocationsDone = .FALSE., Found, GotIt, CalvingFront 
  LOGICAL :: Newton

  INTEGER :: i, n, m, t, istat, DIM, p, STDOFs
  INTEGER :: NonlinearIter, NewtonIter, iter, other_body_id
          
  INTEGER, POINTER :: Permutation(:), &
       ZsPerm(:), ZbPerm(:), &
       NodeIndexes(:)

  REAL(KIND=dp), POINTER :: ForceVector(:)
  REAL(KIND=dp), POINTER :: VariableValues(:), Zs(:), Zb(:)
                            
  REAL(KIND=dp) :: UNorm, cn, dd, NonlinearTol, NewtonTol, MinSRInv, rhow, sealevel, &
                   PrevUNorm, relativeChange,minv

  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:), &
           NodalGravity(:), NodalViscosity(:), NodalDensity(:), &
           NodalZs(:), NodalZb(:),   &
           NodalU(:), NodalV(:), NodalBeta(:),LocalLinVelo(:)

  INTEGER :: iFriction
  REAL(KIND=dp) :: fm
  CHARACTER(LEN=MAX_NAME_LEN) :: Friction
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
#ifdef USE_ISO_C_BINDINGS
    REAL(KIND=dp) :: at, at0
#else
    REAL(KIND=dp) :: at, at0, CPUTime, RealTime
#endif 
       
  SAVE rhow,sealevel
  SAVE STIFF, LOAD, FORCE, AllocationsDone, DIM, SolverName, ElementNodes
  SAVE NodalGravity, NodalViscosity, NodalDensity, &
           NodalZs, NodalZb,   &
           NodalU, NodalV, NodeIndexes, NodalBeta,LocalLinVelo

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


        ZbSol => VariableGet( Solver % Mesh % Variables, 'Zb' )
        IF (ASSOCIATED(ZbSol)) THEN
           Zb => ZbSol % Values
           ZbPerm => ZbSol % Perm
        ELSE
           CALL FATAL(SolverName,'Could not find variable >Zb<')
        END IF

        ZsSol => VariableGet( Solver % Mesh % Variables, 'Zs' )
        IF (ASSOCIATED(ZsSol)) THEN
           Zs => ZsSol % Values
           ZsPerm => ZsSol % Perm
        ELSE
           CALL FATAL(SolverName,'Could not find variable >Zs<')
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

     sealevel = GetConstReal( Model % Constants, 'Sea Level', Found )
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
                       NodalBeta,LocalLinVelo, ElementNodes % x, &
                       ElementNodes % y, ElementNodes % z )

     ALLOCATE( FORCE(STDOFs*N), LOAD(N), STIFF(STDOFs*N,STDOFs*N), &
          NodalGravity(N), NodalDensity(N), NodalViscosity(N), &
          NodalZb(N), NodalZs(N) ,&
          NodalU(N), NodalV(N), NodalBeta(N),LocalLinVelo(N), &
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


     NodalDensity=0.0_dp
     NodalDensity(1:n) = ListGetReal( Material, 'SSA Mean Density',n,NodeIndexes,Found)
     IF (.NOT.Found) &
           CALL FATAL(SolverName,'Could not find Material prop.  >SSA Mean Density<')

     NodalViscosity=0.0_dp
     NodalViscosity(1:n) = ListGetReal( Material, 'SSA Mean Viscosity',n, NodeIndexes,Found)
     IF (.NOT.Found) &
          CALL FATAL(SolverName,'Could not find Material prop. >SSA Mean Viscosity<')

     !NodalSliding = 0.0_dp
     !NodalSliding(1,1:n) = ListGetReal( &
     !      Material, 'SSA Slip Coefficient 1', n, NodeIndexes(1:n), Found )
     !IF (STDOFs==2) THEN
     !   NodalSliding(2,1:n) = ListGetReal( &
     !        Material, 'SSA Slip Coefficient 2', n, NodeIndexes(1:n), Found )  
     !END IF
     Friction = GetString(Material, 'SSA Friction Law', Found)
     IF (.NOT.Found) &
        CALL FATAL(SolverName,'Could not find Material keyword >SSA Friction Law<')

    SELECT CASE(Friction)
       CASE('linear')
         iFriction = 1
         fm = 1.0_dp
       CASE('weertman')
        iFriction = 2
       CASE DEFAULT
         CALL FATAL(SolverName,'Friction should be linear or Weertman')
   END SELECT


     NodalBeta(1:n) = GetReal( Material, 'SSA Friction Parameter', Found, Element)
     IF (.NOT.Found) &
        CALL FATAL(SolverName,'Could not find Material prop. >SSA Friction Parameter<')
     IF (iFriction > 1) THEN
           fm = ListGetConstReal( Material, 'SSA Friction Exponent', Found )
           IF (.NOT.Found) &
               CALL FATAL(SolverName,'Could not find Material prop. >SSA Friction Exponent<')
           LocalLinVelo = 0.0_dp
           LocalLinVelo(1:n) = GetReal(Material, 'SSA Friction Linear Velocity', Found, Element)
           IF (.NOT.Found) &
                    CALL FATAL(SolverName,'Could not find  Material prop. >SSA Friction Linear Velocity<')
     END IF



     ! Get the Nodal value of Zb and Zs
     NodalZb(1:n) = Zb(ZbPerm(NodeIndexes(1:n)))
     NodalZs(1:n) = Zs(ZsPerm(NodeIndexes(1:n)))

     ! Previous Velocity 
     NodalU(1:n) = VariableValues(STDOFs*(Permutation(NodeIndexes(1:n))-1)+1)
     NodalV = 0.0
     IF (STDOFs.EQ.2) NodalV(1:n) = VariableValues(STDOFs*(Permutation(NodeIndexes(1:n))-1)+2)
      

     CALL LocalMatrixUVSSA (  STIFF, FORCE, Element, n, ElementNodes, NodalGravity, &
        NodalDensity, NodalViscosity, NodalZb, NodalZs, &
        NodalU, NodalV, NodalBeta,iFriction,fm,LocalLinVelo, cn, MinSRInv , STDOFs, Newton)

     CALL DefaultUpdateEquations( STIFF, FORCE )

  END DO
  CALL DefaultFinishBulkAssembly()
  
!  
! Neumann condition
!
  DO t=1,GetNOFBoundaryElements()
     BoundaryElement => GetBoundaryElement(t)
     IF ( .NOT. ActiveBoundaryElement() ) CYCLE
     IF ( GetElementFamily() == 1 ) CYCLE

     NodeIndexes => BoundaryElement % NodeIndexes
     IF (ParEnv % myPe .NE. BoundaryElement % partIndex) CYCLE

     n = GetElementNOFNodes()
     FORCE = 0.0e0
     STIFF = 0.0e0

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
        NodalDensity(1:n) = ListGetReal( Material, 'SSA Mean Density',n, NodeIndexes,Found)
        IF (.NOT.Found) &
           CALL FATAL(SolverName,'Could not find Material prop.  >SSA Mean Density<')

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
               NodalDensity, NodalGravity, NodalZb, NodalZs, rhow, sealevel )
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

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixUVSSA(  STIFF, FORCE, Element, n, Nodes, gravity, &
           Density, Viscosity, LocalZb, LocalZs, LocalU, &
           LocalV, LocalBeta,iFriction,fm,LocalLinVelo, cm, MinSRInv, STDOFs , Newton )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), gravity(:), Density(:), &
                     Viscosity(:), LocalZb(:), LocalZs(:), &
                     LocalU(:), LocalV(:) , LocalBeta(:), LocalLinVelo(:)
    INTEGER :: n, cp , STDOFs
    INTEGER :: iFriction
    REAL(KIND=dp) :: cm,fm
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: Newton
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), detJ 
    REAL(KIND=dp) :: g, rho, eta, h, dhdx, dhdy , muder
    REAL(KIND=dp) :: gradS(2),Slip,A(2,2), StrainA(2,2),StrainB(2,2), Exx, Eyy, Exy, Ezz, Ee, MinSRInv 
    REAL(KIND=dp) :: beta,Slip2,Velo(2),LinVelo,ub
    REAL(KIND=dp) :: Jac(2*n,2*n),SOL(2*n)
    LOGICAL :: Stat, NewtonLin,fNewtonLin
    INTEGER :: i, j, t, p, q , dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
!------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()

    STIFF = 0.0d0
    FORCE = 0.0d0
    Jac=0.0d0

! Use Newton Linearisation
    NewtonLin=(Newton.AND.(cm.NE.1.0_dp))
    fNewtonLin = (Newton.AND.(fm.NE.1.0_dp))


    IP = GaussPoints( Element )
    DO t=1,IP % n
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
        IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

! Needed Integration Point value

       g = ABS(SUM( Gravity(1:n) * Basis(1:n) ))
       rho = SUM( Density(1:n) * Basis(1:n) )
       eta = SUM( Viscosity(1:n) * Basis(1:n) )
       gradS = 0._dp
       gradS(1) = SUM( LocalZs(1:n) * dBasisdx(1:n,1) )
       if (STDOFs == 2) gradS(2) = SUM( LocalZs(1:n) * dBasisdx(1:n,2) )
       h = SUM( (LocalZs(1:n)-LocalZb(1:n)) * Basis(1:n) )
       
       beta = SUM( LocalBeta(1:n) * Basis(1:n) )
      ! slip = 0.0_dp
      ! DO i=1,STDOFs
      !    slip(i) = SUM( LocalSliding(i,1:n) * Basis(1:n) )
      ! END DO

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
   IF (iFriction==1) THEN
      Slip = beta
      fNewtonLin = .FALSE.
   ELSE IF (iFriction==2) THEN
      Slip = beta * ub**(fm-1.0_dp) 
      Slip2 = Slip2*Slip*(fm-1.0_dp)/(ub*ub)
   END IF


!------------------------------------------------------------------------------
! In the non-linear case, effective viscosity       
       IF (cm.NE.1.0_dp) THEN
           Exx = SUM(LocalU(1:n)*dBasisdx(1:n,1))
           Eyy = 0.0
           Exy = 0.0
           IF (STDOFs.EQ.2) THEN
              Eyy = SUM(LocalV(1:n)*dBasisdx(1:n,2))
              Ezz = -Exx - Eyy
              Exy = SUM(LocalU(1:n)*dBasisdx(1:n,2))
              Exy = 0.5*(Exy + SUM(LocalV(1:n)*dBasisdx(1:n,1)))
              Ee = 0.5*(Exx**2.0 + Eyy**2.0 + Ezz**2.0) + Exy**2.0
              !Ee = SQRT(Ee)
           ELSE
              !Ee = ABS(Exx)
              Ee = Exx * Exx
           END IF
           muder = eta * 0.5 * (2**cm) * ((cm-1.0)/2.0) *  Ee**((cm-1.0)/2.0 - 1.0)
           IF (sqrt(Ee) < MinSRInv) then
                Ee = MinSRInv*MinSRInv
                muder = 0.0_dp
           Endif
           eta = eta * 0.5 * (2**cm) * Ee**((cm-1.0)/2.0)
       END IF 

       StrainA=0.0_dp
       StrainB=0.0_dp
       If (NewtonLin) then
          StrainA(1,1)=SUM(2.0*dBasisdx(1:n,1)*LocalU(1:n))

          IF (STDOFs.EQ.2) THEN
             StrainB(1,1)=SUM(0.5*dBasisdx(1:n,2)*LocalU(1:n))

             StrainA(1,2)=SUM(dBasisdx(1:n,2)*LocalV(1:n))
             StrainB(1,2)=SUM(0.5*dBasisdx(1:n,1)*LocalV(1:n))

             StrainA(2,1)=SUM(dBasisdx(1:n,1)*LocalU(1:n))
             StrainB(2,1)=SUM(0.5*dBasisdx(1:n,2)*LocalU(1:n))

             StrainA(2,2)=SUM(2.0*dBasisdx(1:n,2)*LocalV(1:n))
             StrainB(2,2)=SUM(0.5*dBasisdx(1:n,1)*LocalV(1:n))

          End if
       Endif

       A = 0.0_dp
       DO p=1,n
         DO q=1,n
         A(1,1) = 2.0*dBasisdx(q,1)*dBasisdx(p,1)  
           IF (STDOFs.EQ.2) THEN
           A(1,1) = A(1,1) + 0.5*dBasisdx(q,2)*dBasisdx(p,2)
           A(1,2) = dBasisdx(q,2)*dBasisdx(p,1) + &
                             0.5*dBasisdx(q,1)*dBasisdx(p,2)
           A(2,1) = dBasisdx(q,1)*dBasisdx(p,2) + &
                             0.5*dBasisdx(q,2)*dBasisdx(p,1)
           A(2,2) = 2.0*dBasisdx(q,2)*dBasisdx(p,2) +&
                             0.5*dBasisdx(q,1)*dBasisdx(p,1)  
          END IF
           A = 2.0 * h * eta * A
           DO i=1,STDOFs
             STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+i) = STIFF((STDOFs)*(p-1)+i,(STDOFs)*(q-1)+i) +&
                  slip * Basis(q) * Basis(p) * IP % S(t) * detJ
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

           If (NewtonLin) then
            ! Maybe a more elegant formulation to get the Jacobian??.......
            IF (STDOFs.EQ.1) THEN
                  Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+1) = Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+1) +&
                        IP % S(t) * detJ * 2.0 * h * StrainA(1,1)*dBasisdx(p,1) * &
                         muder * 2.0 * Exx*dBasisdx(q,1) 

             ELSE IF (STDOFs.EQ.2) THEN
                  Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+1) = Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+1) +&
             IP % S(t) * detJ * 2.0 * h * ((StrainA(1,1)+StrainA(1,2))*dBasisdx(p,1)+(StrainB(1,1)+StrainB(1,2))*dBasisdx(p,2)) * &
             muder *((2.0*Exx+Eyy)*dBasisdx(q,1)+Exy*dBasisdx(q,2)) 

                  Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+2) = Jac((STDOFs)*(p-1)+1,(STDOFs)*(q-1)+2) +&
             IP % S(t) * detJ * 2.0 * h * ((StrainA(1,1)+StrainA(1,2))*dBasisdx(p,1)+(StrainB(1,1)+StrainB(1,2))*dBasisdx(p,2)) * &
             muder *((2.0*Eyy+Exx)*dBasisdx(q,2)+Exy*dBasisdx(q,1)) 

                  Jac((STDOFs)*(p-1)+2,(STDOFs)*(q-1)+1) = Jac((STDOFs)*(p-1)+2,(STDOFs)*(q-1)+1) +&
             IP % S(t) * detJ * 2.0 * h * ((StrainA(2,1)+StrainA(2,2))*dBasisdx(p,2)+(StrainB(2,1)+StrainB(2,2))*dBasisdx(p,1)) * &
             muder *((2.0*Exx+Eyy)*dBasisdx(q,1)+Exy*dBasisdx(q,2)) 

                  Jac((STDOFs)*(p-1)+2,(STDOFs)*(q-1)+2) = Jac((STDOFs)*(p-1)+2,(STDOFs)*(q-1)+2) +&
             IP % S(t) * detJ * 2.0 * h * ((StrainA(2,1)+StrainA(2,2))*dBasisdx(p,2)+(StrainB(2,1)+StrainB(2,2))*dBasisdx(p,1)) * &
             muder *((2.0*Eyy+Exx)*dBasisdx(q,2)+Exy*dBasisdx(q,1)) 
             End if
           Endif

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

    If (NewtonLin) then
         SOL(1:STDOFs*n:STDOFs)=LocalU(1:n)
         If (STDOFs.EQ.2) SOL(2:STDOFs*n:STDOFs)=LocalV(1:n)

         STIFF(1:STDOFs*n,1:STDOFs*n) = STIFF(1:STDOFs*n,1:STDOFs*n) + &
                                        Jac(1:STDOFs*n,1:STDOFs*n)
         FORCE(1:STDOFs*n) = FORCE(1:STDOFs*n) + &
                             MATMUL(Jac(1:STDOFs*n,1:STDOFs*n),SOL(1:STDOFs*n))
    Endif
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixUVSSA
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBCSSA(  STIFF, FORCE, Element, n, ENodes, Density, & 
                      Gravity, LocalZb, LocalZs, rhow, sealevel)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) ::  ENodes
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:),  density(:), Gravity(:), LocalZb(:),&
                         LocalZs(:),rhow, sealevel
    INTEGER :: n
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3), &
                      DetJ,Normal(3), rhoi, g, alpha, h, h_im,norm
    LOGICAL :: Stat
    INTEGER :: t, i
    TYPE(GaussIntegrationPoints_t) :: IP

!------------------------------------------------------------------------------
    STIFF = 0.0d0
    FORCE = 0.0d0

! The front force is a concentrated nodal force in 1D-SSA and
! a force distributed along a line in 2D-SSA    

! 1D-SSA Case : concentrated force at each nodes
    IF (STDOFs==1) THEN  !1D SSA but should be 2D problem (does elmer work in 1D?)
      DO i = 1, n
         g = ABS( Gravity(i) )
         rhoi = Density(i)
         h = LocalZs(i)-LocalZb(i) 
         h_im=max(0._dp,sealevel-LocalZb(i))
         alpha=0.5 * g * (rhoi * h**2.0 - rhow * h_im**2.0)
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
             alpha=0.5 * g * (rhoi * h**2.0 - rhow * h_im**2.0)

! Normal in the (x,y) plane
             Normal = NormalVector( Element, ENodes, IP % U(t), IP % V(t), .TRUE.)
             norm=SQRT(normal(1)**2.0+normal(2)**2.0)
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
END SUBROUTINE AdjointSSA_SSASolver
!------------------------------------------------------------------------------

