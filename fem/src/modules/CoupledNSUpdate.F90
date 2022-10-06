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
! *  Authors: Juha Ruokolainen, Mika  Malinen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: Dec 10th, 2003
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> This subroutine performs a single GCR step to solve the fully coupled incompressible
!> Navier-Stokes system. The search direction is obtained by solving decoupled
!> equations which arise from the consistent splitting algorithm. This solver
!> assumes that such equations are already solved. The search direction 
!> can thus be constructed without performing linear solves in this module. 
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE OptimalSolutionUpdate( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE SolverUtils
  USE ElementUtils

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
  TYPE(Matrix_t), POINTER :: AMatrix, MMatrix, PMatrix

  LOGICAL :: AllocationsDone = .FALSE., Newton = .FALSE., Hybrid = .FALSE., Found, Convect, &
       GotIt, ConvectionStabilization, BlockPreconditioning = .FALSE., &
       NormalTractionBoundary, GradDivStabilization, Parallel, &
       ConstantSystem, ConstantBulkMatrixInUse

  TYPE(Element_t),POINTER :: Element, Parent
  TYPE(Solver_t),POINTER :: PressureSolver, ProjectionSolver, VelocitySolver
  INTEGER, POINTER :: NodeIndexes(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: OuterIterationMethod, NonlinearIterationMethod, eq
  INTEGER :: i, j, k, n, nb, nd, t, istat, dim, m, p, q, &
       MaxIterations, NumberOfNormalTractionNodes, Active
  REAL(KIND=dp) :: Norm = 0, PrevNorm, RelC, Tolerance, ToleranceRatio, &
       atime, stime, at0, NonLinError

  TYPE(ValueList_t), POINTER :: BodyForce, Material, BC
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:,:), Mass(:,:), &
       FORCE(:), rho(:), mu(:), Velocity(:,:), &
       ALocal(:,:), MLocal(:,:), PLocal(:,:), PrevSol(:), NonLinRes(:), &
       Snew(:), R(:), S(:,:), V(:,:), Residual(:)

  INTEGER, ALLOCATABLE :: TractionBCIndeces(:)

  LOGICAL :: NewTimeStep, OutflowBC, PicardIteration
  INTEGER :: CurrentDoneTime = 0, NonlinearIter, iter, nlen, &
       MaxPicardIterations, Round = 0
  REAL(KIND=dp) :: NonlinearTol, RelaxationFactor, GradDivParam, PresErr, EK, StNumber, &
       NewtonThreshold

  INTEGER, ALLOCATABLE :: Indexes(:)
  TYPE(Variable_t), POINTER :: FlowSol

  REAL(KIND=dp), POINTER CONTIG :: Mx(:),Mb(:),Mr(:)
  TYPE(ValueList_t), POINTER :: SolverParams

  SAVE STIFF, LOAD, FORCE, rho, mu, Velocity, AllocationsDone, &
       ALocal, Mass, TractionBCIndeces, NumberOfNormalTractionNodes, &
       CurrentDoneTime, Indexes, MLocal, Norm, &
       PLocal, PrevSol, NonLinRes, Round, Snew, R, S, V
  !------------------------------------------------------------------------------
  
  NewTimeStep = .FALSE.  
  IF (CurrentDoneTime /= Solver % DoneTime) THEN
     CurrentDoneTime = CurrentDoneTime + 1
     NewTimeStep = .TRUE.
     Round = 0
  END IF
  Round = Round + 1 

  SolverParams => GetSolverParams()
  ConstantSystem = GetLogical( SolverParams, 'Constant System', Found )  
  ConstantBulkMatrixInUse = ConstantSystem .AND. &
       ASSOCIATED(Solver % Matrix % BulkValues) .AND. ( .NOT. NewTimeStep ) .AND. &
       ASSOCIATED(Solver % Matrix % BulkRHS)


  dim = CoordinateSystemDimension()
  ! Check whether convection stabilization is used...
  ConvectionStabilization = ListGetLogical( Solver % Values, 'Stabilize', GotIt )
  IF (.NOT. GotIt) ConvectionStabilization = .FALSE.
  ConvectionStabilization = .FALSE.

  ! The option for using grad-div stabilization is not available currently... 
  GradDivStabilization = ListGetLogical( Solver % Values, 'Grad-Div Stabilization', GotIt )
  IF ( GotIt .AND. GradDivStabilization) THEN
     GradDivParam =  ListGetConstReal( Solver % Values, &
          'Grad-Div Stabilization Parameter', GotIt )
     IF ( .NOT. GotIt ) GradDivParam = 0.0d0
  ELSE
     GradDivStabilization = .FALSE.
     GradDivParam = 0.0d0
  END IF

  NonlinearIterationMethod = ListGetString(Solver % Values, 'Nonlinear Iteration Method', GotIt)
  NonlinearIterationMethod = 'picard'   ! This is the current default...

  IF ( NonlinearIterationMethod == 'hybrid') THEN
     Newton = .TRUE.
     Hybrid = .TRUE.
     PicardIteration = .FALSE.
  ELSE
     IF ( NonlinearIterationMethod /= 'picard') THEN
        PicardIteration = .FALSE.
        Newton = .TRUE.
        Hybrid = .FALSE.        
     ELSE
        PicardIteration = .TRUE.
        Newton = .FALSE.
        Hybrid = .FALSE.
     END IF
  END IF

  !print *, NonlinearIterationMethod
  !print *, 'Newton = ', Newton
  !print *, 'Hybrid = ', Hybrid
  !print *, 'Picard = ', PicardIteration   


  !--------------------------------------------------------------------------
  ! Check whether the block preconditioning is to be used.
  !---------------------------------------------------------------------------
  BlockPreconditioning = ListGetLogical( Solver % Values, 'Block Preconditioning', GotIt )
  IF ( .NOT. GotIt) BlockPreconditioning = .TRUE.
  ! IF (BlockPreconditioning) PRINT *, 'OuterIteration residuals for time level', CurrentDoneTime

  MaxIterations = ListGetInteger( Model % Simulation, &
       'Steady State Max Iterations' )

  !----------------------------------------------------------------------------
  !Allocate some permanent storage, this is done first time only:
  !----------------------------------------------------------------------------
  IF ( .NOT. AllocationsDone ) THEN

     p = Solver % Mesh % MaxElementDOFs       ! The size of the one-component stiffness matrix
     n = (dim+1)*p                            ! The size of the complete stiffness matrix
     m = dim * Solver % Mesh % MaxElementDOFs ! The size of the (1,1)-block (the velocity part)
     !print *, Solver % Matrix % NumberOfRows

     ALLOCATE( &
          FORCE(n), &
          STIFF(n,n), &
          Mass(n,n),&
          LOAD(p,4), &
          rho(p), &
          mu(p), &
          Velocity(dim+1,p), &
          ALocal(m,m), & 
          Indexes(p), &
          MLocal(p,p), &
          PLocal(p,p), &
          PrevSol(Solver % Matrix % NumberOfRows), &
          NonLinRes(Solver % Matrix % NumberOfRows), & 
          Snew( Solver % Matrix % NumberOfRows ), &        ! Consistent splitting update
          R( Solver % Matrix % NumberOfRows ), &         ! Residual
          S( Solver % Matrix % NumberOfRows, MaxIterations), &      ! Search directions
          V( Solver % Matrix % NumberOfRows, MaxIterations), &      ! The range of coefficient matrix 
          STAT=istat )
     
     IF ( istat /= 0 ) THEN
        CALL Fatal( 'NavierStokesSolver', 'Memory allocation error.' )
     END IF

     AllocationsDone = .TRUE.
  END IF


  Convect = GetLogical( GetSolverParams(), 'Convective', Found )
  IF ( .NOT. Found ) Convect = .TRUE.
  IF ( .NOT. Convect ) ConvectionStabilization = .FALSE.
  
  atime = CPUTime()
  at0 = RealTime()

!  NonlinearIter = ListGetInteger( Solver % Values, &
!       'Nonlinear System Max Iterations', minv=0 )
!  NonlinearTol = ListGetConstReal( Solver % Values, &
!       'Nonlinear System Convergence Tolerance',minv=0.0d0 )
!  if ( .not. Convect ) NonlinearIter = 0

  ! Initialize the system matrices and vectors...


  IF ( ConstantBulkMatrixInUse ) THEN

    Solver % Matrix % Values = Solver % Matrix % BulkValues 
    Solver % Matrix % RHS  = Solver % Matrix % BulkRHS

  ELSE

     CALL DefaultInitialize()

     DO t=1,Solver % NumberOfActiveElements
        IF ( RealTime() - at0 > 1.0 ) THEN
           WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
                (Solver % NumberOfActiveElements - t) / &
                (Solver % NumberOfActiveElements)), ' % done'
           CALL Info( 'NavierStokesSolver', Message, Level=5 )
           at0 = RealTime()
        END IF
        Element => GetActiveElement(t)
        n  = GetElementNOFNodes()

        nd = GetElementDOFs( Indexes )
        nb = GetElementNOFBDOFs()

        !-----------------------------------------------
        ! Body forces:
        !-----------------------------------------------
        BodyForce => GetBodyForce()
        LOAD = 0.0d0
        IF ( ASSOCIATED(BodyForce) ) THEN
           Load(1:n,1) = GetReal( BodyForce, 'Force 1', Found )
           Load(1:n,2) = GetReal( BodyForce, 'Force 2', Found )
           Load(1:n,3) = GetReal( BodyForce, 'Force 3', Found )
        END IF

        !-----------------------------------------------
        ! Material parameters:
        !-----------------------------------------------
        Material => GetMaterial()
        !rho(1:n) = GetReal( Material, 'Strouhal' )
        !StNumber = rho(1)
        rho(1:n) = GetReal( Material, 'Density' )
        mu(1:n)  = GetReal( Material, 'Viscosity' )
        !-------------------------------------------------------
        ! Get previous elementwise velocity iterate:
        !-------------------------------------------------------
        Velocity = 0.0d0
        IF ( .NOT. Convect ) THEN
           velocity = 0.0d0
        ELSE
           IF (TransientSimulation) THEN
              DO i=1,dim
                 Velocity(i,1:nd) = Solver % Variable % PrevValues( &
                      Solver % Variable % DOFs*(Solver % Variable % &
                      Perm(Indexes(1:nd))-1)+i,1)
              END DO
           ELSE
              CALL Fatal( 'NavierStokesSolver', 'Steady state problems cannot be handled yet' )
           END IF
        END IF

        !------------------------------------------------------
        ! Get element local matrix and rhs vector:
        !-------------------------------------------------------
        CALL LocalMatrix(  STIFF, Mass, FORCE, LOAD, rho, mu, &
             Velocity, Element, n, nd, dim, ConvectionStabilization, &
             Convect, GradDivParam)

        IF (TransientSimulation) THEN
           CALL Default1stOrderTime( Mass, STIFF, FORCE )
        END IF
        !-----------------------------------------------------------------
        ! Update global matrix and rhs vector from local matrix & vector:
        !----------------------------------------------------------------
        CALL DefaultUpdateEquations( STIFF, FORCE )

     END DO

     IF(.NOT. ConstantBulkMatrixInUse ) THEN
       CALL DefaultFinishBulkAssembly()
     END IF


     ! Terms on outflow boundary
     !--------------------------------------------------------------
     IF (.FALSE.) THEN    ! These are needed only for the rotational form...
        !if (Convect) then
        DO t=1, Solver % Mesh % NumberOfBoundaryElements
           Element => GetBoundaryElement(t)
           IF ( .NOT. ActiveBoundaryElement() ) CYCLE

           n  = GetElementNOFNodes()
           nd = GetElementDOFs( Indexes )

           IF ( GetElementFamily() == 1 ) CYCLE

           BC => GetBC()
           IF ( ASSOCIATED( BC ) ) THEN
              OutflowBC= ListGetLogical( BC, 'Outflow Boundary', Found )

              IF (OutflowBC) THEN

                 IF (TransientSimulation) THEN
                    IF ( iter ==1 ) THEN
                       DO i=1,dim
                          Velocity(i,1:nd) = Solver % Variable % PrevValues( &
                               Solver % Variable % DOFs*(Solver % Variable % &
                               Perm(Indexes(1:nd))-1)+i,1)
                       END DO
                    ELSE
                       DO i=1,dim
                          Velocity(i,1:nd) = Solver % Variable % Values( &
                               Solver % Variable % DOFs*(Solver % Variable % &
                               Perm(Indexes(1:nd))-1)+i)
                       END DO
                    END IF
                 ELSE
                    DO i=1,dim
                       Velocity(i,1:nd) = Solver % Variable % Values( &
                            Solver % Variable % DOFs*(Solver % Variable % &
                            Perm(Indexes(1:nd))-1)+i)
                    END DO
                 END IF


                 CALL LocalMatrixBoundary( STIFF, FORCE, &
                      Velocity, rho(1), Element, nd, dim, BlockPreconditioning, ALocal )
                 CALL DefaultUpdateEquations( STIFF, FORCE )

                 IF (BlockPreconditioning) THEN
                    CALL UpdateGlobalPreconditioner( AMatrix, ALocal, nd, dim, &
                         Solver % Variable % Perm( Indexes(1:nd) ) )
                 END IF

              END IF

           END IF
        END DO
     END IF

  END IF

  CALL DefaultFinishAssembly()

  Parallel = ParEnv % Pes > 1

  IF (Parallel) THEN
     CALL DefaultDirichletBCs()
  ELSE
     CALL SetBoundaryConditions(Model, Solver % Matrix, ComponentName(Solver % Variable % name, 1), &
          1, dim+1,  Solver % Variable % Perm, Solver % Matrix % RHS)
     CALL SetBoundaryConditions(Model, Solver % Matrix, ComponentName(Solver % Variable % name, 2), &
          2, dim+1,  Solver % Variable % Perm, Solver % Matrix % RHS)     
     CALL SetBoundaryConditions(Model, Solver % Matrix, ComponentName(Solver % Variable % name, 3), &
          3, dim+1,  Solver % Variable % Perm, Solver % Matrix % RHS)          
     IF ( dim > 2 ) CALL SetBoundaryConditions(Model, Solver % Matrix, ComponentName(Solver % Variable % name, 4), &
          4, dim+1,  Solver % Variable % Perm, Solver % Matrix % RHS)               
  END IF


  atime = CPUTime() - atime
  CALL Info( 'CoupledNSUpdate', 'Assembly done', Level=4 )

  
  !-------------------------------------------------------------------------
  ! Substitute consistent splitting iterate into the search direction variable  
  !------------------------------------------------------------------------
  Snew = 0.0d0

  FlowSol => VariableGet( Model % Variables, 'VelocityTot' )
  IF( .NOT. ASSOCIATED(FlowSol) ) THEN
     CALL Fatal('CoupledNSUpdate','Could not find velocity solution VelocityTot')
  END IF

  Active = GetNOFActive()
  DO t=1,Active
     Element => GetActiveElement(t)
     !n  = GetElementNOFNodes()
     nd = GetElementDOFs( Indexes )
    
     DO i=1,nd
        j = Solver % Variable % Perm( Indexes(i) )
        k = FlowSol % Perm( Indexes(i) )
        DO p=1,dim
           Snew( (dim+1)*(j-1)+p ) = &
                FlowSol % Values( dim*(k-1)+p )
        END DO
     END DO
  END DO

  FlowSol => VariableGet( Model % Variables, 'TotalPressure' )
  IF( .NOT. ASSOCIATED(FlowSol) ) THEN
     CALL Fatal('CoupledNSUpdate','Could not pressure solution TotalPressure')
  END IF

  DO t=1,Active
     Element => GetActiveElement(t)
     !n  = GetElementNOFNodes()
     nd = GetElementDOFs( Indexes )
    
     DO i=1,nd
        j = Solver % Variable % Perm( Indexes(i) )
        k = FlowSol % Perm( Indexes(i) )

        Snew( (dim+1)*j ) = &
             FlowSol % Values( k )

     END DO
  END DO

  ! PRINT *, 'Search direction was found: Starting minimal residual update...'       

  !---------------------------------------------------------------------------------
  ! Find minimal residual update
  !---------------------------------------------------------------------------------
  IF (Parallel) THEN

    IF ( .NOT. ASSOCIATED(Solver % Matrix % ParMatrix) ) &
        CALL ParallelInitMatrix( Solver, Solver % Matrix )

    n =  Solver % Matrix % NumberOfRows
    ALLOCATE( Residual(n) )
    Residual = 0.0d0    

    CALL ParallelInitSolve( Solver % Matrix, Solver % Variable % Values, &
         Solver % Matrix % RHS, Residual )

    CALL ParallelUpdateSolve( Solver % Matrix,  & 
                     Solver % Variable % Values, Residual )

    MMatrix => ParallelMatrix( Solver % Matrix, Mx, Mb, Mr )

    n = MMatrix % NumberOfRows
    CALL GCRUpdate(n, Solver % Matrix, MMatrix, Mx, Mb, Mr, Snew, S, V, R, Round)
    CALL ParallelUpdateResult( Solver % Matrix, Solver % Variable % Values, Residual )

    DEALLOCATE( Residual )

  ELSE

    CALL Orthogonalize(Solver % Matrix % NumberOfRows, Solver % Matrix, &
        Solver % Variable % Values, Solver % Matrix % RHS, Snew, S, V, R, Round, MaxIterations)

  END IF

  !PRINT *, 'Exiting the solver...'



CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE GCRUpdate( n, A, M, x, b, r, Snew, S, V, RR, Round)
!------------------------------------------------------------------------------
    INTEGER :: n, Round
    TYPE(Matrix_t), POINTER :: A, M
    REAL(KIND=dp) :: x(:), b(:), r(:), Snew(:)  
    REAL(KIND=dp) :: S(:,:), V(:,:), RR(:)
!--------------------------------------------------------------------------------
    REAL(KIND=dp) :: T1(n), T2(n), beta, alpha, res !, Snew(n)
    INTEGER :: i,j,k
!--------------------------------------------------------------------------------

    IF ( Parallel ) CALL ParallelVector(A,Snew)   

    !j = 0
    !DO i=1,Solver % Matrix % NumberofRows
    !  k = Solver % Matrix % INVPerm(i)
    !  IF ( Solver % Matrix % ParallelInfo % Neighbourlist(k) % &
    !      Neighbours(1)==Parenv % Mype ) THEN
    !    j=j+1
    !    snew(j) = inSnew(i)
    !  END IF
    !END DO

    IF ( Round == 1) THEN
      CALL Mymv( A, x, r ) 
      r(1:n) = b(1:n)-r(1:n)
    ELSE
      r(1:n) = RR(1:n)
    END IF

    T1(1:n) = Snew(1:n) - x(1:n)
    CALL Mymv( A, T1, T2 )
    !--------------------------------------------------------------
    ! Perform the orthogonalisation of the search directions...
    !--------------------------------------------------------------
    DO i=1,Round-1
       beta = Mydot( n, V(1:n,i), T2(1:n) )
       T1(1:n) = T1(1:n) - beta * S(1:n,i)
       T2(1:n) = T2(1:n) - beta * V(1:n,i)    
    END DO

    alpha = Mynorm(n,T2)
    T1(1:n) = 1.0d0/alpha * T1(1:n)
    T2(1:n) = 1.0d0/alpha * T2(1:n)

    !-------------------------------------------------------------
    ! The update of the solution and save the search data...
    !------------------------------------------------------------- 
    beta = Mydot(n, T2, r)
    !PRINT *, 'beta = ', beta

    x(1:n) = x(1:n) + beta * T1(1:n)
    r(1:n) = r(1:n) - beta * T2(1:n)

    S(1:n,Round) = T1(1:n)
    V(1:n,Round) = T2(1:n)

    !res = ComplexNorm(m,r)/ComplexNorm(m,f)
    RR(1:n) = r(1:n)
    res = MyNorm(n,r)/Mynorm(n,b)
    !PRINT *, 'Residual norm after minimization = ', res
    WRITE(Message,'(a,I4,ES12.3)') 'OuterIteration residual for iterate', &
        Round, res
    CALL Info('NavierStokesSolver',Message,Level=4)

  !-------------------------------------------
  END SUBROUTINE GCRUpdate
  !---------------------------------------------


  !------------------------------------------------------------------------------
  SUBROUTINE Orthogonalize(n, A, x, b, Snew, S, V, R, Round, MaxRounds)
  !---------------------------------------------------------------------------
    INTEGER :: n, Round, MaxRounds
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp) :: x(n), b(n), Snew(n), &
         S(:,:), V(:,:), R(:)
    !--------------------------------------------------------------------
    INTEGER :: m, i
    REAL(KIND=dp) :: y(n),f(n), T1(n), T2(n), beta
    REAL(KIND=dp) :: res0, res, alpha
    !---------------------------------------------------------------------
    m = n
    IF ( Round == 1) THEN
       CALL Mymv( A, x, r )
       r(1:n) = b(1:n) - r(1:n)
       !res0 = ComplexNorm(m,f)
       res = Mynorm(n,r)/Mynorm(n,b)
       PRINT *, 'Residual norm for initial guess = ', res
    END IF

    !print *, 'Residual before minimization = ', res

    T1(1:n) = Snew(1:n) - x(1:n)
    CALL Mymv( A, T1, T2 )

    !--------------------------------------------------------------
    ! Perform the orthogonalisation of the search directions...
    !--------------------------------------------------------------
    DO i=1,Round-1
       beta = Mydot( n, V(1:n,i), T2(1:n) )
       T1(1:n) = T1(1:n) - beta * S(1:n,i)
       T2(1:n) = T2(1:n) - beta * V(1:n,i)    
    END DO

    alpha = Mynorm(n,T2)
    T1(1:n) = 1.0d0/alpha * T1(1:n)
    T2(1:n) = 1.0d0/alpha * T2(1:n)

    !-------------------------------------------------------------
    ! The update of the solution and save the search data...
    !------------------------------------------------------------- 
    beta = Mydot(n, T2, r)
    x(1:n) = x(1:n) + beta * T1(1:n)
    r(1:n) = r(1:n) - beta * T2(1:n)
    S(1:n,Round) = T1(1:n)
    V(1:n,Round) = T2(1:n)

    !res = ComplexNorm(m,r)/ComplexNorm(m,f)
    res = MyNorm(n,r)/Mynorm(n,b)
    WRITE(Message,'(a,I4,ES12.3)') 'OuterIteration residual for iterate', &
        Round, res
    CALL Info('NavierStokesSolver',Message,Level=4)

  !-------------------------------------------
  END SUBROUTINE Orthogonalize
  !---------------------------------------------



  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, Mass, FORCE, LOAD, Nodalrho, &
       Nodalmu, NodalVelo, Element, n, nd, dim, Stabilization, &
       Convect, GradDivParam)
  !------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: STIFF(:,:), Mass(:,:), FORCE(:), LOAD(:,:)
    REAL(KIND=dp) :: Nodalmu(:), Nodalrho(:), NodalVelo(:,:), GradDivParam
    INTEGER :: dim, n, nd
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: Stabilization, Convect
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd), dBasisdx(nd,3), &
         DetJ, LoadAtIP(dim+1), Velo(dim), Grad(dim,dim), AK, w1, w2, w3, StabTerms(nd)
    REAL(KIND=dp), POINTER :: A(:,:), F(:), M(:,:)
    LOGICAL :: Stat
    INTEGER :: t, i, j, k, l, p, q
    TYPE(GaussIntegrationPoints_t) :: IP

    REAL(KIND=dp) :: mu = 1.0d0, rho = 1.0d0, s, c, c2, ch, rotterm, mK, hK, Re, VNorm, &
         a1, a2

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !------------------------------------------------------------------------------

    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0
    Mass = 0.0d0

    !----------------------
    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )

    StabTerms = 0.0d0
    AK = 0.0d0
    ch = ElementDiameter(Element, Nodes)
    rotterm = 0.0d0

    w1 = REAL(0,dp)
    w2 = REAL(0,dp)
    w3 = REAL(0,dp)
    DO t=1,IP % n
       !--------------------------------------------------------------
       ! Basis function values & derivatives at the integration point:
       !--------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis, dBasisdx )

       s = IP % s(t) * detJ

       AK = AK + s

       Velo = MATMUL( NodalVelo(1:dim,1:nd), Basis(1:nd) )
       Grad = MATMUL( NodalVelo(1:dim,1:nd), dBasisdx(1:nd,1:dim) )

       IF (.FALSE.) THEN
          w3 = SUM( NodalVelo(2,1:nd) * dBasisdx(1:nd,1) ) - &
               SUM( NodalVelo(1,1:nd) * dBasisdx(1:nd,2) )

          IF ( (dim > 2) .AND. ( Newton .OR. Stabilization ) ) THEN
             w1 = SUM( NodalVelo(3,1:nd) * dBasisdx(1:nd,2) ) - &
                  SUM( NodalVelo(2,1:nd) * dBasisdx(1:nd,3) )             
             w2 = SUM( NodalVelo(1,1:nd) * dBasisdx(1:nd,3) ) - &
                  SUM( NodalVelo(3,1:nd) * dBasisdx(1:nd,1) )          
          END IF
       
          IF (Stabilization ) THEN
             IF ( dim > 2 ) THEN
                rotterm = rotterm + s * ( ( w2*Velo(3) - w3*Velo(2) )**2 + &
                     ( w3*Velo(1) - w1*Velo(3) )**2 + &
                     ( w1*Velo(2) - w2*Velo(1) )**2 )
             ELSE
                rotterm = rotterm + s * ( w3*w3 * (Velo(1) * Velo(1) + Velo(2) * Velo(2)) )
             END IF
          END IF
       END IF

       !Grad(1,1) = SUM( NodalVelo(1,1:nd) * dBasisdx(1:nd,1) )
       !Grad(1,2) = SUM( NodalVelo(1,1:nd) * dBasisdx(1:nd,2) )
       !Grad(2,1) = SUM( NodalVelo(2,1:nd) * dBasisdx(1:nd,1) )
       !Grad(2,2) = SUM( NodalVelo(2,1:nd) * dBasisdx(1:nd,2) )      
     
       !----------------------------------------------
       ! Material parameters at the integration point:
       !----------------------------------------------
       mu  = SUM( Basis(1:n) * Nodalmu(1:n) )
       rho = SUM( Basis(1:n) * Nodalrho(1:n) )
       !--------------------------------------------
       ! The source term at the integration point:
       !------------------------------------------
       LoadAtIP(1:dim+1) = MATMUL( Basis(1:n), LOAD(1:n,1:dim+1) )

       !----------------------------------------------------------------------------------
       ! The system matrix with only the velocity space augmented by bubbles  
       !---------------------------------------------------------------------------------
       DO p=1,nd
          DO q=1,nd
             i = (dim+1) * (p-1) + 1
             j = (dim+1) * (q-1) + 1
             A => STIFF(i:i+dim,j:j+dim)
             M => Mass(i:i+dim,j:j+dim)
             DO i=1,dim
                DO j = 1,dim
                   A(i,i) = A(i,i) + s * mu * dBasisdx(q,j) * dBasisdx(p,j)
                   A(i,j) = A(i,j) + s * mu * dBasisdx(q,i) * dBasisdx(p,j)
                   
                   IF (.FALSE.) THEN
                      A(i,j) = A(i,j) + s * 1.0d-0 * dBasisdx(q,j) * dBasisdx(p,i)
                   END IF

                END DO
                M(i,i) = M(i,i) + s * rho * Basis(p) * Basis(q)            
                !M(i,i) = M(i,i) + s * StNumber * Basis(p) * Basis(q)    
                IF ( Stabilization ) THEN
                   A(i,dim+1) = A(i,dim+1) - s * Basis(q) * dBasisdx(p,i)
                ELSE
                   ! Testing w.r.t. all velocity test functions and omitting pressure bubbles...
                   IF (q <= n) &
                        A(i,dim+1) = A(i,dim+1) - s * Basis(q) * dBasisdx(p,i)
                END IF

                IF (p <= n) &   
                     ! Testing w.r.t  standard pressure test functions...                
                     ! Pressure bubbles are constructed elsewhere... 
                     A(dim+1,i) = A(dim+1,i) - s * dBasisdx(q,i) * Basis(p)
             END DO


             ! Testing projection stabilization
             !A(dim+1,dim+1) = A(dim+1,dim+1) - s*1.0d-3/mu * Basis(p) * Basis(q)
             !A(dim+1,dim+1) = A(dim+1,dim+1) - s * ch**2 * 1.0d-2 * sum( dBasisdx(q,1:dim) * dBasisdx(p,1:dim) )
             StabTerms(p) = StabTerms(p) + s * Basis(p)
             ! End Testing...



             IF ( Convect ) THEN


                IF ( .TRUE. ) THEN
                   ! The standard convection form in 2d for testing purposes 

                   Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1 ) + &
                        rho * s * Velo(2) * dBasisdx(q,2) * Basis(p) 
                   Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1 ) + &
                        rho * s * Velo(1) * dBasisdx(q,1) * Basis(p) 
                   
                   Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2 ) + &
                        rho * s * Velo(1) * dBasisdx(q,1) * Basis(p) 
                   Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2 ) + &
                        rho * s * Velo(2) * dBasisdx(q,2) * Basis(p)


                   IF (dim > 2) THEN

                      Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1 ) + &
                           rho * s * Velo(3) * dBasisdx(q,3) * Basis(p) 
                      Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2 ) + &
                           rho * s * Velo(3) * dBasisdx(q,3) * Basis(p)                       

                      Stiff( (dim+1)*(p-1)+3, (dim+1)*(q-1)+3 ) = Stiff( (dim+1)*(p-1)+3, (dim+1)*(q-1)+3 ) + &
                           rho * s * Velo(1) * dBasisdx(q,1) * Basis(p) 
                      Stiff( (dim+1)*(p-1)+3, (dim+1)*(q-1)+3 ) = Stiff( (dim+1)*(p-1)+3, (dim+1)*(q-1)+3 ) + &
                           rho * s * Velo(2) * dBasisdx(q,2) * Basis(p)                  
                      Stiff( (dim+1)*(p-1)+3, (dim+1)*(q-1)+3 ) = Stiff( (dim+1)*(p-1)+3, (dim+1)*(q-1)+3 ) + &
                           rho * s * Velo(3) * dBasisdx(q,3) * Basis(p)                       

                   END IF




                ELSE
                   ! The convection term in the rotation form

                   IF ( dim > 2 ) THEN

                      Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1 ) + &
                           s * rho * Velo(2) * dBasisdx(q,2) * Basis(p) + &
                           s * rho * Velo(3) * dBasisdx(q,3) * Basis(p)
                      Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+2 ) - &
                           s * rho * Velo(2) * dBasisdx(q,1) * Basis(p) 
                      Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+3 ) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+3 ) - &
                           s * rho * Velo(3) * dBasisdx(q,1) * Basis(p)                   

                      Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+1 ) - &
                           s * rho * Velo(1) * dBasisdx(q,2) * Basis(p)
                      Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2 ) + &
                           s * rho * Velo(1) * dBasisdx(q,1) * Basis(p) + &
                           s * rho * Velo(3) * dBasisdx(q,3) * Basis(p)
                      Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+3 ) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+3 ) - &
                           s * rho * Velo(3) * dBasisdx(q,2) * Basis(p) 


                      Stiff( (dim+1)*(p-1)+3, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*(p-1)+3, (dim+1)*(q-1)+1 ) - &
                           s * rho * Velo(1) * dBasisdx(q,3) * Basis(p)

                      Stiff( (dim+1)*(p-1)+3, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*(p-1)+3, (dim+1)*(q-1)+2 ) - &
                           s * rho * Velo(2) * dBasisdx(q,3) * Basis(p) 

                      Stiff( (dim+1)*(p-1)+3, (dim+1)*(q-1)+3 ) = Stiff( (dim+1)*(p-1)+3, (dim+1)*(q-1)+3 ) + &
                           s * rho * Velo(2) * dBasisdx(q,2) * Basis(p) + &
                           s * rho * Velo(1) * dBasisdx(q,1) * Basis(p)                  


                      IF ( Newton ) THEN

                         Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+2 ) - &
                              s * rho * w3 * Basis(q) * Basis(p) 
                         Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+3 ) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+3 ) + &
                              s * rho * w2 * Basis(q) * Basis(p) 

                         Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+1 ) + &
                              s * rho * w3 * Basis(q) * Basis(p)
                         Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+3 ) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+3 ) - &
                              s * rho * w1 * Basis(q) * Basis(p) 

                         Stiff( (dim+1)*(p-1)+3, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*(p-1)+3, (dim+1)*(q-1)+1 ) - &
                              s * rho * w2 * Basis(q) * Basis(p)
                         Stiff( (dim+1)*(p-1)+3, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*(p-1)+3, (dim+1)*(q-1)+2 ) + &
                              s * rho * w1 * Basis(q) * Basis(p) 


                      END IF


                   ELSE

                      
                      Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1 ) + &
                           s * rho * Velo(2) * dBasisdx(q,2) * Basis(p) 
                      Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+2 ) - &
                           s * rho * Velo(2) * dBasisdx(q,1) * Basis(p) 

                      Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+1 ) - &
                           s * rho * Velo(1) * dBasisdx(q,2) * Basis(p) 
                      Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2 ) + &
                           s * rho * Velo(1) * dBasisdx(q,1) * Basis(p) 

                      
                      IF ( Newton ) THEN

                         Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+2 ) - &
                              s * rho * w3 * Basis(q) * Basis(p) 
                         Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+1 ) + &
                              s * rho * w3 * Basis(q) * Basis(p)                      
                      END IF


                      ! Testing a hybrid linearization strategy...

                      IF ( Hybrid ) THEN
                         Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1 ) - &
                              rho * s * Grad(1,1) * Basis(q) * Basis(p) 
                         Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+2 ) - &
                              rho * s * Grad(1,2) * Basis(q) * Basis(p) 
                   
                         Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+1 ) - &
                              rho * s * Grad(2,1) * Basis(q) * Basis(p) 
                         Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2 ) - &
                              rho * s * Grad(2,2) * Basis(q) * Basis(p) 
                      END IF



                   END IF

                END IF

             END IF

          END DO

          i = (dim+1) * (p-1) + 1
          F => FORCE(i:i+dim)
          F = F + s * LoadAtIP * Basis(p)

          IF ( Newton ) THEN

             F(1) = F(1) - s * rho * Basis(p) * Velo(2) * w3
             F(2) = F(2) + s * rho * Basis(p) * Velo(1) * w3                    

             ! Testing a hybrid linearization strategy...
             IF ( Hybrid .AND. (dim==2) ) THEN
                F(1) = F(1) - s * rho * Basis(p) * ( Grad(1,1) * Velo(1) + Grad(1,2) * Velo(2) )
                F(2) = F(2) - s * rho * Basis(p) * ( Grad(2,1) * Velo(1) + Grad(2,2) * Velo(2) )
             END IF
             
             IF ( dim > 2 ) THEN
                
                F(1) = F(1) + s * rho * Basis(p) * Velo(3) * w2
                F(2) = F(2) - s * rho * Basis(p) * Velo(3) * w1                                    
                F(3) = F(3) - s * rho * Basis(p) * Velo(1) * w2
                F(3) = F(3) + s * rho * Basis(p) * Velo(2) * w1        

             END IF
          END IF

       END DO
    END DO


    ! Testing projection stabilization
      DO p=1,n
        DO q=1,n
          i = (dim+1) * (p-1) + 1
          j = (dim+1) * (q-1) + 1
          A => STIFF(i:i+dim,j:j+dim)
          !A(dim+1,dim+1) = A(dim+1,dim+1) + 1.0d-3/(AK*mu) * StabTerms(p) * StabTerms(q)
        END DO
      END DO
    ! End Testing...    


    
    !-------------------------------------------------------------------------------------------
    ! The system matrix may have been allocated for the case where both velocities and pressure 
    ! are augmented by bubbles. This nullifies the effect of the pressure bubbles.
    !-------------------------------------------------------------------------------------------
    DO p = n+1,nd
       i = (dim+1) * p
       FORCE(i)   = 0.0d0
       !MASS(:,i)  = 0.0d0
       !MASS(i,:)  = 0.0d0
       STIFF(i,:) = 0.0d0
       !STIFF(:,i) = 0.0d0
       IF ( .NOT. Stabilization ) STIFF(i,i) = 1.0d0
    END DO


    IF ( .FALSE. ) THEN
!    IF ( Stabilization ) THEN
       ! This is the convection stabilization part...
       rotterm = SQRT(rotterm)/SQRT(AK)  
       ch = ch/3.0d0 * SQRT(rotterm)


       hK = element % hK
       mK = element % StabilizationMK
       

       DO t=1,IP % n
          !--------------------------------------------------------------
          ! Basis function values & derivatives at the integration point:
          !--------------------------------------------------------------
          stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
               IP % W(t),  detJ, Basis, dBasisdx )

          s = IP % s(t) * detJ

          Velo = MATMUL( NodalVelo(1:dim,1:nd), Basis(1:nd) )
          w3 = SUM( NodalVelo(2,1:nd) * dBasisdx(1:nd,1) ) - &
               SUM( NodalVelo(1,1:nd) * dBasisdx(1:nd,2) )     
          !c = -Velo(2)*dBasisdx(nd,1)+Velo(1)*dBasisdx(nd,2)

           
          !VNorm = MAX( abs(Velo(1)), abs(Velo(2)), 1.0d-12 )
          a1 = MAXVAL( ABS( NodalVelo(1,1:nd) ) )
          a2 = MAXVAL( ABS( NodalVelo(2,1:nd) ) )
          VNorm = MAX( a1, a2 )
          Re = MIN( 1.0d0, rho * mK * hK * VNorm / (4 * mu) )
          Vnorm = 1.0d0   ! The velocity norm over the entire domain 
          ch = SQRT(3.0d0/100.0d0) * SQRT( hK * Re * rho / VNorm )


          !a1 = maxval( abs( NodalVelo(1,1:nd) ) )
          !a2 = maxval( abs( NodalVelo(2,1:nd) ) )
          !VNorm = MAX( a1, a2 )
          !Re = hK * VNorm / mu
          !ch = sqrt( 2.0d-2 * hK * 2.0d0 * Re / (1 + Re) )
          

          !----------------------------------------------
          ! The construction of the pressure bubbles 
          !------------------------------------------------
          DO p=n+1,nd
             DO q=n+1,nd
                Stiff( (dim+1)*p, (dim+1)*q ) = Stiff( (dim+1)*p, (dim+1)*q ) + 1.0d0 * &                 
                     s * SUM( dBasisdx(p,1:dim) * dBasisdx(q,1:dim) )                
             END DO

             IF ( .TRUE. ) THEN
                ! The robust version of Picard linearization                
                c = -Velo(2)*dBasisdx(p,1)+Velo(1)*dBasisdx(p,2)
                DO q=1,nd

                   Stiff( (dim+1)*p, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*p, (dim+1)*(q-1)+2 ) + rho * &          
                        s * dBasisdx(q,1) * c

                   Stiff( (dim+1)*p, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*p, (dim+1)*(q-1)+1 ) - rho * &          
                        s * dBasisdx(q,2) * c
                
                END DO
             ELSE
                ! The alternative Picard
                DO q=1,nd
                   
                   Stiff( (dim+1)*p, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*p, (dim+1)*(q-1)+2 ) - w3 * &          
                        s * Basis(q) * dBasisdx(p,1)
                   
                   Stiff( (dim+1)*p, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*p, (dim+1)*(q-1)+1 ) + w3 * &          
                        s * Basis(q) * dBasisdx(p,2)
                   
                END DO
             END IF


          END DO

          !---------------------------------------------------------
          ! The projection part of the equation of motion...
          !--------------------------------------------------------

          DO p=1,nd

             IF ( .TRUE. ) THEN
                ! The robust version of Picard linearization   
                DO q=n+1,nd
                   Stiff( (dim+1)*(p-1)+1, (dim+1)*q) = Stiff( (dim+1)*(p-1)+1, (dim+1)*q) + &   
                        ch**2 * s * dBasisdx(q,1) * dBasisdx(p,2) * Velo(2)
             
                   Stiff( (dim+1)*(p-1)+1, (dim+1)*q) = Stiff( (dim+1)*(p-1)+1, (dim+1)*q) - &   
                        ch**2 * s * dBasisdx(q,2) * dBasisdx(p,2) * Velo(1)
             
                   Stiff( (dim+1)*(p-1)+2, (dim+1)*q) = Stiff( (dim+1)*(p-1)+2, (dim+1)*q) - &   
                        ch**2 * s * dBasisdx(q,1) * dBasisdx(p,1) * Velo(2)

                   Stiff( (dim+1)*(p-1)+2, (dim+1)*q) = Stiff( (dim+1)*(p-1)+2, (dim+1)*q) + &   
                        ch**2 * s * dBasisdx(q,2) * dBasisdx(p,1) * Velo(1)

                END DO

                c2 = rho * ch**2 * dBasisdx(p,2) * ( Velo(1)*Velo(1) + Velo(2)*Velo(2) )
                DO q = 1,nd

                   Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1) + &
                        s * c2 * dBasisdx(q,2)
                
                   Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+2) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+2) - &
                        s * c2 * dBasisdx(q,1)   

                END DO

                c2 = rho * ch**2 * dBasisdx(p,1) * ( Velo(1)*Velo(1) + Velo(2)*Velo(2) )
                DO q = 1,nd
                
                   Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+1) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+1) - &
                        s * c2 * dBasisdx(q,2)

                   Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2) + &
                        s * c2 * dBasisdx(q,1)   
                
                END DO
                
             ELSE
                ! The alternative Picard

                DO q=n+1,nd       
        
                   Stiff( (dim+1)*(p-1)+1, (dim+1)*q) = Stiff( (dim+1)*(p-1)+1, (dim+1)*q) + &   
                        ch**2 * w3 * s * dBasisdx(q,2) * Basis(p)

                   Stiff( (dim+1)*(p-1)+2, (dim+1)*q) = Stiff( (dim+1)*(p-1)+2, (dim+1)*q) - &   
                        ch**2 * w3 * s * dBasisdx(q,1) * Basis(p)

                END DO

                DO q = 1,nd
                   
                   Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1) + &
                        w3 * w3 * ch**2 * s * Basis(p) * Basis(q) 

                   Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2) + &
                        w3 * w3 * ch**2 * s * Basis(p) * Basis(q) 
                        
                END DO

             END IF
                
          END DO


       END DO


       IF (.FALSE.) THEN
          
          ! Condensate the pressure stabilization bubble dof
          CALL LCondensateStabilizationBubble( n, nd, dim, Stiff)          

          ! Finally nullify the effect of pressure bubbles
          
          DO q=n+1,nd
             i = (dim+1) * q
             STIFF(i,:) = 0.0d0
             STIFF(:,i) = 0.0d0
             STIFF(i,i) = 1.0d0
          END DO

       END IF
       
    END IF
  !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
  !------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LCondensateStabilizationBubble( n, nd, dim, K )
!------------------------------------------------------------------------------
    USE LinearAlgebra
!------------------------------------------------------------------------------
    INTEGER :: n, nd, dim
    REAL(KIND=dp) :: K(:,:), Kbb(nd-n,nd-n), Kbl(nd-n,(dim+1)*nd-nd+n), &
         Klb((dim+1)*nd-nd+n,nd-n)

    INTEGER :: i, Ldofs((dim+1)*nd-nd+n), Bdofs

    Ldofs = (/ (i, i=1,(dim+1)*n-1) /)
    Bdofs = n*(dim+1)

    Kbb(1,1) = K(Bdofs,Bdofs)
    !Fb(1) = F(Bdofs)
    Kbl(1,1:dim*n-1) = K(Bdofs,Ldofs)
    Klb(1:dim*n-1,1) = K(Ldofs,Bdofs)

    CALL InvertMatrix( Kbb,1 )

    !F(1:dim*n) = F(1:dim*n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    K(1:dim*n-1,1:dim*n-1) = &
         K(1:dim*n-1,1:dim*n-1) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )
!------------------------------------------------------------------------------
  END SUBROUTINE LCondensateStabilizationBubble
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBoundary( STIFF, FORCE, &
      Velocity, rho, Element, nd, dim, BlockPreconditioning, ABlock )
!------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: STIFF(:,:), FORCE(:), ABlock(:,:)
    REAL(KIND=dp) :: Velocity(:,:), rho
    INTEGER :: dim, nd
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: BlockPreconditioning
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3), &
        DetJ, LoadAtIP(dim), Velo(dim), Normal(3), SquaredVelo
    LOGICAL :: Stat
    INTEGER :: t, i, j, k, l, p, q
    TYPE(GaussIntegrationPoints_t) :: IP

    REAL(KIND=dp) :: s, c

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    STIFF = 0.0d0
    FORCE = 0.0d0
    ABlock = 0.0d0
!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    IP = GaussPoints( Element )
!------------------------------------------------------------------------------
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t),  detJ, Basis, dBasisdx )
      s = IP % s(t) * detJ
      
      Normal = Normalvector(Element, Nodes, IP % U(t), IP % V(t), .TRUE.)
      Velo = MATMUL( Velocity(1:dim,1:nd), Basis(1:nd) )
      SquaredVelo = rho * SUM( Velo(1:dim) * Velo(1:dim) ) 

      IF ( Newton ) THEN
         c = 1.0d0
      ELSE
         c = 0.5d0
      END IF

      IF ( .TRUE. ) THEN
         
         DO p=1,nd
            DO i=1,dim
               
               ! The force for the Newton iteration
               IF (Newton) &
                    FORCE( (dim+1)*(p-1)+i ) = FORCE( (dim+1)*(p-1)+i) + s * Normal(i) * Basis(p) * &
                    0.5d0 * SquaredVelo

               DO q=1,nd

                  Stiff( (dim+1)*(p-1)+i, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*(p-1)+i, (dim+1)*(q-1)+1 )  + &
                       c * s * Basis(p) * Normal(i) * Velo(1) * Basis(q) * rho

                  Stiff( (dim+1)*(p-1)+i, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*(p-1)+i, (dim+1)*(q-1)+2 )  + &
                       c * s * Basis(p) * Normal(i) * Velo(2) * Basis(q) * rho             

                  IF ( dim > 2 ) &
                       Stiff( (dim+1)*(p-1)+i, (dim+1)*(q-1)+3 ) = Stiff( (dim+1)*(p-1)+i, (dim+1)*(q-1)+3 )  + &
                       c * s * Basis(p) * Normal(i) * Velo(3) * Basis(q) * rho                 
                  
               END DO

            END DO
         END DO

         ! Change this so that entries are copied only!  
         IF ( BlockPreconditioning) THEN
            DO p=1,nd
               DO q=1,nd
                  DO i=1,dim            
                     
                     ABlock( dim*(p-1)+i, dim*(q-1)+1 ) = ABlock( dim*(p-1)+i, dim*(q-1)+1 )  + &
                          c * s * Basis(p) * Normal(i) * Velo(1) * Basis(q) * rho
                     
                     ABlock( dim*(p-1)+i, dim*(q-1)+2 ) = ABlock( dim*(p-1)+i, dim*(q-1)+2 )  + &
                          c * s * Basis(p) * Normal(i) * Velo(2) * Basis(q) * rho            

                     IF ( dim > 2 ) &
                          ABlock( dim*(p-1)+i, dim*(q-1)+3 ) = ABlock( dim*(p-1)+i, dim*(q-1)+3 )  + &
                          c * s * Basis(p) * Normal(i) * Velo(3) * Basis(q) * rho               

                  END DO
               END DO
            END DO
         END IF



      ELSE
         ! This is the explicit approximation (not in use)
         DO p=1,nd
            DO i=1,dim
               FORCE( (dim+1)*(p-1)+i ) = FORCE( (dim+1)*(p-1)+i) - s * Normal(i) * Basis(p) * &
                    0.5d0 * SquaredVelo
            END DO
         END DO
      END IF


   END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBoundary
!------------------------------------------------------------------------------









!------------------------------------------------------------------------------
      FUNCTION Mydot( n, x, y ) RESULT(s)
!------------------------------------------------------------------------------
        INTEGER :: n
        REAL(KIND=dp) :: s,x(:),y(:)
!------------------------------------------------------------------------------
        IF ( .NOT. Parallel ) THEN
          s = DOT_PRODUCT( x(1:n), y(1:n) )
        ELSE
          s = ParallelDot( n, x, y )
        END IF
!------------------------------------------------------------------------------
      END FUNCTION Mydot
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
    SUBROUTINE Mymv( A, x, b, Update )
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: x(:), b(:)
       LOGICAL, OPTIONAL :: Update
       TYPE(Matrix_t), POINTER :: A
!------------------------------------------------------------------------------
       IF ( .NOT. Parallel ) THEN
         CALL CRS_MatrixVectorMultiply( A, x, b )
       ELSE
         CALL ParallelMatrixVector( A,x,b,Update )
       END IF
!------------------------------------------------------------------------------
     END SUBROUTINE Mymv
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
      FUNCTION Mynorm( n, x ) RESULT(s)
!------------------------------------------------------------------------------
        INTEGER :: n
        REAL(KIND=dp) :: s,x(:)
!------------------------------------------------------------------------------
        IF ( .NOT. Parallel ) THEN
          s = SQRT( DOT_PRODUCT( x(1:n), x(1:n) ) )
        ELSE
          s = ParallelNorm( n, x )
        END IF
!------------------------------------------------------------------------------
      END FUNCTION Mynorm
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
   SUBROUTINE SetBoundaryConditions( Model, StiffMatrix, &
                   Name, DOF, NDOFs, Perm, rhs )
!------------------------------------------------------------------------------
!
! Set dirichlet boundary condition for given dof
!
! TYPE(Model_t) :: Model
!   INPUT: the current model structure
!
! TYPE(Matrix_t), POINTER :: StiffMatrix
!   INOUT: The global matrix
!
! CHARACTER(LEN=*) :: Name
!   INPUT: name of the dof to be set
!
! INTEGER :: DOF, NDOFs
!   INPUT: The order number of the dof and the total number of DOFs for
!          this equation
!
! INTEGER :: Perm(:)
!   INPUT: The node reordering info, this has been generated at the
!          beginning of the simulation for bandwidth optimization
!******************************************************************************
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(Matrix_t), POINTER :: StiffMatrix

    CHARACTER(LEN=*) :: Name 
    INTEGER :: DOF, NDOFs, Perm(:)
    REAL(KIND=dp), OPTIONAL :: rhs(:)    
!------------------------------------------------------------------------------

    TYPE(Element_t), POINTER :: CurrentElement
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER :: i,j,k,n,t,k1,k2,nd
    LOGICAL :: GotIt, periodic
    REAL(KIND=dp) :: Work(Model % MaxElementNodes),s

!------------------------------------------------------------------------------

    DO t = Model % NumberOfBulkElements + 1, &
        Model % NumberOfBulkElements + Model % NumberOfBoundaryElements

      CurrentElement => Model % Elements(t)
!------------------------------------------------------------------------------
!      Set the current element pointer in the model structure to
!      reflect the element being processed
!------------------------------------------------------------------------------
      Model % CurrentElement => Model % Elements(t)
!------------------------------------------------------------------------------
      !n = CurrentElement % TYPE % NumberOfNodes

      n  = GetElementNOFNodes()
      nd = GetElementDOFs( Indexes )
      NodeIndexes => CurrentElement % NodeIndexes(1:n)

      DO i=1,Model % NumberOfBCs
        IF ( CurrentElement % BoundaryInfo % Constraint == &
            Model % BCs(i) % Tag ) THEN

          Work(1:n) = ListGetReal( Model % BCs(i) % Values, &
              Name,n,NodeIndexes, gotIt )
          IF ( gotIt ) THEN
            DO j=1,n
              k = Perm(NodeIndexes(j))
              IF ( k > 0 ) THEN
                k = NDOFs * (k-1) + DOF
                CALL ZeroRow( StiffMatrix,k )
                CALL SetMatrixElement( StiffMatrix,k,k, 1.0d0 )
                IF ( PRESENT(rhs) ) rhs(k) = work(j) 
              END IF
            END DO

            DO j=n+1,nd
               k = Perm(Indexes(j))
               k = NDOFs * (k-1) + DOF              
               CALL ZeroRow( StiffMatrix,k )
               CALL SetMatrixElement( StiffMatrix,k,k, 1.0d0 )
               IF ( PRESENT(rhs) ) rhs(k) = 0.0d0  
            END DO

          END IF
        END IF
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE SetBoundaryConditions
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
END SUBROUTINE OptimalSolutionUpdate
!------------------------------------------------------------------------------
