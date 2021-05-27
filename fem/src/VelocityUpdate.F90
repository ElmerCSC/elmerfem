!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation; either
! * version 2.1 of the License, or (at your option) any later version.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! * 
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library (in file ../LGPL-2.1); if not, write 
! * to the Free Software Foundation, Inc., 51 Franklin Street, 
! * Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen, Mika  Malinen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 12 Dec 2003
! *
! *****************************************************************************/


!------------------------------------------------------------------------------
!> A consistent splitting scheme for the incompressible Navier-Stokes equations:
!> The computation of the velocity update.
!
!> It is assumed by default that this solver is used for preconditioning purposes,
!> i.e. this module is used to generate search directions for minimal residual
!> methods. In this connection the convection term is expressed in the standard
!> form.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE VelocitySolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
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
  LOGICAL :: AllocationsDone = .FALSE., Newton = .FALSE., Hybrid, PicardIteration, &
       Found, Convect, OutflowBC, &
       NewTimeStep, Stabilization, PrintOutput = .FALSE., Parallel, &
       BlockPreconditioning, ConstantBulkMatrix, ConstantBulkMatrixInUse
  TYPE(Element_t), POINTER :: Element

  CHARACTER(LEN=MAX_NAME_LEN) :: NonlinearIterationMethod, ConvectionForm, str


  INTEGER :: i, j, k, n, nb, nd, t, istat, dim, active, NonlinearIter, iter, &
       CurrentDoneTime = 0
  REAL(KIND=dp) :: Norm, PrevNorm, RelC, NonlinearTol, NonLinError, RealTime, at0
  
  INTEGER, ALLOCATABLE :: Indexes(:)
  
  TYPE(Variable_t), POINTER :: PresVar, FlowVar

  TYPE(ValueList_t), POINTER :: BodyForce, Material, BC, SolverParams
  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:,:), MASS(:,:), &
       FORCE(:), rho(:), mu(:), Velocity(:,:), Pressure(:,:), &
       PrevSol(:), NonLinRes(:)

  SAVE STIFF, LOAD, FORCE, rho, mu, Velocity, Pressure, AllocationsDone, Indexes, MASS, &
       CurrentDoneTime, PrevSol, NonLinRes
!------------------------------------------------------------------------------

  Parallel = ParEnv % PEs > 1
  SolverParams => GetSolverParams()  
  Convect = GetLogical( SolverParams, 'Convective', Found )
  IF ( .NOT. Found ) Convect = .TRUE.
  BlockPreconditioning = GetLogical( SolverParams, 'Block Preconditioning', Found )
  IF ( .NOT. Found )  BlockPreconditioning = .TRUE.

  NewTimeStep = .FALSE.
  IF (CurrentDoneTime /= Solver % DoneTime) THEN
     NewTimeStep = .TRUE.
     ! IF (NewTimeStep) PRINT *, 'New time step begins...'
     CurrentDoneTime = CurrentDoneTime + 1
  END IF

  ConstantBulkMatrix = GetLogical( SolverParams, 'Constant Bulk Matrix', Found )
  ConstantBulkMatrixInUse = ConstantBulkMatrix .AND. &
       ASSOCIATED(Solver % Matrix % BulkValues) .AND. ( .NOT. NewTimeStep )  
 


  IF ( ConstantBulkMatrix) THEN
     ! Optimize performance if ILU preconditioning is used...
     str = ListGetString( Solver % Values, &
          'Linear System Preconditioning', Found )          
     IF ( Found .AND. SEQL(str, 'ilu') ) THEN
        IF (NewTimeStep) THEN
           CALL ListAddLogical(Solver % Values, 'No Precondition Recompute', .FALSE.) 
        ELSE
           CALL ListAddLogical(Solver % Values, 'No Precondition Recompute', .TRUE.)     
        END IF
     END IF
  END IF
 

  IF (Convect) THEN

     ConvectionForm = ListGetString(Solver % Values, 'Convection Form', Found)
     IF ( .NOT. Found ) ConvectionForm = 'standard' 
     IF ( ConvectionForm /= 'rotational') THEN
        IF ( ConvectionForm == 'standard') THEN
           PicardIteration = .TRUE.
           Newton = .FALSE.
           Hybrid = .FALSE.
        ELSE
           CALL Fatal( 'VelocitySolver', 'Convection Form must be rotational or standard' )   
        END IF

     ELSE

        NonlinearIterationMethod = ListGetString(Solver % Values, &
             'Nonlinear Iteration Method', Found)
        IF ( .NOT. Found )  NonlinearIterationMethod = 'picard'

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

     END IF

     IF (BlockPreconditioning) THEN
        NonlinearIter = 1
     ELSE
        NonlinearIter = ListGetInteger( Solver % Values, &
             'Nonlinear System Max Iterations', minv=0 )
        NonlinearTol = ListGetConstReal( Solver % Values, &
             'Nonlinear System Convergence Tolerance',minv=0.0d0 )
     END IF
  ELSE
     
     NonlinearIter = 1

  END IF
   

  dim = CoordinateSystemDimension()
  Mesh => GetMesh()

  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  IF ( .NOT. AllocationsDone ) THEN
     n = dim * Mesh % MaxElementDOFs  ! just big enough for elemental arrays
     t = Mesh % MaxElementDOFs
     ALLOCATE( &
          MASS(n,n), &
          FORCE(n), &
          STIFF(n,n), &
          LOAD(t,4), &
          Indexes(t), &
          rho(t), &
          mu(t), &
          Velocity(dim,t), &
          Pressure(1,t), &
          STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( 'VelocitySolver', 'Memory allocation error.' )
     END IF

     IF ( .NOT. Parallel) THEN
        ALLOCATE(   &  
             PrevSol(Solver % Matrix % NumberOfRows), &
             NonLinRes(Solver % Matrix % NumberOfRows), &    
             STAT=istat )
        IF ( istat /= 0 ) &
             CALL Fatal( 'VelocitySolver', 'Memory allocation error.' )
     END IF
    
     AllocationsDone = .TRUE.

  END IF



  Stabilization = ListGetLogical( Solver % Values, 'Stabilize', Found )
  IF (.NOT. Found) Stabilization = .FALSE.
  ! Currently there is no option for the convection stabilization...
  Stabilization = .FALSE.
  Velocity = 0.0d0

  IF (BlockPreconditioning) THEN
     FlowVar => VariableGet( Mesh % Variables, 'Flow' )
     IF ( .NOT. ASSOCIATED(FlowVar) ) &
          CALL Fatal( 'VelocitySolver', 'The coupled flow variable Flow was not found' )     
  ELSE
     PresVar => VariableGet( Mesh % Variables, 'TotalPressure' )
     IF ( .NOT. ASSOCIATED(PresVar) ) &
          CALL Fatal( 'VelocitySolver', 'The pressure variable TotalPressure was not found' )
  END IF



  DO iter=1,NonlinearIter
     at0 = RealTime()
     !Initialize the system and do the assembly:
     !------------------------------------------
     Active = GetNOFActive()
     IF ( ConstantBulkMatrixInUse ) THEN
        Solver % Matrix % Values = Solver % Matrix % BulkValues        
        Solver % Matrix % RHS = 0.0_dp
     ELSE
        CALL DefaultInitialize()
     END IF

     DO t=1,Active 

        !IF ( RealTime() - at0 > 2.0 ) THEN
        !   WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
        !        (Solver % NumberOfActiveElements - t) / &
        !        (Solver % NumberOfActiveElements)), ' % done'
        !   CALL Info( 'VelocitySolver', Message, Level=5 )
        !   at0 = RealTime()
        !END IF

        Element => GetActiveElement(t)
        n  = GetElementNOFNodes()
        nb = GetElementNOFBDOFs()
        nd = GetElementDOFs( Indexes )

        ! Get volume forces, if any:
        !---------------------------
        BodyForce => GetBodyForce()
        LOAD = 0.0d0

        IF ( ASSOCIATED(BodyForce) ) THEN
           Load(1:n,1) = GetReal( BodyForce, 'Force 1', Found )
           Load(1:n,2) = GetReal( BodyForce, 'Force 2', Found )
           Load(1:n,3) = GetReal( BodyForce, 'Force 3', Found )
        END IF

        ! Get material parameters:
        !-------------------------
        Material => GetMaterial()
        rho(1:n) = GetReal( Material, 'Density' )
        mu(1:n)  = GetReal( Material, 'Viscosity' )

        ! Get previous elementwise pressure and velocity:
        !-------------------------------------------
        IF (BlockPreconditioning) THEN
           DO i=1,dim
              Velocity(i,1:nd) = FlowVar % PrevValues( &
                   FlowVar % DOFs * (FlowVar % &
                   Perm(Indexes(1:nd))-1)+i,1)
           END DO
           Pressure(1,1:nd) = FlowVar % Values( &
                FlowVar % DOFs*(FlowVar % &
                Perm(Indexes(1:nd))-1)+dim+1)      

        ELSE
           IF ( NewTimeStep ) THEN
              Pressure(1,1:nd) = PresVar % PrevValues( &
                   PresVar % Perm( Indexes(1:nd) ),1)
              IF (iter == 1) THEN
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

              Pressure(1,1:nd) = PresVar % Values( &
                   PresVar % Perm( Indexes(1:nd) ) )

              DO i=1,dim
                 Velocity(i,1:nd) = Solver % Variable % Values( &
                      Solver % Variable % DOFs*(Solver % Variable % &
                      Perm(Indexes(1:nd))-1)+i)
              END DO

           END IF
        END IF


        ! Get element local matrix and rhs vector:
        !-----------------------------------------
        CALL LocalMatrix( MASS, STIFF, FORCE, LOAD, rho, &
             mu,  Velocity, Pressure, Element, n, nd, nd+nb, dim, Stabilization)
      
        ! Update global matrix and rhs vector from local matrix & vector:
        !----------------------------------------------------------------
        !CALL DefaultUpdateEquations( STIFF, FORCE )
        CALL Default1stOrderTime( MASS, STIFF, FORCE )
        IF ( .NOT. ConstantBulkMatrixInUse ) THEN
           CALL DefaultUpdateEquations( STIFF, FORCE )
        ELSE
           CALL DefaultUpdateForce( FORCE ) 
        END IF

     END DO

     ! This routine is used to save the bulk values after bulk assembly 
     IF(.NOT. ConstantBulkMatrixInUse ) THEN
       CALL DefaultFinishBulkAssembly()
     END IF
    

     ! Terms on outflow boundary
     !--------------------------------------------------------------
     IF (Convect .AND. ConvectionForm == 'rotational') THEN

        Active = GetNOFBoundaryElements()
        DO t=1, Active       
           Element => GetBoundaryElement(t)
           IF ( .NOT. ActiveBoundaryElement() ) CYCLE

           n  = GetElementNOFNodes()
           nd = GetElementDOFs( Indexes )

           IF ( GetElementFamily() == 1 ) CYCLE

           BC => GetBC()
           IF ( ASSOCIATED( BC ) ) THEN
              OutflowBC= ListGetLogical( BC, 'Outflow Boundary', Found )

              IF (OutflowBC) THEN
                 IF (NewTimeStep .AND. iter == 1) THEN
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

                 CALL LocalMatrixBoundary( STIFF, FORCE, rho(1), &
                      Velocity, Element, nd, dim )
                 CALL DefaultUpdateEquations( STIFF, FORCE )
              END IF

           END IF
        END DO
     END IF

     CALL DefaultFinishAssembly()

     IF (Parallel) THEN

        CALL DefaultDirichletBCs()     
        Norm = DefaultSolve()

     ELSE
        !------------------------------------------------------------------------------------
        ! In serial runs the P2P1/Q2Q1 approximation utilizing the shape functions of
        ! p-elements may be used. The boundary conditions are therefore handled by 
        ! the following nonstandard subroutines, so that the mid-edge/face dofs
        ! are simply set to be zero. This may cause non-optimal accuracy, but this is the 
        ! way we do this. Better ways could be developed... 
        !------------------------------------------------------------------------------------
        CALL SetBoundaryConditions(Model, Solver % Matrix, ComponentName(Solver % Variable % name, 1), &
             1, dim,  Solver % Variable % Perm, Solver % Matrix % RHS)
        CALL SetBoundaryConditions(Model, Solver % Matrix, ComponentName(Solver % Variable % name, 2), &
             2, dim,  Solver % Variable % Perm, Solver % Matrix % RHS)  
        IF ( dim > 2) &   
             CALL SetBoundaryConditions(Model, Solver % Matrix, ComponentName(Solver % Variable % name, 3), &
             3, dim,  Solver % Variable % Perm, Solver % Matrix % RHS)

        IF (.NOT. Convect) THEN

           Norm = DefaultSolve()

        ELSE

           !PrevNorm = Norm
           
           !n = Solver % Matrix % NumberOfRows
           !IF ( iter == 1) THEN
           !   PrevSol(1:n) = Solver % Variable % PrevValues(1:n,1)
           !ELSE
           !   PrevSol(1:n) = Solver % Variable % Values(1:n)
           !END IF

           ! Test whether the nonlinear residual is small enough to terminate the nonlinear iteration
           !CALL CRS_MatrixVectorMultiply( Solver % Matrix, PrevSol, NonLinRes )
           !NonLinRes(1:n) = Solver % Matrix % RHS(1:n) - NonLinRes(1:n)
           !NonLinError = PNorm(n, NonLinRes)/PNorm(n, Solver % Matrix % RHS)
           ! IF ( NonLinError < 1.0d-2 ) Newton = .TRUE.
           !WRITE(Message, '(a,E12.4)') 'Nonlinear iteration residual: ', NonLinError
           !CALL Info( 'VelocitySolver', Message, Level=4)
           
           !IF (NonLinError < NonLinearTol ) THEN !.OR. iter == (NonlinearIter+1) ) THEN
           !   EXIT
           !ELSE
           !   Norm = DefaultSolve()
           !END IF
    
           Norm = DefaultSolve()

        END IF

     END IF

  END DO


CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  MASS, STIFF, FORCE, LOAD, Nodalrho, &
     Nodalmu, NodalVelo, NodalPressure, Element, n, nd, ntot, dim, StabilizedMethod)
!------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: MASS(:,:), STIFF(:,:), FORCE(:), LOAD(:,:)
    REAL(KIND=dp) :: Nodalmu(:), Nodalrho(:), NodalVelo(:,:), NodalPressure(:,:)
    INTEGER :: dim, n, nd, ntot
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: StabilizedMethod
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(ntot),dBasisdx(ntot,3), &
               DetJ,LoadAtIP(dim),Velo(dim), Grad(dim,dim), Pressure, w1, w2, w3, &
               StabTerms(dim,dim*ntot), AK, ch, rotterm
    REAL(KIND=dp), POINTER :: M(:,:), A(:,:), F(:)
    LOGICAL :: Stat
    INTEGER :: t, i, j, k, l, p, q
    TYPE(GaussIntegrationPoints_t) :: IP

    REAL(KIND=dp) :: mu = 1.0d0, rho = 1.0d0, s, c, c2

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0
    MASS  = 0.0d0

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )

    ! ch = ElementDiameter(Element, Nodes)
    ! rotterm = 0.0d0    
    ! AK = 0.0d0

    DO t=1,IP % n
       !--------------------------------------------------------------
       ! Basis function values & derivatives at the integration point:
       !--------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t),  detJ, Basis, dBasisdx )

       s = IP % s(t) * detJ
       ! AK = AK + s

       ! Material parameters at the integration point:
       !----------------------------------------------
       mu  = SUM( Basis(1:n) * Nodalmu(1:n) )
       rho = SUM( Basis(1:n) * Nodalrho(1:n) )

       ! Previous velocity & pressure at the integration point:
       !-------------------------------------------------------
       Pressure = SUM( NodalPressure(1,1:n) * Basis(1:n) )
       Velo = MATMUL( NodalVelo(1:dim,1:nd), Basis(1:nd) )
       Grad = MATMUL( NodalVelo(1:dim,1:nd), dBasisdx(1:nd,1:dim) )

       IF ( dim > 2) THEN
          w1 = SUM( NodalVelo(3,1:nd) * dBasisdx(1:nd,2) ) - &
               SUM( NodalVelo(2,1:nd) * dBasisdx(1:nd,3) )             
          w2 = SUM( NodalVelo(1,1:nd) * dBasisdx(1:nd,3) ) - &
               SUM( NodalVelo(3,1:nd) * dBasisdx(1:nd,1) )
          w3 = SUM( NodalVelo(2,1:nd) * dBasisdx(1:nd,1) ) - &
               SUM( NodalVelo(1,1:nd) * dBasisdx(1:nd,2) )         

          !rotterm = rotterm + s * ( ( w2*Velo(3) - w3*Velo(2) )**2 + &
          !     ( w3*Velo(1) - w1*Velo(3) )**2 + &
          !     ( w1*Velo(2) - w2*Velo(1) )**2 )

       ELSE
          w3 = SUM( NodalVelo(2,1:nd) * dBasisdx(1:nd,1) ) - &
               SUM( NodalVelo(1,1:nd) * dBasisdx(1:nd,2) )     

          !rotterm = rotterm + s * ( w3*w3 * (Velo(1) * Velo(1) + Velo(2) * Velo(2)) )
       END IF

       ! The source term at the integration point:
       !------------------------------------------
       LoadAtIP = MATMUL( Basis(1:n), LOAD(1:n,1:dim) )

       ! Finally, the elemental matrix & vector:
       !----------------------------------------
       DO p=1,ntot
          IF ( .NOT. ConstantBulkMatrixInUse ) THEN
             DO q=1,ntot
                i = dim * (p-1)
                j = dim * (q-1)
                A => STIFF(i+1:i+dim,j+1:j+dim)
                M => MASS(i+1:i+dim,j+1:j+dim)
                DO i=1,dim
                   M(i,i) = M(i,i) + s * rho * Basis(q) * Basis(p)
                   DO j = 1,dim
                      A(i,i) = A(i,i) + s * mu * dBasisdx(q,j) * dBasisdx(p,j)
                      A(i,j) = A(i,j) + s * mu * dBasisdx(q,i) * dBasisdx(p,j)
                   END DO
                END DO


                IF ( Convect ) THEN

                   IF ( ConvectionForm == 'standard' ) THEN

                      ! Picard iteration is the only option if the convective term is expressed in the
                      ! standard form

                      Stiff( dim*(p-1)+1, dim*(q-1)+1 ) = Stiff( dim*(p-1)+1, dim*(q-1)+1 ) + &
                           rho * s * Velo(2) * dBasisdx(q,2) * Basis(p) 
                      Stiff( dim*(p-1)+1, dim*(q-1)+1 ) = Stiff( dim*(p-1)+1, dim*(q-1)+1 ) + &
                           rho * s * Velo(1) * dBasisdx(q,1) * Basis(p) 

                      Stiff( dim*(p-1)+2, dim*(q-1)+2 ) = Stiff( dim*(p-1)+2, dim*(q-1)+2 ) + &
                           rho * s * Velo(1) * dBasisdx(q,1) * Basis(p) 
                      Stiff( dim*(p-1)+2, dim*(q-1)+2 ) = Stiff( dim*(p-1)+2, dim*(q-1)+2 ) + &
                           rho * s * Velo(2) * dBasisdx(q,2) * Basis(p)                  

                      IF (dim > 2) THEN

                         Stiff( dim*(p-1)+1, dim*(q-1)+1 ) = Stiff( dim*(p-1)+1, dim*(q-1)+1 ) + &
                              rho * s * Velo(3) * dBasisdx(q,3) * Basis(p) 
                         Stiff( dim*(p-1)+2, dim*(q-1)+2 ) = Stiff( dim*(p-1)+2, dim*(q-1)+2 ) + &
                              rho * s * Velo(3) * dBasisdx(q,3) * Basis(p)                       

                         Stiff( dim*(p-1)+3, dim*(q-1)+3 ) = Stiff( dim*(p-1)+3, dim*(q-1)+3 ) + &
                              rho * s * Velo(1) * dBasisdx(q,1) * Basis(p) 
                         Stiff( dim*(p-1)+3, dim*(q-1)+3 ) = Stiff( dim*(p-1)+3, dim*(q-1)+3 ) + &
                              rho * s * Velo(2) * dBasisdx(q,2) * Basis(p)                  
                         Stiff( dim*(p-1)+3, dim*(q-1)+3 ) = Stiff( dim*(p-1)+3, dim*(q-1)+3 ) + &
                              rho * s * Velo(3) * dBasisdx(q,3) * Basis(p)                       

                      END IF

                   ELSE

                      ! The convective terms in the rotational form

                      IF ( dim > 2) THEN

                         Stiff( dim*(p-1)+1, dim*(q-1)+1 ) = Stiff( dim*(p-1)+1, dim*(q-1)+1 ) + &
                              rho * s * Velo(2) * dBasisdx(q,2) * Basis(p) + &
                              rho * s * Velo(3) * dBasisdx(q,3) * Basis(p)
                         Stiff( dim*(p-1)+1, dim*(q-1)+2 ) = Stiff( dim*(p-1)+1, dim*(q-1)+2 ) - &
                              rho * s * Velo(2) * dBasisdx(q,1) * Basis(p) 
                         Stiff( dim*(p-1)+1, dim*(q-1)+3 ) = Stiff( dim*(p-1)+1, dim*(q-1)+3 ) - &
                              rho * s * Velo(3) * dBasisdx(q,1) * Basis(p)                   

                         Stiff( dim*(p-1)+2, dim*(q-1)+1 ) = Stiff( dim*(p-1)+2, dim*(q-1)+1 ) - &
                              rho * s * Velo(1) * dBasisdx(q,2) * Basis(p)
                         Stiff( dim*(p-1)+2, dim*(q-1)+2 ) = Stiff( dim*(p-1)+2, dim*(q-1)+2 ) + &
                              rho * s * Velo(1) * dBasisdx(q,1) * Basis(p) + &
                              rho * s * Velo(3) * dBasisdx(q,3) * Basis(p)
                         Stiff( dim*(p-1)+2, dim*(q-1)+3 ) = Stiff( dim*(p-1)+2, dim*(q-1)+3 ) - &
                              rho * s * Velo(3) * dBasisdx(q,2) * Basis(p) 


                         Stiff( dim*(p-1)+3, dim*(q-1)+1 ) = Stiff( dim*(p-1)+3, dim*(q-1)+1 ) - &
                              rho * s * Velo(1) * dBasisdx(q,3) * Basis(p)

                         Stiff( dim*(p-1)+3, dim*(q-1)+2 ) = Stiff( dim*(p-1)+3, dim*(q-1)+2 ) - &
                              rho * s * Velo(2) * dBasisdx(q,3) * Basis(p) 

                         Stiff( dim*(p-1)+3, dim*(q-1)+3 ) = Stiff( dim*(p-1)+3, dim*(q-1)+3 ) + &
                              rho * s * Velo(2) * dBasisdx(q,2) * Basis(p) + &
                              rho * s * Velo(1) * dBasisdx(q,1) * Basis(p)


                         IF ( Newton) THEN

                            Stiff( dim*(p-1)+1, dim*(q-1)+2 ) = Stiff( dim*(p-1)+1, dim*(q-1)+2 ) - &
                                 rho * s * w3 * Basis(q) * Basis(p) 
                            Stiff( dim*(p-1)+1, dim*(q-1)+3 ) = Stiff( dim*(p-1)+1, dim*(q-1)+3 ) + &
                                 rho * s * w2 * Basis(q) * Basis(p) 

                            Stiff( dim*(p-1)+2, dim*(q-1)+1 ) = Stiff( dim*(p-1)+2, dim*(q-1)+1 ) + &
                                 rho * s * w3 * Basis(q) * Basis(p)
                            Stiff( dim*(p-1)+2, dim*(q-1)+3 ) = Stiff( dim*(p-1)+2, dim*(q-1)+3 ) - &
                                 rho * s * w1 * Basis(q) * Basis(p) 

                            Stiff( dim*(p-1)+3, dim*(q-1)+1 ) = Stiff( dim*(p-1)+3, dim*(q-1)+1 ) - &
                                 rho * s * w2 * Basis(q) * Basis(p)
                            Stiff( dim*(p-1)+3, dim*(q-1)+2 ) = Stiff( dim*(p-1)+3, dim*(q-1)+2 ) + &
                                 rho * s * w1 * Basis(q) * Basis(p) 

                         END IF


                      ELSE

                         ! Nonlinear terms in rotation form, 2d-case
                         !-----------------------------------------

                         Stiff( dim*(p-1)+1, dim*(q-1)+1 ) = Stiff( dim*(p-1)+1, dim*(q-1)+1 ) + &
                              rho * s * Velo(2) * dBasisdx(q,2) * Basis(p) 
                         Stiff( dim*(p-1)+1, dim*(q-1)+2 ) = Stiff( dim*(p-1)+1, dim*(q-1)+2 ) - &
                              rho * s * Velo(2) * dBasisdx(q,1) * Basis(p) 

                         Stiff( dim*(p-1)+2, dim*(q-1)+1 ) = Stiff( dim*(p-1)+2, dim*(q-1)+1 ) - &
                              rho * s * Velo(1) * dBasisdx(q,2) * Basis(p) 
                         Stiff( dim*(p-1)+2, dim*(q-1)+2 ) = Stiff( dim*(p-1)+2, dim*(q-1)+2 ) + &
                              rho * s * Velo(1) * dBasisdx(q,1) * Basis(p) 


                         IF (Newton) THEN

                            Stiff( dim*(p-1)+1, dim*(q-1)+2 ) = Stiff( dim*(p-1)+1, dim*(q-1)+2 ) - &
                                 rho * s * w3 * Basis(q) * Basis(p) 
                            Stiff( dim*(p-1)+2, dim*(q-1)+1 ) = Stiff( dim*(p-1)+2, dim*(q-1)+1 ) +  &
                                 rho * s * w3 * Basis(q) * Basis(p) 

                         END IF

                      END IF

                   END IF

                END IF

             END DO
          ELSE
             ! Compute mass terms to generate rhs vector
             DO q=1,ntot
                i = dim * (p-1)
                j = dim * (q-1)
                M => MASS(i+1:i+dim,j+1:j+dim)
                DO i=1,dim
                   M(i,i) = M(i,i) + s * rho * Basis(q) * Basis(p)
                END DO
             END DO
          END IF
          i = dim * (p-1)
          F => FORCE(i+1:i+dim)
          F(1:dim) = F(1:dim) +  &
               s * ( LoadAtIP(1:dim) * Basis(p) + Pressure * dBasisdx(p,1:dim) )
 
          ! This should be revised for generality...
          IF (Convect .AND. Newton) THEN

             F(1) = F(1) - rho * s * Basis(p) * Velo(2) * w3
             F(2) = F(2) + rho * s * Basis(p) * Velo(1) * w3        

             IF ( dim > 2 ) THEN
                
                F(1) = F(1) + rho * s * Basis(p) * Velo(3) * w2
                F(2) = F(2) - rho * s * Basis(p) * Velo(3) * w1                                    
                F(3) = F(3) - rho * s * Basis(p) * Velo(1) * w2
                F(3) = F(3) + rho * s * Basis(p) * Velo(2) * w1        

             END IF

          END IF

       END DO
    END DO


!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBoundary( STIFF, FORCE, rho, &
      Velocity, Element, nd, dim )
!------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: STIFF(:,:), FORCE(:)
    REAL(KIND=dp) :: Velocity(:,:)
    INTEGER :: dim, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3), &
        DetJ,LoadAtIP(dim),Velo(dim), Normal(3), SquaredVelo, rho
    LOGICAL :: Stat
    INTEGER :: t, i, j, k, l, p, q
    TYPE(GaussIntegrationPoints_t) :: IP

    REAL(KIND=dp) :: s, c

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    STIFF = 0.0d0
    FORCE = 0.0d0
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
         c = 1.0d0 * rho
      ELSE
         c = 0.5d0 * rho
      END IF


      !DO p=1,nd
      !  DO i=1,dim
      !    FORCE( dim*(p-1)+i ) = FORCE( dim*(p-1)+i) - s * Normal(i) * Basis(p) * &
      !        0.5d0 * SquaredVelo
      !  END DO
      !END DO

      ! This is Picard iteration

      DO p=1,nd
         DO i=1,dim

            IF (Newton) &
                 FORCE( dim*(p-1)+i ) = FORCE( dim*(p-1)+i) + s * Normal(i) * Basis(p) * &
                    0.5d0 * SquaredVelo


            DO q=1,nd

               Stiff( dim*(p-1)+i, dim*(q-1)+1 ) = Stiff( dim*(p-1)+i, dim*(q-1)+1 )  + &
                    s * c * Basis(p) * Normal(i) * Velo(1) * Basis(q)

               Stiff( dim*(p-1)+i, dim*(q-1)+2 ) = Stiff( dim*(p-1)+i, dim*(q-1)+2 )  + &
                    s * c * Basis(p) * Normal(i) * Velo(2) * Basis(q)                  
 
               IF (dim > 2) &
                    Stiff( dim*(p-1)+i, dim*(q-1)+3 ) = Stiff( dim*(p-1)+i, dim*(q-1)+3 )  + &
                   s * c * Basis(p) * Normal(i) * Velo(3) * Basis(q)                    

            END DO

         END DO
      END DO
 
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBoundary
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
    FUNCTION Pnorm( n, x ) RESULT(s)
!------------------------------------------------------------------------------
       INTEGER :: n
       REAL(KIND=dp) :: s,x(:)
!------------------------------------------------------------------------------
       s = SQRT( DOT_PRODUCT( x(1:n), x(1:n) ) )
!------------------------------------------------------------------------------
    END FUNCTION Pnorm
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
END SUBROUTINE VelocitySolver
!------------------------------------------------------------------------------
