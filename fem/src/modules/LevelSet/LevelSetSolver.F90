!/*****************************************************************************
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
! *  Authors: Juha Ruokolainen, Peter Råback
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 16.11.2005
! *
! *****************************************************************************/
 !------------------------------------------------------------------------------
!>  Solve the advection equation for the Levelset-function using stabilization
!>  or bubbless. For the accuracy it is advisable to use 2nd order time-stepping
!>  and courant number smaller than one. 
!> \ingroup Solvers
!------------------------------------------------------------------------------
   SUBROUTINE LevelSetSolver( Model,Solver,Timestep,TransientSimulation )
!------------------------------------------------------------------------------
     USE Types
     USE DefUtils
     USE SolverUtils
     USE MaterialModels
     USE Integration

     IMPLICIT NONE
!------------------------------------------------------------------------------ 
     TYPE(Model_t), TARGET :: Model
     TYPE(Solver_t) :: Solver
     REAL(KIND=dp) :: Timestep
     LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     INTEGER :: i,j,k,n,t,iter,istat,bf_id,CoordinateSystem
 
     TYPE(Matrix_t),POINTER  :: StiffMatrix
     TYPE(Nodes_t)   :: ElementNodes
     TYPE(Element_t),POINTER :: CurrentElement
     TYPE(ValueList_t), POINTER :: Material
 
     REAL(KIND=dp) :: Norm,RelativeChange
     LOGICAL :: Stabilize = .FALSE., GotIt, AllocationsDone = .FALSE.
     INTEGER, POINTER :: SurfPerm(:), NodeIndexes(:)
     REAL(KIND=dp), POINTER :: Surface(:), ForceVector(:), Surf(:), PrevSurface(:) 
     REAL(KIND=dp), ALLOCATABLE :: LocalMassMatrix(:,:),&
       LocalStiffMatrix(:,:),Load(:),LocalForce(:),TimeForce(:)
     INTEGER :: NSDOFs,NonlinearIter,body_id,mat_id,dim
     REAL(KIND=dp) :: dt, r, ct, DsMax
     REAL(KIND=dp), ALLOCATABLE :: ElemVelo(:,:), SurfaceFlux(:)

     SAVE LocalMassMatrix,LocalStiffMatrix,Load, &
         TimeForce, LocalForce, ElementNodes,AllocationsDone, &
         ElemVelo, Surf, SurfaceFlux

#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: at,totat,st,totst
#else
     REAL(KIND=dp) :: at,totat,st,totst,CPUTime
#endif


!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------

     CALL Info( 'LevelSetSolver','-------------------------------------', Level=4 )
     CALL Info( 'LevelSetSolver', 'Solving for the levelset function', Level=4 )
     CALL Info( 'LevelSetSolver','-------------------------------------', Level=4 )

     IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

     CoordinateSystem = CurrentCoordinateSystem()
     dim = CoordinateSystemDimension()

     Surface => Solver % Variable % Values
     IF ( SIZE( Surface ) == 0 ) RETURN
     SurfPerm => Solver % Variable % Perm
     PrevSurface => Solver % Variable % PrevValues(:,1)

     StiffMatrix => Solver % Matrix
     ForceVector => StiffMatrix % RHS
     Norm = Solver % Variable % Norm

!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone ) THEN
       N = Solver % Mesh % MaxElementNodes

       ALLOCATE( ElemVelo(3, N), &
           ElementNodes % x( N ),   &
           ElementNodes % y( N ),   &
           ElementNodes % z( N ),   &
           LocalForce( 2*N ),       &
           TimeForce( 2*N ),        &
           LocalMassMatrix( 2*N,2*N ),  &
           LocalStiffMatrix( 2*N,2*N ), &
           Surf( N ), &
           SurfaceFlux( N ), &
           STAT=istat )
 
       IF ( istat /= 0 ) THEN
         CALL Fatal( 'LevelSetSolver', 'Memory allocation error.' )
       END IF

       AllocationsDone = .TRUE.
     END IF
!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------

     Stabilize = ListGetLogical( Solver % Values,'Stabilize',GotIt )

     NonlinearIter = ListGetInteger( Solver % Values, &
        'Nonlinear System Max Iterations',GotIt )
     IF ( .NOT.GotIt ) NonlinearIter = 1

     dt = Timestep

!------------------------------------------------------------------------------

     totat = 0.0d0
     totst = 0.0d0

     CALL Info( 'LevelSetSolver','-------------------------------------', Level=4 )
     CALL Info( 'LevelSetSolver', 'Updating the levelset function', Level=4 )
     CALL Info( 'LevelSetSolver','-------------------------------------', Level=4 )

     DO iter = 1, NonlinearIter
       
       at = CPUTime()

!------------------------------------------------------------------------------
       CALL DefaultInitialize()
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!      Do the assembly for bulk elements
!------------------------------------------------------------------------------

       DO t=1,Solver % NumberOfActiveElements

         CurrentElement => Solver % Mesh % Elements(Solver % ActiveElements(t))
         body_id = CurrentElement % Bodyid    
         mat_id = ListGetInteger( Model % Bodies( body_id ) % Values, 'Material' )
         Material => Model % Materials(mat_id) % Values

!------------------------------------------------------------------------------
!        Set the current element pointer in the model structure to
!        reflect the element being processed
!------------------------------------------------------------------------------
         Model % CurrentElement => CurrentElement
         n = CurrentElement % TYPE % NumberOfNodes
         NodeIndexes => CurrentElement % NodeIndexes
 
!------------------------------------------------------------------------------
!        Get element nodal coordinates
!------------------------------------------------------------------------------
         ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
         ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
         ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)

         Surf = Surface( SurfPerm(NodeIndexes) )

!------------------------------------------------------------------------------
!         Computed velocity field
!------------------------------------------------------------------------------

         ElemVelo(1,1:n) = ListGetReal( Material,'Levelset Velocity 1',n,NodeIndexes,GotIt)
         ElemVelo(2,1:n) = ListGetReal( Material,'Levelset Velocity 2',n,NodeIndexes,GotIt)
         IF(dim == 3) ElemVelo(3,1:n) = ListGetReal( Material,'Levelset Velocity 3',n,NodeIndexes,GotIt)

!------------------------------------------------------------------------------
!          Given surface flux
!------------------------------------------------------------------------------
         bf_id = ListGetInteger( Model % Bodies(CurrentElement % BodyId) % &
             Values, 'Body Force', GotIt )
         IF ( GotIt ) THEN
           SurfaceFlux(1:n) = ListGetReal( Model % BodyForces(bf_id) % Values,  &
               'Levelset Flux',n,NodeIndexes,gotIt )
         ELSE           
           SurfaceFlux(1:n) = 0.0d0
         END IF

!------------------------------------------------------------------------------
         CALL LocalMatrix( LocalMassMatrix, LocalStiffMatrix, LocalForce, Surf, &
             SurfaceFlux, ElemVelo, Stabilize, CurrentElement, n, ElementNodes )
       
!------------------------------------------------------------------------------
!        If time dependent simulation add mass matrix to stiff matrix
!------------------------------------------------------------------------------
         TimeForce = 0.0_dp
         IF ( TransientSimulation ) THEN
!------------------------------------------------------------------------------
!          NOTE: This will replace LocalStiffMatrix and LocalForce with the
!                combined information...
!------------------------------------------------------------------------------
           CALL Default1stOrderTime( LocalMassMatrix, LocalStiffMatrix, &
               LocalForce )
         END IF
!------------------------------------------------------------------------------
!      Update global matrix and rhs vector from local matrix & vector
!------------------------------------------------------------------------------
         IF ( .NOT.Stabilize ) THEN
           CALL Condensate( N, LocalStiffMatrix,  LocalForce, TimeForce )
         END IF
         
         CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce )
         
       END DO     !  Of ActiveElements
       CALL DefaultFinishBulkAssembly()

!------------------------------------------------------------------------------
!    FinishAssemebly must be called after all other assembly steps, but before
!    Dirichlet boundary settings. Actually no need to call it except for
!    transient simulations.
!------------------------------------------------------------------------------
       CALL DefaultFinishAssembly()
!------------------------------------------------------------------------------
!    Dirichlet boundary conditions
!------------------------------------------------------------------------------
       CALL DefaultDirichletBCs()

!------------------------------------------------------------------------------
     CALL Info( 'LevelSetSolver', 'Assembly done', Level=4 )
     at = CPUTime() - at
     totat = totat + at
    
!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
     st = CPUTime()

     Norm = DefaultSolve()
     RelativeChange = Solver % Variable % NonlinChange

     st = CPUTIme()-st
     totst = totst + st
      
     IF(NonlinearIter > 1) THEN
       WRITE( Message, * ) 'Iteration   : ',iter
       CALL Info( 'LevelSetSolver', Message, Level=4 )      
     END IF
     WRITE( Message, * ) 'Result Norm   : ',Norm
     CALL Info( 'LevelSetSolver', Message, Level=4 )
     WRITE( Message, * ) 'Relative Change : ',RelativeChange
     CALL Info( 'LevelSetSolver', Message, Level=4 )
     
     IF( Solver % Variable % NonlinConverged == 1 ) EXIT

!------------------------------------------------------------------------------
   END DO ! of the nonlinear iteration
!------------------------------------------------------------------------------

   WRITE(Message,'(a,F8.2)') 'Assembly done in time (s):',totat
   CALL Info( 'LevelsetSolver',Message, Level=4 )
   
   WRITE(Message,'(a,F8.2)') 'Solution done in time (s):',totst
   CALL Info( 'LevelsetSolver',Message, Level=4 )

   DsMax = MAXVAL( ABS(Surface - PrevSurface) )
   WRITE(Message,'(a,ES12.3)') 'Maximum Levelset Change',dsmax
   CALL Info( 'LevelSetSolver',Message, Level=4 )     
   CALL ListAddConstReal(Model % Simulation,'res: LevelSet Max Change',dsmax)


!------------------------------------------------------------------------------

CONTAINS

   SUBROUTINE LocalMatrix( MassMatrix,StiffMatrix,ForceVector,  &
      Surf, Flux, NodalVelo, Stabilize,Element,n,Nodes )
!------------------------------------------------------------------------------

     REAL(KIND=dp), DIMENSION(:)   :: ForceVector,Surf,Flux
     REAL(KIND=dp), DIMENSION(:,:) :: MassMatrix,StiffMatrix,NodalVelo
     LOGICAL :: Stabilize
     INTEGER :: n
     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t), POINTER :: Element

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
!
     REAL(KIND=dp) :: ddBasisddx(n,3,3)
     REAL(KIND=dp) :: Basis(2*n)
     REAL(KIND=dp) :: dBasisdx(2*n,3),SQRTElementMetric
     REAL(KIND=dp) :: Velo(3),Force,Grad(3),GradAbs
     REAL(KIND=dp) :: A,M,Load, FL
     REAL(KIND=dp) :: VNorm,hK,mK
     REAL(KIND=dp) :: Lambda=1.0,Pe,Pe1,Pe2,Tau,x,y,z

     INTEGER :: i,j,k,c,p,q,t,dim,N_Integ,NBasis
     REAL(KIND=dp) :: s,u,v,w
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     REAL(KIND=dp) :: SU(n),SW(n)
     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     LOGICAL :: stat,Bubbles

!------------------------------------------------------------------------------

     dim = CoordinateSystemDimension()
     c = dim + 1

     ForceVector = 0.0D0
     StiffMatrix = 0.0D0
     MassMatrix  = 0.0D0
     Load = 0.0d0
     Grad = 0.0d0

     NBasis = n
     Bubbles = .FALSE.
     IF ( .NOT. Stabilize ) THEN
       NBasis = 2*n
       Bubbles = .TRUE.
     END IF

!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IF ( Bubbles ) THEN
        IntegStuff = GaussPoints( element, Element % TYPE % GaussPoints2 )
     ELSE
        IntegStuff = GaussPoints( element )
     END IF
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!    Stabilization parameters: hK, mK (take a look at Franca et.al.)
!    If there is no convection term we don t need stabilization.
!------------------------------------------------------------------------------
     IF ( Stabilize ) THEN
       hK = element % hK
       mK = element % StabilizationMK
     END IF

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
       stat = ElementInfo( Element,Nodes,u,v,w,SQRTElementMetric, &
           Basis,dBasisdx,ddBasisddx,Stabilize,Bubbles )

       s = SQRTElementMetric * S_Integ(t)
!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
       IF ( Coordinates == AxisSymmetric .OR. &
           Coordinates == CylindricSymmetric ) THEN
         s = s * SUM( Nodes % x(1:n)*Basis(1:n) )
       END IF

!------------------------------------------------------------------------------
         
       FL = SUM( Flux(1:n) * Basis(1:n) )
       
       DO i=1,dim
         Grad(i) = SUM( dBasisdx(1:n,i) * Surf(1:n) )
       END DO
       GradAbs = SQRT( SUM( Grad(1:dim) * Grad(1:dim) ) )
       IF ( GradAbs > 10*AEPS ) THEN
         Grad = Grad / GradAbs
       END IF

!------------------------------------------------------------------------------
!      Velocity from previous iteration at the integration point
!------------------------------------------------------------------------------

       Velo = 0.0D0
       Velo(1) = SUM( NodalVelo(1,1:n)*Basis(1:n) )
       Velo(2) = SUM( NodalVelo(2,1:n)*Basis(1:n) )
       IF ( dim > 2 ) Velo(3) = SUM( NodalVelo(3,1:n)*Basis(1:n) )


       IF ( Stabilize ) THEN
!------------------------------------------------------------------------------
!           Stabilization parameter Tau
!------------------------------------------------------------------------------
         VNorm = SQRT( SUM(Velo(1:dim)**2) ) + FL
         Pe  = 1.0d0
         Tau = 0.0D0
         IF ( VNorm /= 0.0 ) THEN
           Tau = hK * Pe / (2 * VNorm)
           Tau = 1.0d0 / SQRT( (2.0d0 / dt)**2 + 1.0d0/Tau**2 )
         ELSE
           Tau = dt / 2.0d0 
         END IF

!------------------------------------------------------------------------------
!           Compute residual & stablization vectors
!------------------------------------------------------------------------------
         DO p=1,N
           SU(p) = 0.0d0
           DO i = 1,dim
             SU(p) = SU(p) + dBasisdx(p,i) * (Velo(i) + FL * Grad(i) )
           END DO

           SW(p) = 0.0d0
           DO i = 1,dim
             SW(p) = SW(p) + dBasisdx(p,i) * (Velo(i) + FL * Grad(i) )
           END DO
         END DO
       END IF


!------------------------------------------------------------------------------
!      Loop over basis functions of both unknowns and weights
!------------------------------------------------------------------------------
       DO p=1,NBasis
         DO q=1,NBasis
!------------------------------------------------------------------------------
!         The diffusion term
!------------------------------------------------------------------------------

           M = Basis(q) * Basis(p)
           A = 0.0d0

           IF(GradAbs > AEPS) THEN
             DO i=1,dim
               A = A + FL * (Grad(i) / GradAbs) * dBasisdx(q,i) * Basis(p)
             END DO
           END IF

!------------------------------------------------------------------------------
!           The convection term
!------------------------------------------------------------------------------

           DO i=1,dim
             A = A + Velo(i) * dBasisdx(q,i) * Basis(p)
           END DO
!------------------------------------------------------------------------------
!           Next we add the stabilization...
!------------------------------------------------------------------------------
           IF ( Stabilize ) THEN
             A = A + Tau * SU(q) * SW(p)
             M = M + Tau * Basis(q) * SW(p)
           END IF

           StiffMatrix(p,q) = StiffMatrix(p,q) + s * A
           MassMatrix(p,q)  = MassMatrix(p,q)  + s * M

         END DO
       END DO

!------------------------------------------------------------------------------
!      The righthand side...
!------------------------------------------------------------------------------
       Force = 0.0d0

       DO p=1,NBasis
          Load = Basis(p)
          IF ( Stabilize ) Load = Load + Tau * SW(p)
          ForceVector(p) = ForceVector(p) + s * Force * Load
       END DO
     END DO

!------------------------------------------------------------------------------
   END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 END SUBROUTINE LevelSetSolver
!------------------------------------------------------------------------------


