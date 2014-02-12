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
! *  Subroutines for moving and renormalizing the level set function, 
! *  and computing the curvature of the level set field. 
! *
! ******************************************************************************
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

     REAL(KIND=dp) :: at,totat,st,totst,CPUTime


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


!------------------------------------------------------------------------------
!>  Renormalizes the levelset function using straight-forward geometric search
!>  Also includes an option to do the covection at the same time as an alternavtive
!>  for using a separate solver for the advection.
!> \ingroup Solvers
!------------------------------------------------------------------------------
   SUBROUTINE LevelSetDistance( Model,Solver,Timestep,TransientSimulation )
!------------------------------------------------------------------------------
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
     TYPE(Nodes_t)   :: ElementNodes
     TYPE(Element_t),POINTER :: CurrentElement, Element
     TYPE(ValueList_t), POINTER :: Material
     TYPE(Variable_t), POINTER :: SurfSol, DistanceSol
 
     INTEGER :: i,j,k,l,n,t,iter,istat,body_id,mat_id,bf_id,&
          CoordinateSystem, TimesVisited = 0
     REAL(KIND=dp) :: Norm,x,y,s,a,b,c,d,x0,x1,y0,y1
     INTEGER, POINTER :: NodeIndexes(:)
     LOGICAL :: GotIt, Convect, ExtractAllocated = .FALSE., DistanceAllocated=.FALSE.
     INTEGER, POINTER :: SurfPerm(:)
     REAL(KIND=dp), POINTER :: Surface(:),Distance(:), Surf(:)
     REAL(KIND=dp), ALLOCATABLE :: ZeroNodes(:,:,:), Direction(:)
     INTEGER, POINTER :: DistancePerm(:)
     INTEGER :: ZeroLevels, ReinitializeInterval, ExtractInterval
     LOGICAL :: Reinitialize, Extrct
     REAL(KIND=dp) :: Relax, dt, r, NarrowBand, DsMax
     REAL(KIND=dp), ALLOCATABLE :: ElemVelo(:,:), SurfaceFlux(:)
     REAL(KIND=dp) :: at,totat,st,totst,CPUTime
     CHARACTER(LEN=MAX_NAME_LEN) :: LevelSetVariableName

     SAVE ElementNodes, ElemVelo, Direction, ZeroNodes, TimesVisited, &
         Distance, DistancePerm, ExtractAllocated, DistanceAllocated

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------

     TimesVisited = TimesVisited + 1
     ReinitializeInterval = ListGetInteger(Solver % Values,&
          'Reinitialize Interval',GotIt) 
     IF(.NOT. GotIt) ReinitializeInterval = 1

     ExtractInterval = ListGetInteger(Solver % Values,&
         'Extract Interval',GotIt) 
     IF(.NOT. GotIt) ExtractInterval = ReinitializeInterval
     
     IF( ReinitializeInterval == 0) THEN
       Reinitialize = .FALSE.
     ELSE       
       Reinitialize = ( MOD(TimesVisited, ReinitializeInterval) == 0 )
     END IF

     IF( ExtractInterval == 0) THEN
       Extrct = Reinitialize
     ELSE
       Extrct = Reinitialize .OR. ( MOD(TimesVisited, ExtractInterval) == 0 )
     END IF
     
     IF(.NOT. Extrct) THEN
       CALL Info( 'LevelSetDistance','--------------------------------------', Level=4 )
       CALL Info( 'LevelSetDistance','Doing nothing this time', Level=4 )
       CALL Info( 'LevelSetDistance','--------------------------------------', Level=4 )          
       RETURN
     END IF

     ! The variable that should be reinitialized
     LevelSetVariableName = ListGetString(Solver % Values,'LevelSet Variable',GotIt) 
     IF(.NOT. GotIT) LevelSetVariableName = 'Surface'
     SurfSol => VariableGet( Solver % Mesh % Variables, TRIM(LevelSetVariableName) )
     IF(ASSOCIATED(SurfSol)) THEN
       Surface  => SurfSol % Values
       SurfPerm => SurfSol % Perm
     ELSE
       CALL Warn('LevelSetDistance','SurfSol does not exist: '//TRIM(LevelSetVariableName))
       RETURN
     END IF

     CoordinateSystem = CurrentCoordinateSystem() 
     Convect = ListGetLogical(Solver % Values,'Levelset Convect',GotIt)
     NarrowBand = ListGetConstReal(Solver % Values,'Narrow Band',GotIt)
     IF(.NOT. GotIt) NarrowBand = HUGE(NarrowBand)
     dsMax = 0.0d0
     dt = Timestep
     
 
!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. ExtractAllocated ) THEN
       N = Solver % Mesh % MaxElementNodes
       ALLOCATE( ElementNodes % x( N ), ElementNodes % y( N ), ElementNodes % z( N ),   &
           ElemVelo( 2, N), ZeroNodes(Solver % Mesh % NumberOfBulkElements,2,2), &
           STAT=istat )
       IF ( istat /= 0 ) THEN
         CALL Fatal( 'LevelSetDistance', 'Memory allocation error 1.' )
       END IF
       IF( Convect ) THEN
         ALLOCATE( Direction( Solver % Mesh % NumberOfBulkElements), STAT=istat )
       END IF
       IF ( istat /= 0 ) THEN
         CALL Fatal( 'RerormalizeSolver', 'Memory allocation error 2.' )
       END IF
       ExtractAllocated = .TRUE.
     END IF

!------------------------------------------------------------------------------
!    Extract the zero levelset 
!------------------------------------------------------------------------------

     CALL Info( 'LevelSetDistance','--------------------------------------', Level=4 )
     CALL Info( 'LevelSetDistance','Extracting the zero levelset', Level=4 )
     CALL Info( 'LevelSetDistance','--------------------------------------', Level=4 )

     st = CPUTime()

     CALL ExtractZeroLevel()

     st = CPUTIme()-st
     WRITE(Message,'(a,F8.2)') 'Zero level extracted in time (s):',st
     CALL Info( 'LevelSetDistance',Message, Level=4 )

     IF( ZeroLevels == 0) THEN
       CALL Warn('LevelSetDistance','The does not seem to be a zero level-set present, exiting...')
       RETURN
     END IF

     IF(.NOT. Reinitialize) THEN
       CALL Info('LevelSetDistance','Exiting without reinitialization')
       RETURN
     END IF

     CALL Info( 'LevelSetDistance','--------------------------------------', Level=4 )
     CALL Info( 'LevelSetDistance','Computing the signed distance function', Level=4 )
     CALL Info( 'LevelSetDistance','--------------------------------------', Level=4 )

!------------------------------------------------------------------------------
!    Allocate some permanent storage for computing the signed distance
!------------------------------------------------------------------------------
     IF ( .NOT. DistanceAllocated ) THEN
       
       ! The variable for computing the distance
       DistanceSol => Solver % Variable
       IF(ASSOCIATED(DistanceSol)) THEN
         DistancePerm => DistanceSol % Perm
         Distance => DistanceSol % Values
       ELSE
         ALLOCATE(Distance (SIZE(Surface)),STAT=istat)
         IF ( istat /= 0 ) THEN
           CALL Fatal( 'LevelSetDistance', 'Memory allocation error 1.' )
         END IF
         DistancePerm => SurfPerm
       END IF
       DistanceAllocated = .TRUE.
     END IF

!------------------------------------------------------------------------------
!    Compute the signed distance
!------------------------------------------------------------------------------
     st = CPUTIme()
     IF(Convect) THEN
       DO i=1,Solver % Mesh % NumberOfNodes
         Distance(DistancePerm(i)) =  ComputeDistanceWithDirection( Solver % Mesh % Nodes % x(i), &
             Solver % Mesh % Nodes % y(i), 0.0d0, Surface(i) )
       END DO
     ELSE
       DO i=1,Solver % Mesh % NumberOfNodes
         Distance(DistancePerm(i)) =  ComputeDistance( Solver % Mesh % Nodes % x(i), &
             Solver % Mesh % Nodes % y(i), 0.0d0 )
       END DO
       WHERE( Surface < 0 ) Distance = -Distance
     END IF

     n = Solver % Mesh % NumberOfNodes
     Solver % Variable % Norm = SQRT( SUM(Distance**2)/n )

!------------------------------------------------------------------------------
!    Apply the reinitialization to the primary levelset field
!------------------------------------------------------------------------------

     IF( ListGetLogical(Solver % Values,'Reinitialize Passive',GotIt) ) THEN
       CALL Info('LevelSetDistance','Reinitialization not applied to Levelset function')
     ELSE
       ! Update also the previous timesteps so that the differentials remain
       ! unchanged. Otherwise spurious effects are introduced. 
       IF(ASSOCIATED(SurfSol % PrevValues)) THEN        
         j = MIN(2, SIZE(SurfSol % PrevValues,2) )
         IF( ReinitializeInterval > j) THEN                 
           DO i=1,j
             SurfSol % PrevValues(:,i) = SurfSol % PrevValues(:,i) + Distance - Surface
           END DO
         END IF
       END IF
       Surface = Distance
     END IF

     st = CPUTIme()-st
     WRITE(Message,'(a,F8.2)') 'Reinitialization done in time (s):',st
     CALL Info( 'LevelSetNormalize',Message, Level=4 )
 
     IF(Convect) THEN
       WRITE(Message,'(a,ES12.3)') 'Maximum Levelset Change',dsmax
       CALL Info( 'LevelSetDistance',Message, Level=4 )     
       CALL ListAddConstReal(Model % Simulation,'res: LevelSet Max Change',dsmax)
     END IF

!------------------------------------------------------------------------------

CONTAINS


!------------------------------------------------------------------------------
!> Extract the zero levelset as defined by the levelset function.
!------------------------------------------------------------------------------
  SUBROUTINE ExtractZeroLevel()
!------------------------------------------------------------------------------

    INTEGER :: i,j,k,l,m,n,div,onetwo,corners, &
        TriangleIndexes(3), mat_id, body_id, LocalInd(3)
    TYPE(Variable_t), POINTER :: Var
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: x0,y0,x1,y1,t,nx(3),ny(3),nz(3),srf(3),w0(3),w1(3),fval,dt, &
        r1x, r1y, r2x, r2y, aid, Maxsrf, ds1, ds2
    LOGICAL :: FileCreated = .FALSE., FileAppend, GotIt, FileSave, Found, FileNumber
    TYPE(ValueList_t),POINTER :: Material
    CHARACTER(LEN=MAX_NAME_LEN) :: Filename, Filename2
    INTEGER :: VisitedTimes = 0, NumberOfFields=0
    TYPE(ValueList_t), POINTER :: Params

    SAVE VisitedTimes, FileAppend, NumberOfFields
!------------------------------------------------------------------------------

    VisitedTimes = VisitedTimes + 1
    Params => Solver % Values
    Filename = ListGetString(Params,'Filename',FileSave )

    IF(FileSave) THEN         
      FileNumber = ListGetLogical(Params,'Filename Numbering',GotIt)     
      FileAppend = ListGetLogical(Params,'File Append',GotIt)

      IF( FileNumber ) THEN
        WRITE( Filename,'(A,I0)') TRIM(Filename),VisitedTimes
        OPEN (10,FILE=Filename)
      ELSE IF(FileAppend .AND. VisitedTimes > 1) THEN 
        OPEN (10, FILE=Filename, POSITION='APPEND')
      ELSE 
        OPEN (10,FILE=Filename)
      END IF      
    END IF
    
    Surface => SurfSol % Values
    SurfPerm => SurfSol % Perm
    
    ZeroLevels = 0
    DO i=1,Solver % Mesh % NumberOfBulkElements
      
      Element => Solver % Mesh % Elements(i)
      n = Element % TYPE % NumberOfNodes
      NodeIndexes => Element % NodeIndexes
      
      IF ( ALL( Surface(SurfPerm(NodeIndexes)) < 0) .OR. &
          ALL( Surface(SurfPerm(NodeIndexes)) > 0) ) CYCLE
      
      corners = Element % TYPE % ElementCode / 100 
      IF(corners < 3 .OR. corners > 4) THEN
        CALL Warn('ExtractZeroLevel','Implemented only for triangles and quads')
      END IF

      IF( Convect ) THEN
        Model % CurrentElement => Element

        body_id = Element % BodyId
        mat_id = ListGetInteger( Model % Bodies(body_id) % Values, 'Material', &
            minv=1, maxv=Model % NumberOfMaterials )
        Material => Model % Materials(mat_id) % Values

        ElemVelo(1,1:n) = ListGetReal( Material,'Levelset Velocity 1',n,NodeIndexes,GotIt)
        ElemVelo(2,1:n) = ListGetReal( Material,'Levelset Velocity 2',n,NodeIndexes,GotIt)
      END IF


      DO div = 1,corners-2

        SELECT CASE (corners) 
        CASE (3)
          LocalInd(1) = 1
          LocalInd(2) = 2
          LocalInd(3) = 3
          
        CASE(4)               
          IF(div == 1) THEN
            LocalInd(1) = 1 
            LocalInd(2) = 2
            LocalInd(3) = 4
          ELSE
            LocalInd(1) = 2
            LocalInd(2) = 3
            LocalInd(3) = 4
          END IF

        END SELECT
        
        TriangleIndexes = NodeIndexes(LocalInd)
        srf = Surface(SurfPerm(TriangleIndexes))                       
        IF ( ALL(srf < 0) .OR. ALL( srf > 0) ) CYCLE

        nx = Solver % Mesh % Nodes % x(TriangleIndexes)
        ny = Solver % Mesh % Nodes % y(TriangleIndexes)
        nz = Solver % Mesh % Nodes % z(TriangleIndexes)

        CALL TriangleIsoLineWeights( nx,ny,nz,srf,w0,w1,found)
        IF( .NOT. Found) CYCLE

        x0 = SUM(w0 * nx)
        y0 = SUM(w0 * ny)
        x1 = SUM(w1 * nx)
        y1 = SUM(w1 * ny)

        r1x = x1 - x0
        r1y = y1 - y0 

        ds1 = SQRT( r1x*r1x + r1y*r1y)
        IF(ds1 < AEPS) CYCLE

        ZeroLevels = ZeroLevels + 1

        IF(Convect) THEN

          ! Find the value that differs most from the zero levelset
          j = 1
          MaxSrf = srf(1)
          DO k=2,3
            IF(ABS(srf(k)) > ABS(Maxsrf)) THEN
              Maxsrf = srf(k)
              j = k
            END IF
          END DO
          
          r2x = nx(j) - x0
          r2y = ny(j) - y0
          
          aid = r1x * r2y - r2x * r1y
          ds2 = SQRT( r2x*r2x + r2y*r2y)
          
          Direction( ZeroLevels ) = aid / (ds1 * ds2)
          IF( Maxsrf < 0.0) THEN
            Direction( ZeroLevels ) = -Direction( ZeroLevels )
          END IF
          
          r1x = SUM(w0 * ElemVelo(1,LocalInd)) * dt
          r1y = SUM(w0 * ElemVelo(2,LocalInd)) * dt
          r2x = SUM(w1 * ElemVelo(1,LocalInd)) * dt
          r2y = SUM(w1 * ElemVelo(2,LocalInd)) * dt

          ds1 = SQRT( r1x*r1x + r1y*r1y)
          ds2 = SQRT( r2x*r2x + r2y*r2y)
          dsmax = MAX(dsmax, MAX(ds1, ds2) )
            
          ZeroNodes(ZeroLevels,1,1) = x0 + r1x
          ZeroNodes(ZeroLevels,1,2) = y0 + r1y
          ZeroNodes(ZeroLevels,2,1) = x1 + r2x
          ZeroNodes(ZeroLevels,2,2) = y1 + r2y
        ELSE
          ZeroNodes(ZeroLevels,1,1) = x0
          ZeroNodes(ZeroLevels,1,2) = y0
          ZeroNodes(ZeroLevels,2,1) = x1
          ZeroNodes(ZeroLevels,2,2) = y1
        END IF


        IF(FileSave) THEN
          
          DO onetwo = 1,2
            
            IF( FileAppend ) THEN
              WRITE(10,'(I4)',ADVANCE='NO') Solver % DoneTime
            END IF            

	    m = 0
            Var => Model % Variables
            DO WHILE( ASSOCIATED( Var ) )
              
              IF ( .NOT. Var % Output .OR. SIZE(Var % Values) == 1 .OR. (Var % DOFs /= 1) ) THEN
                Var => Var % Next        
                CYCLE
              END IF
	      m = m + 1              

              fval = 0.0d0
              DO k=1,3
                l = TriangleIndexes(k)
                IF ( ASSOCIATED(Var % Perm) ) l = Var % Perm(l)
                IF(l > 0) THEN
                  IF(onetwo == 1) THEN
                    fval = fval + w0(k) * (Var % Values(l))
                  ELSE
                    fval = fval + w1(k) * (Var % Values(l))
                  END IF
                END IF
              END DO
              
              WRITE(10,'(ES20.11E3)',ADVANCE='NO') fval
              Var => Var % Next          
            END DO
            WRITE(10,'(A)') ' '
            
          END DO
	END IF        

      END DO

    END DO ! of elements

    


    IF(FileSave) THEN
      CLOSE(10)

      IF( m /= NumberOfFields ) THEN
	IF( NumberOfFields > 0 ) THEN
  	  CALL Warn('ExtractZeroLevel','Mismacth in number of fields')
          PRINT *,j,' vs. ',NumberOfFields
	END IF
	NumberOfFields = m

        OPEN (10, FILE=TRIM(Filename)//TRIM(".names") )
        WRITE(10,'(A,A)') 'Variables in file: ',TRIM(Filename)
        j = 1
        WRITE(10,'(I3,": ",A)') j,'timestep'
        
        Var => Model % Variables
        DO WHILE( ASSOCIATED( Var ) )          
          IF ( .NOT. Var % Output .OR. SIZE(Var % Values) == 1 .OR. (Var % DOFs /= 1) ) THEN
            Var => Var % Next        
            CYCLE 
          END IF          
          j = j + 1
          WRITE(10,'(I3,": ",A)') j,TRIM(Var % Name)
          Var => Var % Next          
        END DO
        CLOSE(10)
      END IF  
    END  IF

!------------------------------------------------------------------------------
   END SUBROUTINE ExtractZeroLevel
!------------------------------------------------------------------------------
 

!------------------------------------------------------------------------------
!> This subroutine extracts the zero line of one triangular element. 
!------------------------------------------------------------------------------
   SUBROUTINE TriangleIsoLineWeights( NX,NY,NZ,S,w0,w1,Found )
!------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=dp) :: NX(:),NY(:),NZ(:),S(:),w0(3),w1(3)
      LOGICAL :: Found
      REAL(KIND=dp) :: t
!------------------------------------------------------------------------------

      Found = .TRUE.

      w0 = 0.0d0
      w1 = 0.0d0

      IF ( ABS(S(1)) < AEPS .AND. ABS(S(2)) < AEPS ) THEN
        w0(1) = 1.0d0
        w1(2) = 1.0d0
      ELSE IF ( ABS(S(1)) < AEPS .AND. ABS(S(3)) < AEPS ) THEN
        w0(1) = 1.0d0
        w1(3) = 1.0d0
      ELSE IF ( ABS(S(2)) < AEPS .AND. ABS(S(3)) < AEPS ) THEN
        w0(2) = 1.0d0
        w1(3) = 1.0d0
      ELSE IF ( ALL(S <= 0) .OR. ALL( S >= 0) ) THEN
        Found = .FALSE.
      ELSE
        IF ( S(1) >= 0 .AND. S(2) >= 0 .OR. &
            S(1) <= 0 .AND. S(2) <= 0 ) THEN          
          t = -S(1) / ( S(3) - S(1) )
          w0(3) = t
          w0(1) = 1-t
          t = -S(2) / ( S(3) - S(2) )
          w1(3) = t
          w1(2) = 1-t
        ELSE IF ( S(1) >= 0 .AND. S(3) >= 0 .OR. &
            S(1) <= 0 .AND. S(3) <= 0 ) THEN
          t = -S(1) / ( S(2) - S(1) )
          w0(2) = t
          w0(1) = 1-t
          t = -S(3) / ( S(2) - S(3) )
          w1(2) = t
          w1(3) = 1-t
          
        ELSE IF ( S(2) >= 0 .AND. S(3) >= 0 .OR. &
            S(2) <= 0 .AND. S(3) <= 0 ) THEN          
          t = -S(2) / ( S(1) - S(2) )
          w0(1) = t
          w0(2) = 1-t
          t = -S(3) / ( S(1) - S(3) )
          w1(1) = t
          w1(3) = 1-t
        ELSE 
          PRINT *,'TriangleIsoLineWeights: this should not occur'
          PRINT *,s(1),s(2),s(3)
          STOP
        END IF
      END IF
!------------------------------------------------------------------------------
    END SUBROUTINE TriangleIsoLineWeights
!------------------------------------------------------------------------------

 
!------------------------------------------------------------------------------
!> Computes the distance from the given zero levelset given by ZeroNodes. 
!------------------------------------------------------------------------------
   FUNCTION ComputeDistance(xp,yp,zp) RESULT(dist)
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: xp,yp,zp,dist
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: x0,y0,x1,y1,a,b,c,d,s
     INTEGER :: i,j,k,n
!------------------------------------------------------------------------------
     dist = HUGE(dist)
     DO i=1,ZeroLevels
       x0 = ZeroNodes(i,1,1)
       y0 = ZeroNodes(i,1,2)
       
       x1 = ZeroNodes(i,2,1)
       y1 = ZeroNodes(i,2,2)
       
       a = xp - x0
       b = x0 - x1
       d = y0 - y1
       c = yp - y0
       s = b**2 + d**2
       
       x = x0
       y = y0
       IF ( s > 10*AEPS ) THEN
         s = MIN( MAX( -(a*b + c*d) / s, 0.0d0), 1.0d0 )
         x = (1-s) * x0 + s * x1
         y = (1-s) * y0 + s * y1
       END IF
       
       dist = MIN( dist, SQRT( (xp - x)**2 + (yp - y)**2 ) )
     END DO
!------------------------------------------------------------------------------
   END FUNCTION ComputeDistance
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Computes the signed distance from the given zero levelset given by ZeroNodes. 
!------------------------------------------------------------------------------
   FUNCTION ComputeDistanceWithDirection(xp,yp,zp,prevdist) RESULT(mindist)
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: xp,yp,zp,prevdist,mindist
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: x0,y0,x1,y1,a,b,c,d,s,dist,r1x,r1y,r2x,r2y,angle,angle0
     INTEGER :: i,j,k,n
!------------------------------------------------------------------------------

     IF(prevdist - dsmax > narrowband) THEN
       mindist = prevdist - dsmax
       RETURN
     ELSE IF(prevdist + dsmax < -narrowband) THEN
       mindist = prevdist + dsmax
       RETURN
     END IF

     mindist = HUGE(mindist)

     DO i=1,ZeroLevels
       x0 = ZeroNodes(i,1,1)
       y0 = ZeroNodes(i,1,2)
       
       x1 = ZeroNodes(i,2,1)
       y1 = ZeroNodes(i,2,2)
       
       a = xp - x0
       b = x0 - x1
       d = y0 - y1
       c = yp - y0
       s = b**2 + d**2
       
       x = x0
       y = y0
       IF ( s > 10*AEPS ) THEN
         s = MIN( MAX( -(a*b + c*d) / s, 0.0d0), 1.0d0 )
         x = (1-s) * x0 + s * x1
         y = (1-s) * y0 + s * y1
       END IF
       
       dist = SQRT( (xp - x)**2 + (yp - y)**2 ) 
       

       IF(dist <= (ABS(mindist) + AEPS) ) THEN
         
         r1x = x1 - x0
         r1y = y1 - y0
         r2x = xp - x0
         r2y = yp - y0
         
         angle = r1x * r2y - r2x * r1y
         
         ! Favor parents with clear angles
         IF( dist < (ABS(mindist) - AEPS) .OR. (ABS(angle) > ABS(angle0)) ) THEN
           IF(Direction(i) * angle < 0.0) THEN
             mindist = -dist
           ELSE
             mindist = dist
           END IF
           angle0 = angle
         END IF
           
       END IF

      END DO
!------------------------------------------------------------------------------
    END FUNCTION ComputeDistanceWithDirection
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
 END SUBROUTINE LevelSetDistance
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
!>  Compute the volume and area in 3D or area and line intgral in 2D over the 
!>  levelset function. This is better done within a dedicated solver since 
!>  it is crucial for the accuracy that the Heaviside and Delta function are
!>  computed at Gaussian integration points.
!> \ingroup Solvers
!------------------------------------------------------------------------------
   SUBROUTINE LevelSetIntegrate( Model,Solver,Timestep,TransientSimulation )
!------------------------------------------------------------------------------
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
     TYPE(Nodes_t)   :: ElementNodes
     TYPE(Element_t),POINTER :: CurrentElement
     TYPE(ValueList_t), POINTER :: Material
     TYPE(Variable_t), POINTER :: SurfSol

     INTEGER :: i,j,k,n,t,istat,bf_id
     INTEGER, POINTER :: NodeIndexes(:), SurfPerm(:)
     LOGICAL :: Visited = .FALSE., GotIt
     REAL(KIND=dp), POINTER :: Surface(:), NodalSurf(:)
     INTEGER :: body_id, dim
     REAL(KIND=dp) :: TotVolume, TotArea, Alpha, Relax, InitVolume, dSurface, dt, &
         Moment(3)
     CHARACTER(LEN=MAX_NAME_LEN) :: LevelSetVariableName
     TYPE(ValueList_t), POINTER :: Params

     SAVE ElementNodes, Visited, NodalSurf, InitVolume

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------

     Params => Solver % Values

     ! The variable that should be renormalized
     LevelSetVariableName = ListGetString(Params,'Level Set Variable',GotIt) 
     IF(GotIt) THEN
       SurfSol => VariableGet( Solver % Mesh % Variables, TRIM(LevelSetVariableName) )
     ELSE  
       SurfSol => VariableGet( Solver % Mesh % Variables, 'Surface' )
     END IF
     IF(ASSOCIATED(SurfSol)) THEN
       Surface  => SurfSol % Values
       SurfPerm => SurfSol % Perm
     ELSE
       CALL Warn('LevelSetIntegrate','Surface variable does not exist')
     END IF

     IF ( ALL( SurfPerm == 0) ) THEN
       CALL Warn('LevelSetIntegrate','Nothing to compute')
       RETURN
     END IF
 
     dim = CoordinateSystemDimension()
 
!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. Visited ) THEN
       N = Solver % Mesh % MaxElementNodes

       ALLOCATE( NodalSurf( N ), &
           ElementNodes % x( N ),   &
           ElementNodes % y( N ),   &
           ElementNodes % z( N ),   &
           STAT=istat )
 
       IF ( istat /= 0 ) THEN
         CALL Fatal( 'LevelSetIntegrate', 'Memory allocation error.' )
       END IF
     END IF
!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------

     TotVolume = 0.0d0
     TotArea = 0.0d0
     Moment = 0.0d0
     
     Alpha = ListGetConstReal(Model % Simulation,'Levelset Bandwidth',GotIt) 
     IF(.NOT. GotIt) Alpha = ListGetConstReal(Params,'Levelset Bandwidth')      
       
     CALL Info( 'LevelSetIntegrate','-------------------------------------', Level=4 )
     CALL Info( 'LevelSetIntegrate', 'Integrating over levelset function', Level=4 )
     CALL Info( 'LevelSetIntegrate','-------------------------------------', Level=4 )

     DO t=1,Solver % Mesh % NumberOfBulkElements
       
       CurrentElement => Solver % Mesh % Elements(t)
       n = CurrentElement % TYPE % NumberOfNodes
       NodeIndexes => CurrentElement % NodeIndexes
       IF( ANY(SurfPerm(NodeIndexes) == 0)) CYCLE
       
       Model % CurrentElement => CurrentElement
       body_id = CurrentElement % Bodyid    
       k = ListGetInteger( Model % Bodies( body_id ) % Values, 'Material' )
       Material => Model % Materials(k) % Values
       
!-----------------------------------------------------------------------------
!        Get element nodal coordinates
!------------------------------------------------------------------------------
       ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
       ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
       ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)
         
       NodalSurf = Surface( SurfPerm(NodeIndexes) )

       CALL HeavisideIntegrate( NodalSurf, CurrentElement, n, ElementNodes, &
           Alpha, TotVolume, TotArea, Moment)      
     END DO
     
     Moment = Moment / TotVolume
          
     IF(dim == 3) THEN
       WRITE(Message,'(a,ES12.3)') 'Center 3',Moment(3)
       CALL Info( 'LevelSetIntegrate',Message, Level=4 )

       WRITE(Message,'(a,ES12.3)') 'Inside Volume',TotVolume
       CALL Info( 'LevelSetIntegrate',Message, Level=4 )
     
       WRITE(Message,'(a,ES12.3)') 'Interface Area',TotArea
       CALL Info( 'LevelSetIntegrate',Message, Level=4 )
       
       CALL ListAddConstReal(Model % Simulation,'res: LevelSet Center 3',Moment(3))
       CALL ListAddConstReal(Model % Simulation,'res: LevelSet Volume',TotVolume)
       CALL ListAddConstReal(Model % Simulation,'res: LevelSet Area',TotArea)
     ELSE       
       WRITE(Message,'(a,ES12.3)') 'Inside Area',TotVolume
       CALL Info( 'LevelSetIntegrate',Message, Level=4 )
     
       WRITE(Message,'(a,ES12.3)') 'Interface length',TotArea
       CALL Info( 'LevelSetIntegrate',Message, Level=4 )
       
       CALL ListAddConstReal(Model % Simulation,'res: LevelSet Area',TotVolume)
       CALL ListAddConstReal(Model % Simulation,'res: LevelSet Length',TotArea)       
     END IF

     WRITE(Message,'(a,ES12.3)') 'Center 2',Moment(2)
     CALL Info( 'LevelSetIntegrate',Message, Level=4 )

     WRITE(Message,'(a,ES12.3)') 'Center 1',Moment(1)
     CALL Info( 'LevelSetIntegrate',Message, Level=4 )
     
     CALL ListAddConstReal(Model % Simulation,'res: LevelSet Center 2',Moment(2))
     CALL ListAddConstReal(Model % Simulation,'res: LevelSet Center 1',Moment(1))

     IF ( ListGetLogical(Params,'Conserve Volume',GotIt) ) THEN
       IF(.NOT. Visited) THEN
         InitVolume = ListGetConstReal(Params,'Initial Volume',GotIt)
         IF(.NOT. GotIt) InitVolume = TotVolume
       END IF

       Relax = ListGetConstReal(Params,'Conserve Volume Relaxation',GotIt)
       IF(.NOT. GotIt) Relax = 1.0d0
      
       dSurface = Relax * (InitVolume - TotVolume) / TotArea
       IF(ABS(dSurface) > AEPS) THEN
          Surface = Surface + dSurface
       END IF

       WRITE(Message,'(a,ES12.3)') 'Levelset correction',dSurface
       CALL Info( 'LevelSetIntegrate',Message, Level=4 )
       CALL ListAddConstReal(Model % Simulation,'res: Levelset Correction',dSurface)
     END IF

     Visited = .TRUE.


!------------------------------------------------------------------------------

 CONTAINS

!------------------------------------------------------------------------------
   SUBROUTINE HeavisideIntegrate( Surf,Element,n,Nodes,Alpha,Volume,Area,Moment)
!------------------------------------------------------------------------------

     REAL(KIND=dp), DIMENSION(:)   :: Surf
     INTEGER :: n
     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t), POINTER :: Element
     REAL(KIND=dp) :: Alpha, Volume, Area, Moment(3)

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: Basis(n)
     REAL(KIND=dp) :: dBasisdx(n,3),detJ
     REAL(KIND=dp) :: Val,Grad(3),Velo(3),NormalVelo,GradAbs,Heavi,Delta
     INTEGER :: t,N_Integ, CoordinateSystem
     REAL(KIND=dp) :: s,u,v,w,x,y,z
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
     LOGICAL :: stat

!------------------------------------------------------------------------------

     dim = CoordinateSystemDimension()
     CoordinateSystem = CurrentCoordinateSystem()

!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IntegStuff = GaussPoints( element )

     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

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
       stat = ElementInfo( Element,Nodes,u,v,w,detJ, Basis,dBasisdx)

       s = detJ * S_Integ(t)

       x = SUM( Nodes % x(1:n) * Basis(1:n) )
       y = SUM( Nodes % y(1:n) * Basis(1:n) )
       IF(dim == 3) z = SUM( Nodes % z(1:n) * Basis(1:n) )

!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
       IF ( CoordinateSystem == AxisSymmetric .OR. CoordinateSystem == CylindricSymmetric ) THEN
         s = s * x * 2.0d0 * PI
       END IF

!------------------------------------------------------------------------------
       Val = SUM( Basis(1:n) * Surf(1:n) )

       IF( Val < -Alpha) THEN
         Heavi = 0.0d0
       ELSE IF(Val > Alpha) THEN
         Heavi = 1.0d0
       ELSE
         Heavi = (1.0d0 + SIN( (Val/Alpha) * (PI/2) ) ) / 2.0d0
         Delta = (1.0d0 + COS( (Val/Alpha) * PI ) ) / (2.0d0 * Alpha)
         
         DO i=1,dim
           Grad(i) = SUM( dBasisdx(1:n,i) * Surf(1:n) )
         END DO
         GradAbs = SQRT( SUM( Grad(1:dim) * Grad(1:dim) ) )          
         
         Area = Area + s * Delta * GradAbs
       END IF

       Volume = Volume + s * Heavi
       Moment(1) = Moment(1) + s * Heavi * x
       Moment(2) = Moment(2) + s * Heavi * y
       IF(dim == 3) Moment(3) = Moment(3) + s * Heavi * z       

     END DO
     
!------------------------------------------------------------------------------
   END SUBROUTINE HeavisideIntegrate
!------------------------------------------------------------------------------

 END SUBROUTINE LevelSetIntegrate
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>  Computes the curvature from the level set function. Additional diffusion may 
!>  be added in order to limit the singular curvature peaks.
!> \ingroup Solvers
!------------------------------------------------------------------------------
   SUBROUTINE LevelSetCurvature( Model,Solver,Timestep,TransientSimulation )
!------------------------------------------------------------------------------
     USE DefUtils
     USE SolverUtils
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
     TYPE(Matrix_t),POINTER  :: StiffMatrix
     TYPE(Nodes_t) :: ElementNodes, ParentNodes
     TYPE(Element_t),POINTER :: Element, Parent
     TYPE(Variable_t), POINTER :: SurfSol

     INTEGER :: i,j,k,n,pn,t,istat,bf_id,CoordinateSystem
     REAL(KIND=dp) :: Norm, Coeff, Diff, Alpha, Val, Delta
     INTEGER, POINTER :: NodeIndexes(:)
     INTEGER, POINTER :: CurvPerm(:), SurfPerm(:)
     REAL(KIND=dp), POINTER :: Curvature(:),ForceVector(:), Curv(:),Surface(:) 
     REAL(KIND=dp), ALLOCATABLE :: LocalStiffMatrix(:,:),LocalForce(:),Surf(:)
     REAL(KIND=dp) :: at,totat,st,totst,CPUTime
     CHARACTER(LEN=MAX_NAME_LEN) :: LevelSetVariableName
     LOGICAL :: GotIt, Stat, AllocationsDone = .FALSE.

     SAVE LocalStiffMatrix,LocalForce, ElementNodes,ParentNodes,Surf,AllocationsDone

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
     IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

     CoordinateSystem = CurrentCoordinateSystem()

     Curvature => Solver % Variable % Values
     CurvPerm => Solver % Variable % Perm
     IF ( SIZE( Curvature ) == 0 ) RETURN

     LevelSetVariableName = ListGetString(Solver % Values,'LevelSet Variable',GotIt) 
     IF(GotIt) THEN
       SurfSol => VariableGet( Solver % Mesh % Variables, TRIM(LevelSetVariableName) )
     ELSE  
       SurfSol => VariableGet( Solver % Mesh % Variables, 'Surface' )
     END IF
     Surface => Surfsol % Values
     SurfPerm => SurfSol % Perm

     StiffMatrix => Solver % Matrix
     ForceVector => StiffMatrix % RHS
     Norm = Solver % Variable % Norm

     Diff = ListGetConstReal(Solver % Values,'Curvature Diffusion',GotIt)

!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone ) THEN
       N = Solver % Mesh % MaxElementNodes

       ALLOCATE( ElementNodes % x( N ),   &
           ElementNodes % y( N ),   &
           ElementNodes % z( N ),   &
           ParentNodes % x( N ),   &
           ParentNodes % y( N ),   &
           ParentNodes % z( N ),   &
           LocalForce( N ),       &
           LocalStiffMatrix( N, N ), &
           Surf( N ), &
           STAT=istat )
 
       IF ( istat /= 0 ) THEN
         CALL Fatal( 'CurvatureSolve', 'Memory allocation error.' )
       END IF

       AllocationsDone = .TRUE.
     END IF
!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------

     totat = 0.0d0
     totst = 0.0d0

     at = CPUTime()
     
     CALL Info( 'LevelSetCurvature','-------------------------------------', Level=4 )
     CALL Info( 'LevelSetCurvature','Solving Level set curvature', Level=4 )
     CALL Info( 'LevelSetCurvature','-------------------------------------', Level=4 )

     CALL DefaultInitialize()

!------------------------------------------------------------------------------
!      Do the assembly for bulk elements
!------------------------------------------------------------------------------

     DO t=1,Solver % NumberOfActiveElements
       
       Element => Solver % Mesh % Elements(Solver % ActiveElements(t))
       n = Element % TYPE % NumberOfNodes
       NodeIndexes => Element % NodeIndexes
       Model % CurrentElement => Element
 
       ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
       ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
       ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)

       Surf = Surface( SurfPerm(NodeIndexes) )

       CALL LocalMatrix( LocalStiffMatrix, LocalForce, &
           Surf, Element, n, ElementNodes )
         
       CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce )

     END DO     
     CALL DefaultFinishBulkAssembly()

!------------------------------------------------------------------------------
!      Do the assembly for boundary elements
!------------------------------------------------------------------------------

     DO t=Solver % Mesh % NumberOfBulkElements + 1, &
         Solver % Mesh % NumberOfBulkElements + &
         Solver % Mesh % NumberOfBoundaryElements
      
       Element => Solver % Mesh % Elements(t)
!------------------------------------------------------------------------------
       DO i=1,Model % NumberOfBCs
         IF ( Element % BoundaryInfo % Constraint == Model % BCs(i) % Tag ) THEN

           n = Element % TYPE % NumberOfNodes
           NodeIndexes => Element % NodeIndexes
           Model % CurrentElement => Element
          
           IF ( ANY( CurvPerm(NodeIndexes) <= 0 ) ) CYCLE         

           ! Check that the dimension of element is suitable for fluxes
           IF( .NOT. PossibleFluxElement(Element) ) CYCLE
           
           IF ( .NOT. ListGetLogical(Model % BCs(i) % Values, &
               'Levelset Curvature BC',gotIt) ) CYCLE

           ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
           ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
           ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)
           
           Parent => Element % BoundaryInfo % Left
          
           stat = ASSOCIATED( Parent )
           IF ( stat ) stat = stat .AND. ALL(CurvPerm(Parent % NodeIndexes) > 0)
          
           IF ( .NOT. stat ) THEN
             Parent => Element % BoundaryInfo % Right
             
            stat = ASSOCIATED( Parent )
            IF ( stat ) stat = ALL(CurvPerm(Parent % NodeIndexes(1:pn)) > 0)
            
            IF ( .NOT. stat )  THEN
              CALL Warn( 'LevelSetCurvature', &
                  'No curvature solution available for specified boundary' )
              CYCLE
            END IF
          END IF
          
          pn = Parent % TYPE % NumberOfNodes           
          ParentNodes % x(1:pn) = Solver % Mesh % Nodes % x(Parent % NodeIndexes)
          ParentNodes % y(1:pn) = Solver % Mesh % Nodes % y(Parent % NodeIndexes)
          ParentNodes % z(1:pn) = Solver % Mesh % Nodes % z(Parent % NodeIndexes)

          Surf(1:pn) = Surface(SurfPerm(Parent % NodeIndexes))
          
!------------------------------------------------------------------------------
!             Get element matrix and rhs due to boundary conditions ...
!------------------------------------------------------------------------------
          
          CALL LocalBoundary( LocalStiffMatrix, LocalForce,  &
              Surf, Element, Parent, n, pn, ElementNodes, ParentNodes )
!------------------------------------------------------------------------------
!             Update global matrices from local matrices
!------------------------------------------------------------------------------
          CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce )
!------------------------------------------------------------------------------
           END IF   
        END DO    
      END DO   
!------------------------------------------------------------------------------

     CALL DefaultFinishAssembly()
     
     at = CPUTime() - at
     WRITE(Message,'(a,F8.2)') 'Assembly done in time (s):',at
     CALL Info( 'LevelSetCurvature',Message, Level=4 )

!------------------------------------------------------------------------------
!    Solve the system and we are done.
!------------------------------------------------------------------------------
     st = CPUTime()
     CALL SolveSystem( StiffMatrix, ParMatrix, ForceVector, &
         Curvature, Norm, 1, Solver )
     
     st = CPUTIme()-st
     WRITE(Message,'(a,F8.2)') 'Solution done in time (s):',st
     CALL Info( 'LevelSetCurvature',Message, Level=4 )

!------------------------------------------------------------------------------

     Coeff = ListGetConstReal(Solver % Values,'Curvature Coefficient',GotIt) 
     IF(GotIt) THEN
       Curvature = Coeff * Curvature
     END IF

     Alpha = ListGetConstReal(Model % Simulation,'Levelset Bandwidth',GotIt) 
     IF(.NOT. GotIt) Alpha = ListGetConstReal(Solver % Values,'Levelset Bandwidth',GotIt)      
     IF(GotIt) THEN
       DO i=1,SIZE(CurvPerm) 
         j = CurvPerm(i)
         k = SurfPerm(i)
         IF(j == 0 .OR. k == 0) CYCLE
         Val = Surface(k)

         IF( Val < -Alpha) THEN
           Delta = 0.0d0
         ELSE IF(Val > Alpha) THEN
           Delta = 0.0d0
         ELSE
           Delta = (1.0d0 + COS( (Val/Alpha) * PI ) ) / (2.0d0 * Alpha)
         END IF
         
         Curvature(j) = Delta * Curvature(j)
       END DO
     END IF


CONTAINS

   SUBROUTINE LocalMatrix( StiffMatrix,ForceVector, Surf, &
       Element, n, Nodes )
!------------------------------------------------------------------------------

     REAL(KIND=dp), DIMENSION(:)   :: ForceVector,Surf
     REAL(KIND=dp), DIMENSION(:,:) :: StiffMatrix
     INTEGER :: n
     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t), POINTER :: Element

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
!
     REAL(KIND=dp) :: Basis(n)
     REAL(KIND=dp) :: dBasisdx(n,3),detJ
     REAL(KIND=dp) :: Grad(3),GradAbs,A,B,xpos
     INTEGER :: i,j,k,c,p,q,t,dim,N_Integ
     REAL(KIND=dp) :: s,u,v,w
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     LOGICAL :: stat

!------------------------------------------------------------------------------

     dim = CoordinateSystemDimension()
     Grad = 0.0d0
     ForceVector = 0.0D0
     StiffMatrix = 0.0D0

!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IntegStuff = GaussPoints( element )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

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
       stat = ElementInfo( Element,Nodes,u,v,w,detJ,Basis,dBasisdx)

       s = detJ * S_Integ(t)

!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
       IF ( Coordinates == AxisSymmetric .OR. Coordinates == CylindricSymmetric ) THEN
         xpos = SUM( Nodes % x(1:n) * Basis(1:n) )
         s = xpos * s
       END IF

!------------------------------------------------------------------------------

       DO i=1,dim
         Grad(i) = SUM( dBasisdx(1:n,i) * Surf(1:n) )
       END DO
       GradAbs = SQRT( SUM( Grad * Grad ) )
       IF ( GradAbs > 10*AEPS ) THEN
         Grad = Grad / GradAbs
       END IF

!------------------------------------------------------------------------------
!      Loop over basis functions of both unknowns and weights
!------------------------------------------------------------------------------
       DO p=1,n

         DO q=1,n
           A = Basis(q) * Basis(p)
           A = A + Diff * SUM( dBasisdx(q,1:dim) * dBasisdx(p,1:dim) )

           StiffMatrix(p,q) = StiffMatrix(p,q) + s * A
         END DO

         B = - SUM( dBasisdx(p, 1:dim) * Grad(1:dim))  
         ForceVector(p) = ForceVector(p) + s * B        

       END DO
     END DO

!------------------------------------------------------------------------------
   END SUBROUTINE LocalMatrix


!------------------------------------------------------------------------------
   SUBROUTINE LocalBoundary( BoundaryMatrix, BoundaryVector, &
        Surf, Element, Parent, n, pn, ElementNodes, ParentNodes )
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: BoundaryMatrix(:,:), BoundaryVector(:), Surf(:)
     TYPE(Nodes_t)   :: ElementNodes, ParentNodes
     TYPE(Element_t), POINTER :: Element, Parent
     INTEGER :: n, pn

     REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), Normal(3), Grad(3), GradAbs, detJ, A
     REAL(KIND=dp) :: ParentBasis(pn), ParentdBasisdx(pn,3)
     REAL(KIND=dp) :: u,v,w,s,x(n),y(n),z(n),xpos, Force
     REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)

     INTEGER :: i,j,k,t,p,q,N_Integ,dim

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: stat
!------------------------------------------------------------------------------

     BoundaryVector = 0.0d0
     BoundaryMatrix = 0.0d0
     dim = CoordinateSystemDimension()
     Grad = 0.0d0

!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IntegStuff = GaussPoints( Element )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------
     DO t=1,N_Integ
       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)
!------------------------------------------------------------------------------
!     Basis function values & derivates at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,ElementNodes,u,v,w,detJ,Basis,dBasisdx )
       s = S_Integ(t) * detJ

!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
      IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
        xpos = SUM( ElementNodes % x(1:n)*Basis(1:n) )
        s = xpos * s
      END IF

      Normal = Normalvector( Element, ElementNodes, u, v, .TRUE. )
      
      !------------------------------------------------------------------------------
      ! Need parent element basis etc., for computing normal derivatives on boundary.
      !------------------------------------------------------------------------------
      
      DO i = 1,n
        DO j = 1,pn
          IF ( Element % NodeIndexes(i) == Parent % NodeIndexes(j) ) THEN
            x(i) = Parent % TYPE % NodeU(j)
            y(i) = Parent % TYPE % NodeV(j)
            z(i) = Parent % TYPE % NodeW(j)
            EXIT
          END IF
        END DO
      END DO
      
      u = SUM( Basis(1:n) * x(1:n) )
      v = SUM( Basis(1:n) * y(1:n) )
      w = SUM( Basis(1:n) * z(1:n) )
      
      stat = ElementInfo( Parent, ParentNodes, u, v, w, detJ, ParentBasis, ParentdBasisdx )
          
      DO j = 1,DIM
        Grad(j) = SUM( ParentdBasisdx(1:pn,j) * Surf(1:pn) )
      END DO
      GradAbs = SQRT( SUM(Grad(1:dim) * Grad(1:dim)) )

      IF ( GradAbs > 10*AEPS ) THEN         
        Grad = Grad / GradAbs
      END IF

!------------------------------------------------------------------------------
      
      DO p=1,n
        DO q=1,n
          A = Diff * SUM( Normal(1:dim) * dBasisdx(q,1:dim)) * Basis(p)
          BoundaryMatrix(p,q) = BoundaryMatrix(p,q) + s * A
        END DO
      END DO

      Force = SUM( Normal(1:dim) * Grad(1:dim) )
      DO q=1,N
        BoundaryVector(q) = BoundaryVector(q) + s * Basis(q) * Force
      END DO

    END DO

  END SUBROUTINE LocalBoundary
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE LevelSetCurvature
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Determines a timestep based on the maximum local Courant number. 
!> \ingroup UDF
!------------------------------------------------------------------------------
FUNCTION LevelSetTimestep( Model ) RESULT( dt )
  USE Types
  USE Lists
  USE Integration
  USE ElementDescription

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt

  TYPE(Solver_t), POINTER :: Solver
  TYPE(Variable_t), POINTER ::  TimeVariable, SurfSol
  REAL(KIND=dp), POINTER :: Surface(:), Surf(:)
  INTEGER, POINTER :: SurfPerm(:), NodeIndexes(:)
  TYPE(Nodes_t)   :: ElementNodes
  TYPE(Element_t),POINTER :: CurrentElement
  TYPE(ValueList_t), POINTER :: Material
  REAL(KIND=dp), POINTER :: TimestepSizes(:,:)
 
  INTEGER :: i,j,k,n,t,elem,N_Integ,body_id,mat_id, dim, TimeIntervals
  REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:), NodalVelo(:,:)
  REAL(KIND=dp) :: Val,Grad(3),Velo(3),NormalVelo,AbsVelo,GradAbs,detJ,&
      MaxNormVelo,MaxAbsVelo
  REAL(KIND=dp) :: s,u,v,w,dt0,prevdt,cumtime,dsmax
  TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
  REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
  LOGICAL :: stat, GotIt, AllocationsDone = .FALSE.
  CHARACTER(LEN=MAX_NAME_LEN) :: LevelSetVariableName

  SAVE AllocationsDone, Basis, dBasisdx, NodalVelo, Surf, ElementNodes, prevdt


  TimestepSizes => ListGetConstRealArray( CurrentModel % Simulation,&
      'Timestep Sizes', GotIt )
  TimeIntervals = SIZE(TimestepSizes)
 
  IF(TimeIntervals > 1) THEN
    CALL Warn('LevelSetTimestep','Implemented only for one Time Interval')
  END IF

  dt0 = TimestepSizes(1,1)

  TimeVariable => VariableGet(Model % Variables, 'Time')
  cumtime = TimeVariable % Values(1)

  dim = CoordinateSystemDimension()

  ! The variable that should be reinitialized
  LevelSetVariableName = 'Surface'
  SurfSol => VariableGet( Model % Variables, TRIM(LevelSetVariableName) )
  IF(ASSOCIATED(SurfSol)) THEN
    Surface  => SurfSol % Values
    SurfPerm => SurfSol % Perm
  ELSE
    CALL Warn('LevelSetTimeStep','SurfSol does not exist: '//TRIM(LevelSetVariableName))
    RETURN
  END IF

  Solver => SurfSol % Solver
  DsMax = ListGetConstReal(Model % Simulation,'LevelSet Courant Number',GotIt)
  IF(.NOT. GotIt) DsMax = 1.0d0

  IF(.NOT. AllocationsDone) THEN
    N = Solver % Mesh % MaxElementNodes
    ALLOCATE( Basis(n), dBasisdx(n,3), NodalVelo(3,n), Surf(n), &
        ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n) )
    prevdt = dt0
    AllocationsDone = .TRUE. 
  END IF

  MaxNormVelo = 0.0d0
  MaxAbsVelo = 0.0d0

  DO elem=1,Solver % Mesh % NumberOfBulkElements
    
    CurrentElement => Solver % Mesh % Elements(elem)
    n = CurrentElement % TYPE % NumberOfNodes
    NodeIndexes => CurrentElement % NodeIndexes
    
    IF(ANY(SurfPerm(NodeIndexes) == 0)) CYCLE
    
    Surf(1:n) = Surface( SurfPerm(NodeIndexes) )
    IF(ALL(Surf(1:n) < 0.0d0) .OR. ALL(Surf(1:n) > 0.0d0) ) CYCLE

    ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
    ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
    ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)
    
    Model % CurrentElement => CurrentElement
    body_id = CurrentElement % Bodyid    
    mat_id = ListGetInteger( Model % Bodies( body_id ) % Values, 'Material' )
    Material => Model % Materials(mat_id) % Values

!------------------------------------------------------------------------------
!         Computed or given velocity field
!------------------------------------------------------------------------------
       
    NodalVelo(1,1:n) = ListGetReal( Material,'Levelset Velocity 1',n,NodeIndexes,GotIt)
    NodalVelo(2,1:n) = ListGetReal( Material,'Levelset Velocity 2',n,NodeIndexes,GotIt)
    IF(dim == 3) NodalVelo(3,1:n) = ListGetReal( Material,'Levelset Velocity 3',n,NodeIndexes,GotIt)

!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
    IntegStuff = GaussPoints( CurrentElement )
    
    U_Integ => IntegStuff % u
    V_Integ => IntegStuff % v
    W_Integ => IntegStuff % w
    S_Integ => IntegStuff % s
    N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!    Maximum at any integration point
!------------------------------------------------------------------------------
    
    DO t=1,N_Integ

      u = U_Integ(t)
      v = V_Integ(t)
      w = W_Integ(t)
         
!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------

      stat = ElementInfo( CurrentElement,ElementNodes,u,v,w,detJ, Basis,dBasisdx)
      
      DO i=1,dim
        Grad(i) = SUM( dBasisdx(1:n,i) * Surf(1:n) )
        Velo(i) = SUM( Basis(1:n) * NodalVelo(i,1:n) )
      END DO
      
      GradAbs = SQRT( SUM( Grad(1:dim) * Grad(1:dim) ) )         
      NormalVelo = SUM( Grad(1:dim) * Velo(1:dim) ) / GradAbs
      NormalVelo = NormalVelo / SQRT(detJ)
      MaxNormVelo = MAX(MaxNormVelo, ABS(NormalVelo))
      
      AbsVelo = SQRT(SUM (Velo(1:dim) * Velo(1:dim)) )
      AbsVelo = AbsVelo / SQRT(detJ)
      MaxAbsVelo = MAX(MaxAbsVelo, ABS(AbsVelo))
    END DO
  END DO

  dt = dt0
  IF( ListGetLogical(Model % Simulation,'LevelSet Timestep Directional',GotIt) ) THEN
    IF( MaxNormVelo * dt0 > dsMax) dt = dsMax / MaxNormVelo
  ELSE
    IF( MaxAbsVelo * dt0 > dsMax) dt = dsMax / MaxAbsVelo
  END IF
  
  IF( dt < dt0) THEN
    IF(dt > prevdt) dt = 0.5d0 * (dt + prevdt)
  END IF
  prevdt = dt
  
  WRITE(Message,'(a,ES12.3)') 'Levelset timestep',dt
  CALL Info( 'LevelSetTimestep',Message, Level=4 )
  CALL ListAddConstReal(Model % Simulation,'res: Levelset timestep',dt)

END FUNCTION LevelSetTimestep

