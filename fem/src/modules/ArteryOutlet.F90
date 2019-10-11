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
! *  Solve for the 1D characteristics equation arising from newtonian flow in an 
! *  elastic tube. The subroutine may be used as the outlet for blood flow simulation.
! *
! ******************************************************************************
! *
! *  Authors: Esko J�rvinen, Mikko Lyly, Peter Råback
! *  Email:   Esko.Jarvinen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 20 Nov 2001
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> Initialization of the primary solver. 
!> If requested create an internal mesh.
!------------------------------------------------------------------------------
SUBROUTINE OutletCompute_Init( Model,Solver,dt,TransientSimulation )

  USE MeshUtils
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Mesh_t), POINTER :: Mesh, PMesh


  Params => GetSolverParams()
  IF( GetLogical( Params,'1D Mesh Create') ) THEN
    CALL Info('OutletCompute_Init','Creating internal 1D mesh')

    Mesh => CreateLineMesh( Params )
    Solver % Mesh => Mesh 

    PMesh => Model % Meshes
    IF( ASSOCIATED( PMesh ) ) THEN
      DO WHILE ( ASSOCIATED( PMesh % Next ) ) 
        Pmesh => PMesh % Next
      END DO
      Pmesh % Next => Mesh
    END IF
  END IF
    
END SUBROUTINE OutletCompute_Init

 
!------------------------------------------------------------------------------
!>  Solve for the 1D characteristics equation arising from newtonian flow in an 
!>  elastic tube. The subroutine may be used as the outlet for blood flow simulation.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE OutletCompute( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE Types
  USE Lists 
  USE Integration
  USE ElementDescription
  USE SolverUtils
  USE MeshUtils
  USE DefUtils
  USE MaterialModels
  USE ElementUtils
  USE ModelDescription

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t):: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
 
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
 
  TYPE(Matrix_t),POINTER  :: StiffMatrix
  TYPE(Nodes_t)   :: ElementNodes
  TYPE(Element_t),POINTER :: CurrentElement, Element
  TYPE(ValueList_t), POINTER :: Material
  TYPE(Variable_t), POINTER :: LVar, FlowSol
  TYPE(Mesh_t), POINTER :: Mesh1D, Mesh3D
  TYPE(Solver_t), POINTER :: FlowSolver

  INTEGER :: t, k,n,m,ie,bf_id,mat_id,prev_mat_id,istat,LocalNodes,i,j, nonliniter, l, &
      bc, joinnode, Connections, SolidConnections, fsstep, fsstepmax, NonlinearIter, &
      mindisti, jf, js, NoActive
  
  INTEGER, POINTER :: NodeIndexes(:), WPerm(:)
  INTEGER, ALLOCATABLE :: LumpedBoundaries(:), SolidEndBoundaries(:)
  
  REAL(KIND=dp) :: Norm, PrevNorm, FlowOut2d
  
  LOGICAL :: GotIt, GotIt2, AllocationsDone = .FALSE. , FirstTime = .TRUE.
  LOGICAL :: CenterPointsComputed = .FALSE., InternalMesh

  CHARACTER(LEN=MAX_NAME_LEN) :: Method, Name
  
  REAL(KIND=dp), POINTER :: Wnodal(:), ForceVector(:), Lnodal(:), Anodal(:), &
      Pnodal(:), Qnodal(:), FlowSolution(:)
  
  REAL(KIND=dp), ALLOCATABLE :: LocalStiffMatrix(:,:),LocalMassMatrix(:,:), &
      Lelem(:), LocalForce(:), Cnodal(:), Unodal(:), Wprev(:), Lprev(:), &
      Aprev(:), FluidicForces(:), FluidicAreas(:), FluidicFluxes(:), SolidEndAreas(:)
  
  REAL(KIND=dp) :: NonlinearTol, Density, PoissonRatio,&
      WallYoungs, ArteryRadius, ArteryWallThickness, BetaCoeff, Aref, W20, &
      WdiffSum, WnodalSum, WprevnodalSum, Werror, fsAlpha, fsTheta, fsdTheta, fsBeta, &
      pres2dout, area2dout, psi, RadiusIn1d, AreaIn1d, BetaCoeff2, &
      dR, Wbcnode, A2, Q2, dist, mindist  
  
  SAVE LocalStiffMatrix,LocalMassMatrix,&
      Lelem, LocalForce,ElementNodes, AllocationsDone, &
      Cnodal, Unodal, Wprev, Lprev, Aprev,&
      Lnodal, Anodal, Qnodal, Pnodal, Wnodal, &
      FluidicForces, FluidicAreas, FluidicFluxes, & 
      LumpedBoundaries, SolidEndBoundaries, SolidEndAreas, &
      CenterPointsComputed

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------

  CALL Info('OutletCompute','Starting')
  
  Mesh1D => Solver % Mesh

  Wnodal => Solver % Variable % Values
  LocalNodes = SIZE( Wnodal )

  IF(LocalNodes == 0) THEN
    CALL Fatal('OutletCompute','No active variables')
  END IF

  WPerm => Solver % Variable % Perm
  StiffMatrix => Solver % Matrix
  IF(.NOT. ASSOCIATED(StiffMatrix)) THEN
    CALL Fatal('OutletCompute','StiffMatrix does not exist')
  END IF
  ForceVector => StiffMatrix % RHS
  Norm = Solver % Variable % Norm

  InternalMesh = GetLogical( Solver % Values,'1D Mesh Create')

!------------------------------------------------------------------------------
! Find the solver and mesh used for the flow simulation
!------------------------------------------------------------------------------
  Name = ListGetString( Solver % Values,'Flow Solution Equation',GotIt)
  IF(.NOT. GotIt) Name = 'navier-stokes'

  DO i=1,Model % NumberOfSolvers
    FlowSolver => Model % Solvers(i)
    IF ( ListGetString( FlowSolver % Values, 'Equation' ) == TRIM(Name) ) EXIT
  END DO    
  Mesh3D => FlowSolver % Mesh

!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
  IF ( .NOT. AllocationsDone ) THEN
    N = Model % MaxElementNodes
    M = Model % NumberOfBCs
    
    ALLOCATE( ElementNodes % x( N ),   &
        ElementNodes % y( N ),   &
        ElementNodes % z( N ),   &
        LocalForce( N ),         &
        LocalStiffMatrix( N,N ), &
        LocalMassMatrix( N,N ),  &
        Lelem( N ),              &
        LumpedBoundaries( M ),   &
        SolidEndBoundaries( M ),   &
        FluidicAreas( M ),       &
        SolidEndAreas( M ),       &
        FluidicForces( M ),      &
        FluidicFluxes( M ),      &
        STAT=istat )
    
    IF ( istat /= 0 ) THEN
      CALL Fatal('ArteryOnedim','Memory allocation error')
    END IF
    
    LVar => VariableGet( Mesh1D % Variables, 'Lnodal' )
    IF ( ASSOCIATED( LVar ) ) THEN
      Lnodal => LVar % Values
    ELSE
      ALLOCATE( Lnodal(LocalNodes))
    END IF
    
    LVar => VariableGet( Mesh1D % Variables, 'Anodal' )
    IF ( ASSOCIATED( LVar ) ) THEN
      Anodal => LVar % Values
    ELSE
      ALLOCATE( Anodal(LocalNodes))
    END IF
    
    LVar => VariableGet( Mesh1D % Variables, 'Pnodal' )
    IF ( ASSOCIATED( LVar ) ) THEN
      Pnodal => LVar % Values
    ELSE
      ALLOCATE( Pnodal(LocalNodes))
    END IF
    
    LVar => VariableGet( Mesh1D % Variables, 'Qnodal' )
    IF ( ASSOCIATED( LVar ) ) THEN
      Qnodal => LVar % Values
    ELSE
      ALLOCATE( Qnodal(LocalNodes))
    END IF
    
    ALLOCATE( Cnodal( LocalNodes ), & 
        Unodal( LocalNodes ), & 
        Wprev( LocalNodes ),  &
        Lprev( LocalNodes ),  &
        Aprev( LocalNOdes ) )
    
    Cnodal = 0
    Unodal = 0
    Wprev = 0
    Lprev = 0
    Aprev = 0
    
    AllocationsDone = .TRUE.
  END IF
  
  !---------------------------------------------------------------------------------------
  ! Compute the lumped values for area, flux and pressure on the corresponding boundaries
  !---------------------------------------------------------------------------------------

  ! Check which boundaries should be lumped and mark them with the corresponding 1D dof
  ! This is the old obsolete method!
  LumpedBoundaries = 0
  SolidEndBoundaries = 0

  DO bc=1,Model % NumberOfBCs       
    
    ! Coupling with what Boundary Condition nro (j)   
    jf = ListGetInteger( Model % BCs(bc) % Values, 'Fluid Coupling With Boundary', GotIt)
    IF(.NOT. GotIt) CYCLE
    
    js = ListGetInteger( Model % BCs(bc) % Values, 'Structure Coupling With Boundary', GotIt)
    
    DO t = Mesh1D % NumberOfBulkElements + 1, &
        Mesh1D % NumberOfBulkElements + Mesh1D % NumberOfBoundaryElements
      
      Element => Mesh1D % Elements(t)
      
      ! In which Boundary Condition coupling is given 
      IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) CYCLE
      
      ! Boundary Element should be 1 node boundary
      Model % CurrentElement => Element
      n = Element % TYPE % NumberOfNodes
      IF(n /= 1) CYCLE 
      
      ! Element % NodeIndexes(1) is boundary elements node number, k is the node nro of
      ! the 1D model
      k = WPerm(Element % NodeIndexes(1))

      LumpedBoundaries(jf) = k     
      IF( js > 0 ) SolidEndBoundaries(js) = k
      
      EXIT
    END DO
  END DO
  
  Connections = COUNT( LumpedBoundaries > 0 )
  SolidConnections = COUNT( SolidEndBoundaries > 0 )

  ! If there are no connections then check the connections from 
  ! the solver section. This may hold only one connection!
  !------------------------------------------------------------
  IF( Connections == 0 ) THEN
    jf = ListGetInteger( Solver % Values, 'Fluid Coupling With Boundary', GotIt)
    js = ListGetInteger( Solver % Values, 'Structure Coupling With Boundary', GotIt)

    IF( jf == 0 ) THEN
      DO bc=1,Model % NumberOfBCs       
        IF( ListGetLogical( Model % BCs(bc) % Values, '1D Fluid Coupling', GotIt) ) THEN
          jf = bc
          EXIT
        END IF
      END DO
    END IF
    IF( js == 0 ) THEN
      DO bc=1,Model % NumberOfBCs       
        IF( ListGetLogical( Model % BCs(bc) % Values, '1D Structure Coupling', GotIt) ) THEN
          js = bc
          EXIT
        END IF
      END DO
    END IF

    IF( jf == 0 ) THEN
      CALL Fatal('OutletCompute','1D boundary not connected!')
    END IF

    mindisti = 0
    mindist = HUGE( mindist ) 
    DO t=1,Mesh1D % NumberOfNodes
      dist = Mesh1D % Nodes % x(t) + &
          Mesh1D % Nodes % y(t) + &
          Mesh1D % Nodes % z(t)
      IF( dist < mindist ) THEN
        mindisti = WPerm( t ) 
        mindist = dist
      END IF
    END DO

    IF( mindisti == 0 ) THEN
      CALL Fatal('OutletCompute','1D boundary not defined!')
    END IF

    k = WPerm( mindisti )
    LumpedBoundaries(jf) = k
    IF( js > 0 ) SolidEndBoundaries(js) = k
  END IF
  
  ! LumpedBoundaries is a vector a_j, a_j = k (1d node nro), if a_i /= 0,
  ! and j is boundary which is coupled with the lumped boundary
  !---------------------------------------------------------------------------
  
  CALL LumpedFluidicForce( LumpedBoundaries, FluidicForces, FluidicAreas,FluidicFluxes)
  
  DO bc=1,Model % NumberOfBCs      
    IF(LumpedBoundaries(bc) == 0) CYCLE
    
    WRITE(Message,'(A,I3)') 'Lumping Boundary:',bc
    CALL Info('OutletCompute',Message)
    WRITE(Message,'(A,E11.4)') 'Force:',FluidicForces(bc)
    CALL Info('OutletCompute',Message)
    WRITE(Message,'(A,E11.4)') 'Area:',FluidicAreas(bc)
    CALL Info('OutletCompute',Message)
    WRITE(Message,'(A,E11.4)') 'Flux:',FluidicFluxes(bc)
    CALL Info('OutletCompute',Message)    
  END DO


!---------------------------------------------------------------------------------------
! Compute the center point coordinates of the solid boundaries
!---------------------------------------------------------------------------------------

  IF ( .NOT. CenterPointsComputed ) THEN
    
    IF ( CurrentCoordinateSystem() == Cartesian ) THEN
      IF( SolidConnections == 0 ) THEN
        CALL Fatal('ArteryOutlet','Structure coupling should be defined!')
      END IF
            
      CALL SurfaceCenterPoints( SolidEndBoundaries, SolidConnections, SolidEndAreas)
      
      DO bc=1,Model % NumberOfBCs  
        IF(SolidEndBoundaries(bc) /= 0) THEN          
          WRITE(Message,'(A,I3)') 'Solid End Area Boundary: ',bc
          CALL Info('OutletCompute',Message)
          WRITE(Message,'(A,E11.4)') 'Solid End Area:',SolidEndAreas(bc)
          CALL Info('OutletCompute',Message)          
        END IF
      END DO      
    ELSE      
      WRITE(Message,'(A)') ' Coordinate System is not Cartesian.'
      CALL Info('OutletCompute',Message)
      WRITE(Message,'(A)') ' No center points of the outlets are computed.'
      CALL Info('OutletCompute',Message)      
    END IF
    
    CenterPointsComputed = .TRUE.
  END IF

  !------------------------------------------------------------------------------
  !    Do some additional initialization
  !------------------------------------------------------------------------------
  CALL InitializeToZero( StiffMatrix, ForceVector )
  
  NonlinearTol = ListGetConstReal( Solver % Values, &
      'Nonlinear System Convergence Tolerance',minv=0.0d0 )
  
  NonlinearIter = ListGetInteger( Solver % Values, &
      'Nonlinear System Max Iterations', minv=0 )
  
  Method = ListGetString( Solver % Values, 'Timestepping Method', GotIt )     
  prev_mat_id = -1
  
  !------------------------------------------------------------------------------
  ! Fractional step (optional) initialization 
  !------------------------------------------------------------------------------
  
  IF ( Method == 'fs') THEN 
    fsTheta = ListGetConstReal( CurrentModel % Simulation,'FS Theta', GotIt )
    
    IF(.NOT. GotIt) fsTheta = 1.0 - SQRT(2.0)/2.0
    fsdTheta = 1.0 - 2.0 * fsTheta
    fsAlpha  = fsdTheta / ( 1.0 - fsTheta )
    fsBeta   = 1.0 - fsAlpha        
    fsstepmax = 3
    
    CALL ListAddConstReal( Solver % Values, 'fsTheta', fsTheta )
    CALL ListAddConstReal( Solver % Values, 'fsdTheta', fsdTheta )
    CALL ListAddConstReal( Solver % Values, 'fsAlpha', fsAlpha )
    CALL ListAddConstReal( Solver % Values, 'fsBeta', fsBeta )
  ELSE
    fsstepmax = 1
  END IF
  
  !------------------------------------------------------------------------------
  ! Fractional step loop
  !------------------------------------------------------------------------------
  
  DO fsstep = 1,fsstepmax
    
    IF ( Method == 'fs') THEN 
      CALL InitializeTimeStep(Solver)
      Solver % Variable % PrevValues(:,1) = Solver % Variable % Values 
      CALL ListAddConstReal( Solver % Values, 'fsstep', 1.0d0*fsstep )
      
      WRITE(Message,'(A,I3)') 'Fractional Step:',fsstep
      CALL Info('OutletCompute',Message)
    END IF
    
    !------------------------------------------------------------------------------
    ! Nonlinear iteration
    !------------------------------------------------------------------------------

    prev_mat_id = -1
    
    DO nonliniter = 1,NonlinearIter
      
      CALL InitializeToZero( StiffMatrix, ForceVector )
      
      WRITE(Message,'(A,I0)') 'Nonlinear iteration: ',nonliniter
      CALL Info('OutletCompute',Message)
      
      !------------------------------------------------------------------------------
      !    Do the assembly
      !------------------------------------------------------------------------------      
      IF( InternalMesh ) THEN
        NoActive = Mesh1D % NumberOfBulkElements 
      ELSE
        NoActive = Solver % NumberOfActiveElements
      END IF


      DO ie = 1, NoActive

        IF( InternalMesh ) THEN        
          CurrentElement => Mesh1D % Elements( ie )
        ELSE
          CurrentElement => Mesh1D % Elements( Solver % ActiveElements(ie) )
        END IF

        Model % CurrentElement => CurrentElement
        n = CurrentElement % TYPE % NumberOfNodes
        
        NodeIndexes => CurrentElement % NodeIndexes
        ElementNodes % x(1:n) = Mesh1D % Nodes % x(NodeIndexes)
        ElementNodes % y(1:n) = Mesh1D % Nodes % y(NodeIndexes)
        ElementNodes % z(1:n) = Mesh1D % Nodes % z(NodeIndexes)
        
        mat_id = ListGetInteger( Model % Bodies( CurrentElement % BodyId ) % Values, 'Material')    
        
        IF(mat_id /= prev_mat_id) THEN
          Material => Model % Materials(mat_id) % Values
          
          Density = ListGetConstReal( Material, 'Density')
          WallYoungs = ListGetConstReal( Material, 'Artery Wall Youngs Modulus')
          ArteryRadius = ListGetConstReal( Material, 'Artery Radius')
          ArteryWallThickness = ListGetConstReal( Material, 'Artery Wall Thickness')
          PoissonRatio = ListGetConstReal( Material, 'Artery Poisson Ratio')
          
          Aref = Pi * ArteryRadius**2
          BetaCoeff = ( 1 / (1 - PoissonRatio**2) ) * &
              ( SQRT(pi) * ArteryWallThickness  * WallYoungs ) / ( Aref)
          W20 = -4 * SQRT( (BetaCoeff/(2.0*Density) ) * SQRT(Aref) ) 
          
          prev_mat_id = mat_id
          
        END IF
        
        IF( FirstTime ) THEN
          Anodal( Wperm(NodeIndexes) ) = Aref
          Unodal( Wperm(NodeIndexes) ) = 0.0 
          Cnodal( Wperm(NodeIndexes) ) = SQRT( (BetaCoeff/(2.0*Density))  * SQRT(Aref) )
          Lnodal( Wperm(NodeIndexes) ) = SQRT( (BetaCoeff/(2.0*Density))  * SQRT(Aref) )
          Pnodal( Wperm(NodeIndexes) ) = 0.0
        END IF
        
        Lelem(1:n) = Lnodal( Wperm(NodeIndexes) )
        
        !------------------------------------------------------------------------------
        !      Get element local matrix, and rhs vector
        !------------------------------------------------------------------------------
        CALL LocalMatrix(  LocalStiffMatrix, LocalMassMatrix, &
            LocalForce, Lelem, CurrentElement, n, ElementNodes)
        
        IF( TransientSimulation ) THEN
          CALL Add1stOrderTime( LocalMassMatrix, LocalStiffMatrix, &
              LocalForce, dt, N, 1, WPerm(NodeIndexes), Solver )
        END IF
        
        !------------------------------------------------------------------------------
        !      Update global matrix and rhs vector from local matrix & vector
        !------------------------------------------------------------------------------
        CALL UpdateGlobalEquations( StiffMatrix, LocalStiffMatrix, &
            ForceVector, LocalForce, n, 1, WPerm(NodeIndexes) )
        
        !------------------------------------------------------------------------------
        
      END DO !Assembly
      
      !------------------------------------------------------------------------------
      !    FinishAssembly must be called after all other assembly steps, but before
      !    Dirichlet boundary settings. Actually no need to call it except for
      !    transient simulations.
      !------------------------------------------------------------------------------
      CALL FinishAssembly( Solver,ForceVector )
      
      FirstTime = .FALSE.
      
      !------------------------------------------------------------------------------
      !  Built-in Dirichlet boundary conditions
      !------------------------------------------------------------------------------
      
      DO bc=1,Model % NumberOfBCs
        
        k = LumpedBoundaries(bc)
        IF(k == 0) CYCLE
        
        psi = FluidicForces(bc) / FluidicAreas(bc)
        A2 = FluidicAreas(bc)
        Q2 = FluidicFluxes(bc)
        
        Wbcnode = (Q2/A2) + 2 * SQRT(2.0/Density) * ( SQRT(psi + BetaCoeff * SQRT(Aref)) )
        
        IF ( StiffMatrix % FORMAT == MATRIX_SBAND ) THEN
          CALL SBand_SetDirichlet( StiffMatrix, ForceVector, k, Wbcnode )
        ELSE IF ( StiffMatrix % FORMAT == MATRIX_CRS .AND. StiffMatrix % Symmetric ) THEN
          CALL CRS_SetSymmDirichlet( StiffMatrix, ForceVector, k, Wbcnode )
        ELSE
          ForceVector(k) = Wbcnode
          CALL ZeroRow( StiffMatrix, k )
          CALL SetMatrixElement( StiffMatrix, k, k, 1.0d0 )
        END IF
        
      END DO
      
      !------------------------------------------------------------------------------
      !    Solve the system and we are done.
      !------------------------------------------------------------------------------
      Wprev = Wnodal
      Lprev = Lnodal
      Aprev = Anodal
      
      CALL SolveSystem( StiffMatrix, ParMatrix, ForceVector, &
          Wnodal, Norm, 1, Solver )
      
      !------------------------------------------------------------------------------
      !  Update the dofs dependent on Wnodal
      !------------------------------------------------------------------------------
      Anodal = ( (Density/BetaCoeff)**2 ) * (Wnodal - W20)**4 / (4**5)
      Unodal = ( Wnodal + W20 ) / 2.0
      Qnodal = Anodal * Unodal
      Cnodal = SQRT( BetaCoeff/(2.0*Density) ) * Anodal**0.25
      Lnodal = Unodal + Cnodal
      Pnodal = BetaCoeff * ( SQRT(Anodal) - SQRT(Aref) )
      
      !------------------------------------------------------------------------------
      !  Check the nonlinear convergence
      !------------------------------------------------------------------------------
      WdiffSum = SUM( ABS(Lnodal - Lprev) )
      WnodalSum = SUM( ABS(Lnodal) )
      WprevnodalSum = SUM( ABS(Lprev) )
      
      IF( WnodalSum + WprevnodalSum > TINY(WnodalSum)) THEN
        Werror = 2 * WdiffSum / (WnodalSum + WprevnodalSum)
      ELSE
        Werror = 0.0d0
      END IF
      
      WRITE(Message,'(A,E11.4)') 'Relative Change:',Werror
      CALL Info('OutletCompute',Message)
      
      IF ( Werror < NonlinearTol ) EXIT
      
    END DO ! Nonlinear iterations
    !------------------------------------------------------------------------------
    
  END DO ! Fractional-Step
  !------------------------------------------------------------------------------
  
  !------------------------------------------------------------------------------
  !  Save the new pressure and change in radius to be used in the outlet
  !  boundary (2D or 3D) of the Navier-Stokes equation.
  !------------------------------------------------------------------------------
  
  DO bc=1,Model % NumberOfBCs
    
    IF(LumpedBoundaries(bc) /= 0) THEN
      joinnode = LumpedBoundaries(bc)

      pres2dout = Pnodal(joinnode)
      AreaIn1d = Anodal(joinnode)
      RadiusIn1d = SQRT( AreaIn1d / PI )         
      dR = RadiusIn1d - ArteryRadius
      
      IF(Connections == 1) THEN
        Name = 'res: pout'
      ELSE
        WRITE(Name,'(A,I0)') 'res: pout',bc
      END IF
      CALL ListAddConstReal( Model % Simulation, Name, pres2dout )
      
      IF(Connections == 1) THEN
        Name = 'res: dRout'
      ELSE 
        WRITE(Name,'(A,I0)') 'res: dRout',bc
      END IF
      CALL ListAddConstReal( Model % Simulation, Name, dR )                 

    END IF
  END DO
  
CONTAINS
  
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix( StiffMatrix,MassMatrix,&
      Force, Lelem, Element, n, Nodes )
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: StiffMatrix(:,:), MassMatrix(:,:), &
        Force(:), Lelem(:)
    INTEGER :: n
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3)
    REAL(KIND=dp) :: SqrtElementMetric,U,V,W,S,Lambda
    LOGICAL :: Stat
    INTEGER :: i,p,q,t,DIM
    
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
    
!------------------------------------------------------------------------------
    DIM = CoordinateSystemDimension()
    Force = 0.0d0
    StiffMatrix = 0.0d0
    MassMatrix = 0.0d0

!------------------------------------------------------------------------------
!      Numerical integration
!------------------------------------------------------------------------------
    IntegStuff = GaussPoints( Element, n+1 )
    
    DO t=1, IntegStuff % n
      
      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)
      S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element,Nodes,U,V,W, SqrtElementMetric, &
          Basis, dBasisdx, ddBasisddx, .FALSE. )
      
      S = S * SqrtElementMetric
      Lambda = SUM( Lelem(1:n) * Basis(1:n) )

!------------------------------------------------------------------------------
!        The Poisson equation
!------------------------------------------------------------------------------
      DO p = 1,N
        Force(p) = Force(p) + L * Basis(p) * S
      END DO
      
      DO p = 1,N
        DO q = 1,N
          MassMatrix(p,q) = MassMatrix(p,q) + Basis(p) * Basis(q)*s
        END DO
      END DO
      
      DO p = 1,N
        DO q = 1,N
          DO i=1,DIM
            StiffMatrix(p,q) = StiffMatrix(p,q) + Lambda * Basis(p) &
                * dBasisdx(q,i) * S
          END DO
        END DO
      END DO
      
!------------------------------------------------------------------------------
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix


!------------------------------------------------------------------------------
!> Computes the center point of the higher dimensional boundary ends. 
!------------------------------------------------------------------------------
  SUBROUTINE SurfaceCenterPoints( SolidEndBoundaries, SolidConnections, SolidEndAreas)
!------------------------------------------------------------------------------

    IMPLICIT NONE
    
    INTEGER :: SolidEndBoundaries(:)
    REAL(KIND=dp) :: SolidEndAreas(:)
    LOGICAL :: Stat, LocalAllocationDone = .FALSE. 
    INTEGER :: i,j,k,n,pn,t,dim, SolidConnections
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Variable_t), POINTER :: Var
    TYPE(Nodes_t) :: ElementNodes, ParentNodes
    TYPE(Element_t), POINTER   :: Element, Parent
    CHARACTER(LEN=MAX_NAME_LEN) :: VariableName, XName, YName, ZName
    
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:), ddBasisddx(:,:,:), &
        ParentBasis(:), ParentdBasisdx(:,:), x(:), y(:), z(:)
    REAL(KIND=dp) :: u, v, w, s, detJ, xpos(3), SumOfWeights, XPosWeighted(3)
    REAL(KIND=dp) :: Normal(3)
    
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    
    SAVE LocalAllocationDone, ElementNodes, Basis, dBasisdx, x, y, z
    
    CALL Info( 'SurfaceCenterPoints', '-------------------------------------',Level=4 )
    CALL Info( 'SurfaceCenterPoints', 'Computing Centerpoints:  ', Level=4 )
    CALL Info( 'SurfaceCenterPoints', '-------------------------------------',Level=4 )
    
    
    IF(.NOT. LocalAllocationDone) THEN
      n = Mesh3D % MaxElementNodes
      
      ALLOCATE( ElementNodes % x(n), &
          ElementNodes % y(n),  &
          ElementNodes % z(n), &
          Basis(n), dBasisdx(n,3), &
          x(n), y(n), z(n) ) 
      LocalAllocationDone = .TRUE.         
    END IF
   
    
    DIM = CoordinateSystemDimension()
    
    DO k=1, Model % NumberOfBCs
      IF( SolidEndBoundaries(k) == 0 ) CYCLE

      SumOfWeights = 0.0_dp
      XposWeighted = 0.0_dp
      
      DO t = Mesh3D % NumberOfBulkElements + 1, &
          Mesh3D % NumberOfBulkElements + &
          Mesh3D % NumberOfBoundaryElements
        
        Element => Mesh3D % Elements(t)
        IF ( Element % TYPE % ElementCode == 101 ) CYCLE
        IF ( Model % BCs(k) % Tag /= Element % BoundaryInfo % Constraint ) CYCLE

    
        Model % CurrentElement => Mesh3D % Elements(t)
        n = Element % TYPE % NumberOfNodes
        NodeIndexes => Element % NodeIndexes
        
        ElementNodes % x(1:n) = Mesh3D % Nodes % x(NodeIndexes)
        ElementNodes % y(1:n) = Mesh3D % Nodes % y(NodeIndexes)
        ElementNodes % z(1:n) = Mesh3D % Nodes % z(NodeIndexes)
        
        IntegStuff = GaussPoints( Element )
        
!------------------------------------------------------------------------------
        DO l=1,IntegStuff % n
!------------------------------------------------------------------------------

          u = IntegStuff % u(l)
          v = IntegStuff % v(l)
          w = IntegStuff % w(l)
          s = IntegStuff % s(l) 

!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
          stat = ElementInfo( Element, ElementNodes, u, v, w, &
              detJ, Basis, dBasisdx)

!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
          s = s * detJ
    
          IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
            s = s * 2 * PI * SUM( ElementNodes % x(1:n) * Basis(1:n) )
          END IF
          
!------------------------------------------------------------------------------

          xpos(1) = SUM( ElementNodes % x(1:n) * Basis(1:n) )
          xpos(2) = SUM( ElementNodes % y(1:n) * Basis(1:n) )
          xpos(3) = SUM( ElementNodes % z(1:n) * Basis(1:n) )

          XPosWeighted = XPosWeighted + s * xpos 
          SumOfWeights = SumOfWeights + s

          SolidEndAreas(k) = SolidEndAreas(k) + s

!------------------------------------------------------------------------------
        END DO ! Integ
        
      END DO ! Boundary elements

!------------------------------------------------------------------------------

! Writing the center points of each solid end boundary to Model % Simulation

      XposWeighted = XposWeighted / SumOfWeights
      
      WRITE(Message,'(A,I0)') 'Center of Boundary ',k 
      CALL Info('SurfaceCenterPoints',Message)

      DO i=1,3
        IF( i == 1 ) THEN
          XName = 'res: xcenterpoint'
        ELSE IF( i == 2 ) THEN
          XName = 'res: ycenterpoint'
        ELSE
          XName = 'res: zcenterpoint'
        END IF
        IF(SolidConnections > 1) THEN
          WRITE(XName,'(A,I0)') TRIM(XName),k
        END IF

        CALL ListAddConstReal( Model % Simulation, XName, XPosWeighted(i) )

        WRITE(Message,'(A,I0,A,E11.4)') 'Coordinate ',i,':',XPosWeighted(i)  
        CALL Info('SurfaceCenterPoints',Message)
      END DO
!------------------------------------------------------------------------------
      
    END DO ! BCs

  END SUBROUTINE SurfaceCenterPoints
  

!---------------------------------------------------------------------------------------------
!> Computes the fluidic forces, fluxes and areas of the higher dimensional boundary.
!---------------------------------------------------------------------------------------------
  SUBROUTINE LumpedFluidicForce( LumpedBoundaries, FluidicForces, FluidicAreas, FluidicFluxes )
!---------------------------------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER :: LumpedBoundaries(:)
    REAL(KIND=dp) :: FluidicForces(:), FluidicAreas(:), FluidicFluxes(:)
    
    REAL(KIND=dp), ALLOCATABLE :: Pressure(:), Velocity(:,:), Viscosity(:)
    LOGICAL :: Stat, ViscousForce, Compressible, LocalAllocationDone = .FALSE. 
    INTEGER :: i,j,k,n,pn,t,dim
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Variable_t), POINTER :: Var
    TYPE(ValueList_t), POINTER :: Material
    TYPE(Nodes_t) :: ElementNodes, ParentNodes
    TYPE(Element_t), POINTER   :: Element, Parent
    CHARACTER(LEN=MAX_NAME_LEN) :: VariableName
    
    REAL(KIND=dp) :: Force(3), SumNormal(3), Area, LForce(3), Flux, NormalVelo
    
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:), &
        ParentBasis(:), ParentdBasisdx(:,:), x(:), y(:), z(:)
    REAL(KIND=dp) :: u, v, w, s, detJ, xpos
    REAL(KIND=dp) :: Grad(3,3), Stress(3,3), Normal(3), Div, Visc
    
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    
    SAVE LocalAllocationDone, ElementNodes, ParentNodes, Basis, dBasisdx, &
     ParentBasis, ParentdBasisdx, x, y, z, Pressure, Viscosity, Velocity
    
    CALL Info( 'LumpedForceCompute', '-------------------------------------',Level=4 )
    CALL Info( 'LumpedForceCompute', 'Computing Lumped Boundary:  ', Level=4 )
    CALL Info( 'LumpedForceCompute', '-------------------------------------',Level=4 )

    FluidicForces = 0.0_dp
    FluidicAreas = 0.0_dp
    FluidicFluxes = 0.0_dp

    CALL SetCurrentMesh( Model, Mesh3D )
    
    VariableName = GetString( Solver % Values, 'Velocity Field Name', stat )
    IF ( .NOT. stat )  THEN
      Var => VariableGet( Mesh3D % Variables, 'Flow Solution', .TRUE. )
    ELSE
      Var => VariableGet( Mesh3D % Variables, VariableName, .TRUE. )
    END IF
    
    IF(.NOT. LocalAllocationDone) THEN
      n = Mesh3D % MaxElementNodes
      
      ALLOCATE( ElementNodes % x(n), &
          ElementNodes % y(n),  &
          ElementNodes % z(n), &
          ParentNodes % x(n), &
          ParentNodes % y(n),  &
          ParentNodes % z(n), &
          Basis(n), dBasisdx(n,3), &
          ParentBasis(n), ParentdBasisdx(n,3), &
          x(n), y(n), z(n), &
          Pressure( n ), &
          Viscosity( n ), &
          Velocity( 3, n ) )
      LocalAllocationDone = .TRUE.         
    END IF
    
    
    DIM = CoordinateSystemDimension()
    
    ViscousForce = ListGetLogical( Solver % Values, 'Consider Viscous Force',stat )
    IF ( .NOT. stat )  ViscousForce = .FALSE.    
    Compressible = ListGetLogical( Solver % Values, 'Consider Compressibility',stat )
        
    DO k=1, Model % NumberOfBCs
      IF( LumpedBoundaries(k) == 0 ) CYCLE
      
      Force = 0.0d0
      SumNormal = 0.0d0
      Area = 0.0d0
      Flux = 0.0d0
      
      DO t = Mesh3D % NumberOfBulkElements + 1, &
          Mesh3D % NumberOfBulkElements + &
          Mesh3D % NumberOfBoundaryElements
        
        Element => Mesh3D % Elements(t)
        IF ( Element % TYPE % ElementCode == 101 ) CYCLE
        IF ( Model % BCs(k) % Tag /= Element % BoundaryInfo % Constraint ) CYCLE
        
        Model % CurrentElement => Mesh3D % Elements(t)
        n = Element % TYPE % NumberOfNodes
        NodeIndexes => Element % NodeIndexes
        
        ElementNodes % x(1:n) = Mesh3D % Nodes % x(NodeIndexes)
        ElementNodes % y(1:n) = Mesh3D % Nodes % y(NodeIndexes)
        ElementNodes % z(1:n) = Mesh3D % Nodes % z(NodeIndexes)
        
        Parent => Element % BoundaryInfo % Left
        
        stat = ASSOCIATED( Parent )
        IF ( stat ) stat = stat .AND. ALL(Var % Perm(Parent % NodeIndexes(1:n)) > 0)
        
        IF ( .NOT. stat ) THEN
          Parent => ELement % BoundaryInfo % Right
          
          stat = ASSOCIATED( Parent )
          IF ( stat ) stat = ALL(Var % Perm(Parent % NodeIndexes(1:n)) > 0)
          
          IF ( .NOT. stat )  THEN
            CALL Warn( 'LumpedForceCompute', &
                'No flow solution available for specified boundary' )
            CYCLE
          END IF
        END IF
        
        pn = Parent % TYPE % NumberOfNodes           
        ParentNodes % x(1:pn) = Mesh3D % Nodes % x(Parent % NodeIndexes)
        ParentNodes % y(1:pn) = Mesh3D % Nodes % y(Parent % NodeIndexes)
        ParentNodes % z(1:pn) = Mesh3D % Nodes % z(Parent % NodeIndexes)
        
        j = ListGetInteger( Model % Bodies(Parent % BodyId) % Values, 'Material', &
            minv=1, maxv=Model % NumberOFMaterials )
        Material => Model % Materials(j) % Values
        
        Viscosity(1:pn) = ListGetReal( Material, 'Viscosity', pn, Parent % NodeIndexes )
        
        Velocity = 0.0d0
        DO i=1,pn
          DO j=1,DIM
            Velocity(j,i) = &
                Var % Values(Var % DOFs * (Var % Perm(Parent % NodeIndexes(i))-1)+j)
          END DO
        END DO
        Pressure(1:pn) = Var % Values(Var % DOFs * Var % Perm(Parent % NodeIndexes))
        
        IntegStuff = GaussPoints( Element )
        
!------------------------------------------------------------------------------
        DO l=1,IntegStuff % n
!------------------------------------------------------------------------------

          u = IntegStuff % u(l)
          v = IntegStuff % v(l)
          w = IntegStuff % w(l)
          s = IntegStuff % s(l) 

!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------

          stat = ElementInfo( Element, ElementNodes, u, v, w, detJ, Basis, dBasisdx)

!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
          s = s * detJ
          
          IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
            xpos = SUM( ElementNodes % x(1:n) * Basis(1:n) )
            s = s * 2.0 * PI * xpos
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
      
    
!------------------------------------------------------------------------------

          Stress = 0.0d0
          Div = 0.0d0
          Visc = 0.0d0
          
          NormalVelo = 0.0d0
          DO j = 1,DIM
            NormalVelo = NormalVelo + Normal(j) * SUM( ParentBasis(1:pn) * Velocity(j,1:pn) )
          END DO
          
          IF ( ViscousForce ) THEN           
            Grad = MATMUL( Velocity(:,1:pn),ParentdBasisdx )
            Visc = SUM( Viscosity(1:pn) * ParentBasis(1:pn) )
            
            IF ( Compressible ) THEN             
              IF ( CurrentCoordinateSystem() == Cartesian ) THEN
                DO i = 1, DIM
                  Div = Div + Grad(i,i)
                END DO
              ELSE
                Div = SUM( Velocity(1,1:pn) * ParentdBasisdx(1:pn,1) ) + &
                    SUM( Velocity(1,1:pn) * ParentBasis(1:pn) ) / xpos + &
                    SUM( Velocity(2,1:pn) * ParentdBasisdx(1:pn,2) )
              END IF
            END IF
            
            Stress = Visc * ( Grad + TRANSPOSE(Grad) )           
          END IF
          
          DO i=1,DIM
            Stress(i,i) = Stress(i,i) - SUM( Pressure(1:pn) * ParentBasis ) &
                -(2.0d0/3.0d0) * Visc * Div
          END DO
          
          LForce = -MATMUL( Stress, Normal )
          
          Force  = Force  + s * LForce
          Flux = Flux + s * NormalVelo
          SumNormal = SumNormal + s * Normal
          
          Area = Area + s
 
!------------------------------------------------------------------------------
        END DO

      END DO
      
      FluidicAreas(k) = SQRT(SUM(SumNormal * SumNormal))
      FluidicForces(k) = SUM(Force*SumNormal)/FluidicAreas(k)       
      FluidicFluxes(k) = Flux
      
    END DO

    CALL SetCurrentMesh( Model, Mesh1D )
    
  END SUBROUTINE LumpedFluidicForce
  

!------------------------------------------------------------------------------
END SUBROUTINE OutletCompute
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Compute the initial guess for the characteristics variable
!> \ingroup UDF
!------------------------------------------------------------------------------
FUNCTION OutletInit( Model,n,t ) RESULT( Winit )
  
  USE Types
  USE Lists
  
  IMPLICIT NONE
  
  TYPE(Model_t) :: Model
  INTEGER :: n,mat_id
  REAL( kind=dp ) :: t, Winit  
  TYPE(ValueList_t), POINTER :: Material
  REAL( kind=dp ) :: Aref, BetaCoeff, W20, rho, Ew, r, hw, nu
  
  mat_id = ListGetInteger( Model % Bodies(  &
      Model % CurrentElement % BodyId ) % Values, 'Material' )
  Material => Model % Materials(mat_id) % Values
  
  rho = ListGetConstReal( Material, 'Density')
  Ew = ListGetConstReal( Material, 'Artery Wall Youngs Modulus')
  R = ListGetConstReal( Material, 'Artery Radius')
  hw = ListGetConstReal( Material, 'Artery Wall Thickness')
  nu = ListGetConstReal( Material, 'Artery Poisson Ratio')
  
  Aref = pi * R**2
  BetaCoeff = ( 1 / (1 - nu**2) ) * ( SQRT(pi) * hw  * Ew ) / ( Aref)

  Winit = 4 * ( SQRT( (BetaCoeff/(2.0*rho)) * SQRT(Aref)) )
  
END FUNCTION OutletInit


!------------------------------------------------------------------------------
!> Return the change in the radius x-component computed by the 1D model
!> \ingroup UDF
!------------------------------------------------------------------------------
FUNCTION OutletdX( Model,n,t ) RESULT(dx)
  
  USE Types
  USE Lists
  USE MeshUtils
  
  IMPLICIT NONE
  TYPE(Model_t) :: Model  
  INTEGER :: n, bc, bc2
  REAL( kind=dp ) :: t, dRout, xcenterpoint, ycenterpoint, zcenterpoint, &
           x, y, z, x0, y0, dx, dy
  LOGICAL :: GotIt,GotIt1,GotIt2,GotIt3
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName, VarName1, VarName2, VarName3
  
  !------------------------------------------------------------------------------
  ! Get the requested change in radius of the higher dimensional boundary
  !------------------------------------------------------------------------------
  dRout = ListGetConstReal( Model % Simulation, 'res: dRout', GotIt)

  IF(.NOT. GotIt) THEN    
    bc = Model % CurrentElement % BoundaryInfo % Constraint
    bc2 = ListGetInteger( Model % BCs(bc) % Values,'Structure Coupling With Boundary', GotIt)
    IF(GotIt) bc = bc2
    WRITE(VarName,'(A,I0)') 'res: dRout ',bc
    dRout = ListGetConstReal( Model % Simulation, VarName, GotIt)  
  END IF

  ! If this is axisymmetric coordinate system we're done since for that the system
  ! always is centered around origin.
  !-------------------------------------------------------------------------------
  IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
    dx = dRout
    RETURN
  END IF


  !------------------------------------------------------------------------------
  ! Getting the center point Model % Simulation 
  !------------------------------------------------------------------------------
  xcenterpoint = ListGetConstReal( Model % Simulation, 'res: xcenterpoint', GotIt1)
  ycenterpoint = ListGetConstReal( Model % Simulation, 'res: ycenterpoint', GotIt2)
  zcenterpoint = ListGetConstReal( Model % Simulation, 'res: zcenterpoint', GotIt3)

  IF(.NOT. (GotIt1 .AND. GotIt2 .AND. GotIt3) ) THEN    
    bc = Model % CurrentElement % BoundaryInfo % Constraint
    bc2 = ListGetInteger( Model % BCs(bc) % Values,'Structure Coupling With Boundary', GotIt)
    IF(GotIt) bc = bc2
    WRITE(VarName1,'(A,I0)') 'res: xcenterpoint ',bc
    WRITE(VarName2,'(A,I0)') 'res: ycenterpoint ',bc
    WRITE(VarName3,'(A,I0)') 'res: zcenterpoint ',bc
    xcenterpoint= ListGetConstReal( Model % Simulation, VarName1, GotIt)  
    ycenterpoint= ListGetConstReal( Model % Simulation, VarName2, GotIt)  
    zcenterpoint= ListGetConstReal( Model % Simulation, VarName3, GotIt)  
  END IF

  !------------------------------------------------------------------------------
  ! Point under inspection.
  !------------------------------------------------------------------------------
  x = Model % Nodes % x(n) 
  y = Model % Nodes % y(n) 
  z = Model % Nodes % z(n) 
  
  x0 = x - xcenterpoint
  y0 = y - ycenterpoint
  
  dx = SQRT( 1./(1 + (y0/x0)**2) ) * dRout
  dy = SQRT( 1./(1 + (x0/y0)**2) ) * dRout
  
  IF( x0 <= 0 ) THEN
    dx = -dx
  END IF
  
  IF( y0 <= 0 ) THEN
    dy = -dy
  END IF


END FUNCTION OutletdX


!------------------------------------------------------------------------------
!> Return the change in the radius y-component computed by the 1D model
!> \ingroup UDF
!------------------------------------------------------------------------------
FUNCTION OutletdY( Model,n,t ) RESULT(dy)
  
  USE Types
  USE Lists
  USE MeshUtils
  
  IMPLICIT NONE
  TYPE(Model_t) :: Model  
  INTEGER :: n, bc, bc2
  REAL( kind=dp ) :: t, dRout, xcenterpoint, ycenterpoint, zcenterpoint, &
           x, y, z, x0, y0, dx, dy
  LOGICAL :: GotIt,GotIt1,GotIt2,GotIt3
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName, VarName1, VarName2, VarName3

  IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
    CALL Fatal('OutletDy','This BC does not make sense in axisymmetric geometry')
  END IF

  !------------------------------------------------------------------------------
  ! Getting the center points and displacement dR from Model % Simulation 

  xcenterpoint = ListGetConstReal( Model % Simulation, 'res: xcenterpoint', GotIt1)
  ycenterpoint = ListGetConstReal( Model % Simulation, 'res: ycenterpoint', GotIt2)
  zcenterpoint = ListGetConstReal( Model % Simulation, 'res: zcenterpoint', GotIt3)

  IF(.NOT. (GotIt1 .AND. GotIt2 .AND. GotIt3) ) THEN    
    bc = Model % CurrentElement % BoundaryInfo % Constraint
    bc2 = ListGetInteger( Model % BCs(bc) % Values,'Structure Coupling With Boundary', GotIt)
    IF(GotIt) bc = bc2
    WRITE(VarName1,'(A,I0)') 'res: xcenterpoint ',bc
    WRITE(VarName2,'(A,I0)') 'res: ycenterpoint ',bc
    WRITE(VarName3,'(A,I0)') 'res: zcenterpoint ',bc
    xcenterpoint= ListGetConstReal( Model % Simulation, VarName1, GotIt)  
    ycenterpoint= ListGetConstReal( Model % Simulation, VarName2, GotIt)  
    zcenterpoint= ListGetConstReal( Model % Simulation, VarName3, GotIt)  
  END IF

  dRout = ListGetConstReal( Model % Simulation, 'res: dRout', GotIt)  

  IF(.NOT. GotIt) THEN    
    bc = Model % CurrentElement % BoundaryInfo % Constraint
    bc2 = ListGetInteger( Model % BCs(bc) % Values,'Structure Coupling With Boundary', GotIt)
    IF(GotIt) bc = bc2
    WRITE(VarName,'(A,I0)') 'res: dRout ',bc
    dRout = ListGetConstReal( Model % Simulation, VarName, GotIt)  
  END IF

!------------------------------------------------------------------------------

  x = Model % Nodes % x(n) 
  y = Model % Nodes % y(n) 
  z = Model % Nodes % z(n) 
  
  x0 = x - xcenterpoint
  y0 = y - ycenterpoint
  
  dx = SQRT( 1./(1 + (y0/x0)**2) ) * dRout
  dy = SQRT( 1./(1 + (x0/y0)**2) ) * dRout
  
  IF( x0 <= 0 ) THEN
    dx = -dx
  END IF
  
  IF( y0 <= 0 ) THEN
    dy = -dy
  END IF
  

END FUNCTION OutletdY



!------------------------------------------------------------------------------
!> Return the pressure computed by the characteristics model.
!> \ingroup UDF
!------------------------------------------------------------------------------
FUNCTION OutletPres( Model,n,t ) RESULT(pout)
  
  USE Types
  USE Lists
  
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER :: n, bc
  REAL( kind=dp ) :: t, pout
  LOGICAL :: GotIt
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName

  pout = ListGetConstReal( Model % Simulation, 'res: pout', GotIt)  
  IF(.NOT. GotIt) THEN
    bc = Model % CurrentElement % BoundaryInfo % Constraint
    WRITE(VarName,'(A,I0)') 'res: pout ',bc
    pout = ListGetConstReal( Model % Simulation, VarName, GotIt)  
  END IF

  pout = -pout

END FUNCTION OutletPres

!------------------------------------------------------------------------------
