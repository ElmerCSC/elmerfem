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
! *  Module containing a solver that normalized the artificial compressibility for
! *  optimal fluid-structure coupling and for for computing the artificial 
! *  compressibility elementwise.
! *
! ******************************************************************************/
! *
! *  Authors: Peter Raback
! *  Email:   Peter.Rabackcsc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 12.02.2002
! *  Modified: 26.4.2022 (for Couldron simulations)
! * 
! ******************************************************************************

!------------------------------------------------------------------------------
!> Subroutine that may be used to normalize the amplitude of artificial compressibility for
!> optimal fluid-structure coupling.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE CompressibilityScale( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver   !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model     !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt        !< Timestep size for time dependent simulations
  LOGICAL :: Transient !< Steady state or transient simulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  REAL(KIND=dp), ALLOCATABLE :: Pressure(:), &
      Displacement(:,:), Compressibility(:)
  TYPE(Solver_t), POINTER :: FlowSolver
  REAL(KIND=dp) :: SideVolume, SidePressure, SideArea, TotalVolume, &
      InitVolume, TotalVolumeCompress, CompressSuggest, &
      CompressScale, CompressScaleOld, Relax, Norm, &
      Timestep, PrevTimestep=-1, MinimumSideVolume=1.0d10, &
      TransitionVolume, dVolume
  LOGICAL :: Stat, GotIt, SubroutineVisited=.False., &
      ScaleCompressibility, WeightByDisplacement, PressureMode
  INTEGER :: i,j,k,n,pn,t,dim,mat_id,TimeStepVisited,istat
  INTEGER, POINTER :: NodeIndexes(:)
  TYPE(Variable_t), POINTER :: FVar,Dvar
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: BC, Material, Params
  TYPE(Nodes_t) :: ElementNodes
  TYPE(Element_t), POINTER   :: Element
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName
!------------------------------------------------------------------------------

  SAVE SubroutineVisited, TimeStepVisited, PrevTimestep, InitVolume, &
      MinimumSideVolume

  CALL Info('CompressibilityScale','Scaling compressibility for optimal FSI convergence',Level=5)
  
  timestep = GetTimestep()
  IF( timestep /= prevTimestep ) THEN
    TimeStepVisited = 0
    PrevTimestep = Timestep
  END IF

  Params => GetSolverParams()
  
  Mesh => Solver % Mesh

  FVar => NULL()
  DO i=1,Model % NumberOfSolvers
    FlowSolver => Model % Solvers(i)
    IF ( ListGetString( FlowSolver % Values, 'Equation' ) == 'navier-stokes' ) THEN
      FVar => FlowSolver % Variable
      EXIT
    END IF
  END DO
  
  IF(.NOT. ASSOCIATED( FVar ) ) THEN
    VarName = ListGetString( Params,'Flow Variable Name',GotIt)
    IF(GotIt) FVar => VariableGet( Mesh % Variables, VarName )
  END IF

  IF(.NOT. ASSOCIATED( FVar ) ) THEN
    CALL Fatal('CompressibilityScale','Please give valid "Flow Variable Name"!')
  END IF

  PressureMode = ( FVar % Dofs == 1 )
  IF( PressureMode ) THEN
    CALL Info('CompressibilityScale','Using fluid solver for just pressure',Level=10)
  END IF
  
  VarName = ListGetString( Params,'Displacement Variable Name',GotIt)
  IF(.NOT. GotIt) VarName = 'Displacement'
  Dvar => VariableGet( Mesh % Variables, VarName, .TRUE.)
  IF(.NOT. ASSOCIATED(DVar)) THEN
    CALL Fatal('CompressibilityScale','Please give valid "Displacement Variable Name"!')
  END IF
   
  n = Mesh % MaxElementNodes
  ALLOCATE( ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n),&
      Compressibility(n), Pressure(n), Displacement( 3,n),STAT=istat ) 
  IF ( istat /= 0 ) CALL Fatal( 'CompressibilityScale', 'Memory allocation error.' )
  
  TotalVolume = 0.0_dp
  TotalVolumeCompress = 0.0_dp
  SidePressure = 0.0_dp
  SideVolume = 0.0_dp
  SideArea = 0.0_dp
  
  dim = CoordinateSystemDimension()

  WeightByDisplacement = ListGetLogical( Params,'Weight By Displacement',GotIt)
  

  IF(.NOT. PressureMode ) THEN
    DO t=1,Mesh % NumberOfBulkElements
      
      Element => Mesh % Elements(t)
      Model % CurrentElement => Element
      
      n = Element % TYPE % NumberOfNodes
      NodeIndexes => Element % NodeIndexes

      IF( ANY ( FVar % Perm( NodeIndexes ) == 0 ) ) CYCLE

      ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
      ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
      ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)

      mat_id = ListGetInteger( Model % Bodies(Element % BodyId) % &
          Values, 'Material',minv=1,maxv=Model % NumberOfMaterials )

      Material => Model % Materials(mat_id) % Values

      Compressibility(1:n) = &
          ListGetReal(Material,'Artificial Compressibility',n,NodeIndexes,gotIt)
      IF(.NOT. GotIt) CYCLE
      
      CALL CompressibilityIntegrate(Element, n, ElementNodes, &
          Compressibility, TotalVolume, TotalVolumeCompress)
    END DO
  END IF

  
  ! Compute the force acting on the boundary
  DO t = Mesh % NumberOfBulkElements + 1, &
      Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
    
!------------------------------------------------------------------------------
    Element => Mesh % Elements(t)
    IF ( Element % TYPE % ElementCode == 101 ) CYCLE
          
!------------------------------------------------------------------------------
!    Set the current element pointer in the model structure to 
!    reflect the element being processed
!------------------------------------------------------------------------------
    Model % CurrentElement => Element
!------------------------------------------------------------------------------
    n = Element % TYPE % NumberOfNodes
    NodeIndexes => Element % NodeIndexes
     
    ! Only integrate over BCs shared with both solvers
    IF( ANY ( DVar % Perm( NodeIndexes ) == 0 ) ) CYCLE
    IF( ANY ( FVar % Perm( NodeIndexes ) == 0 ) ) CYCLE
     
    BC => GetBC(Element)
    IF(.NOT. ASSOCIATED(BC)) CYCLE
    
    ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
    ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
    ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)

    Pressure(1:n) = FVar % Values(FVar % DOFs * FVar % Perm(NodeIndexes))
      
    Displacement = 0.0_dp
    DO j=1,dim
      Displacement(j,1:n) = &
          DVar % Values(DVar % dofs * (DVar % Perm(NodeIndexes(1:n))-1)+j)
    END DO
    
    IF( PressureMode ) THEN
      Material => GetMaterial(Element)
    END IF

    CALL PressureIntegrate(Element, n, ElementNodes)
  END DO

  ! Compute the initial volume from the 1st volume
  IF(ABS(SideVolume) < MinimumSideVolume) THEN
    MinimumSideVolume = ABS(SideVolume)
    InitVolume = TotalVolume-SideVolume 
  END IF

  CompressScaleOld = ListGetConstReal( Model % Simulation, &
      'Artificial Compressibility Scaling',GotIt)
  IF(.NOT. GotIt) CompressScaleOld = 1.0

  Relax = GetCReal( &
      Solver % Values, 'Nonlinear System Relaxation Factor',gotIt )
  IF(.NOT. gotIt) Relax = 1.0;

  TransitionVolume = ListGetConstReal( &
      Solver % Values, 'Artificial Compressibility Critical Volume',gotIt )
  IF(.NOT. gotIt) TransitionVolume = 0.01;
 
  ScaleCompressibility = ListGetLogical( &
      Solver % Values, 'Artificial Compressibility Scale',gotIt )
  IF(.NOT. gotIt) ScaleCompressibility = .TRUE.

  IF(SideVolume/TotalVolume > TransitionVolume) THEN
    dVolume = TotalVolume - InitVolume
  ELSE
    dVolume = SideVolume
  END IF
  CompressSuggest = (dVolume/TotalVolume)/(SidePressure * SideArea) 
  CompressScale = CompressSuggest*TotalVolume/ (TotalVolumeCompress*Relax)


  Norm = CompressScale
  IF(TimeStepVisited == 0) Norm = Norm * 2.0
  Solver % Variable % Norm = Norm

  TimeStepVisited = TimeStepVisited + 1

  IF(ScaleCompressibility) THEN
    CALL ListAddConstReal( Model % Simulation, &
        'Artificial Compressibility Scaling',CompressScale)
  END IF

  CALL ListAddConstReal( Model % Simulation, &
      'res: Relative Volume Change',dVolume/InitVolume)
  CALL ListAddConstReal( Model % Simulation, &
      'res: Mean Pressure on Surface',SidePressure/SideArea)
  CALL ListAddConstReal( Model % Simulation, &
      'res: Suggested Compressibility',CompressSuggest)
  !CALL ListAddConstReal( Model % Simulation, &
  !    'res: Iterations',1.0_dp*TimeStepVisited)

  WRITE(Message,'(A,T25,E15.4)') 'Relative Volume Change',dVolume/InitVolume
  CALL Info('CompressibilityScale',Message,Level=5)
  WRITE(Message,'(A,T25,E15.4)') 'Mean Pressure on Surface',SidePressure/SideArea
  CALL Info('CompressibilityScale',Message,Level=5)
  WRITE(Message,'(A,T25,E15.4)') 'Suggested Compressibility',CompressSuggest
  CALL Info('CompressibilityScale',Message,Level=5)
  WRITE(Message,'(A,T25,E15.4)') 'Compressibility Scaling Factor',&
      CompressScale/CompressScaleOld
  CALL Info('CompressibilityScale',Message,Level=5)

  DEALLOCATE( ElementNodes % x, ElementNodes % y, ElementNodes % z, &
      Pressure, Displacement, Compressibility)

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE PressureIntegrate(Element, n, Nodes)
!------------------------------------------------------------------------------
    INTEGER :: n
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element

!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3)
    REAL(KIND=dp) :: u,v,w,s,x,y,z, Pres
    REAL(KIND=dp) :: Grad(3,3), Stress(3,3), Normal(3), Ident(3,3)
    REAL(KIND=dp) :: NormalDisplacement,Symb(3,3,3),dSymb(3,3,3,3)
    REAL(KIND=dp) :: SqrtElementMetric, SqrtMetric, Metric(3,3)
    REAL(KIND=dp) :: ElemArtif(n), ElemGap(n), ElemPres0(n), AC, Gap, Pres0
    LOGICAL :: SurfAC
    
    INTEGER :: N_Integ, CoordSys
    
    LOGICAL :: stat
    INTEGER :: i,t
    TYPE(GaussIntegrationPoints_t), TARGET :: IP
    
    Ident = 0.0_dp
    DO i=1,3
      Ident(i,i) = 1.0_dp
    END DO
    
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
    IP = GaussPoints( Element )
    
    CoordSys = CurrentCoordinateSystem()

    SurfAC = .FALSE.
    IF( PressureMode ) THEN
      ElemArtif(1:n) = GetReal( Material,'Artificial Compressibility',GotIt)
      IF(GotIt) THEN
        ElemGap(1:n) = GetReal(Material,'Gap Height')
      ELSE
        ElemArtif(1:n) = GetReal( Material,'Surface Compressibility',SurfAC)
      END IF
    END IF
    
    ElemPres0(1:n) = GetReal( Material,'Equilibrium Pressure',GotIt) 

!------------------------------------------------------------------------------
    DO t=1,IP % n
!------------------------------------------------------------------------------
       u = IP % u(t)
       v = IP % v(t)
       w = IP % w(t)
!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, u, v, w, &
           SqrtElementMetric, Basis, dBasisdx )
       
       s = SqrtElementMetric * IP % s(t)

       IF ( CoordSys /= Cartesian ) THEN
         X = SUM( Nodes % X(1:n) * Basis(1:n) )
         Y = SUM( Nodes % Y(1:n) * Basis(1:n) )
         Z = SUM( Nodes % Z(1:n) * Basis(1:n) )
         CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,X,Y,Z )
         s = s * SqrtMetric
       END IF
       
       Normal = Normalvector( Element,Nodes, u,v, .TRUE. )
       
       NormalDisplacement = 0.0
       DO i=1,dim
         NormalDisplacement = NormalDisplacement + &
             SUM(Basis(1:n) * Displacement(i,1:n)) * Normal(i)
       END DO
       
       Pres0 = SUM(Basis(1:n) * ElemPres0(1:n) )
       Pres = SUM( Basis(1:n) * Pressure(1:n))

       IF( WeightByDisplacement ) THEN
         s = s * ABS( NormalDisplacement )
       END IF
              
       SideVolume = SideVolume + s * ABS(NormalDisplacement)
       SidePressure = SidePressure + s * ABS(Pres-Pres0)
       SideArea = SideArea + s

       IF( PressureMode ) THEN
         AC = SUM(Basis(1:n) * ElemArtif(1:n) )
         IF(SurfAC) THEN
           Gap = 1.0_dp
         ELSE
           Gap = SUM(Basis(1:n) * ElemGap(1:n) ) 
         END IF
         TotalVolume = TotalVolume + s*Gap
         TotalVolumeCompress = TotalVolumeCompress + s*AC*Gap
       END IF
       
!------------------------------------------------------------------------------
     END DO
!------------------------------------------------------------------------------
  END SUBROUTINE PressureIntegrate
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE CompressibilityIntegrate(Element, n, Nodes, &
      Compressibility, TotalVolume, TotalVolumeCompress)
!------------------------------------------------------------------------------
    INTEGER :: n
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: Compressibility(:), TotalVolume, TotalVolumeCompress
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3)
    REAL(KIND=dp) :: SqrtElementMetric, U, V, W, S, C
    LOGICAL :: Stat
    INTEGER :: i,p,q,t,dim, NBasis, CoordSys
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: X,Y,Z,Metric(3,3),SqrtMetric,Symb(3,3,3),dSymb(3,3,3,3)

!------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()
    CoordSys = CurrentCoordinateSystem()

    Metric = 0.0_dp
    Metric(1,1) = 1.0_dp
    Metric(2,2) = 1.0_dp
    Metric(3,3) = 1.0_dp

!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------

    NBasis = n
    IP = GaussPoints( Element )

!------------------------------------------------------------------------------
    DO t=1,IP % n
      U = IP % u(t)
      V = IP % v(t)
      W = IP % w(t)

!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, U, V, W, SqrtElementMetric, &
          Basis, dBasisdx )
      s = IP % s(t) * SqrtElementMetric

      IF ( CoordSys /= Cartesian ) THEN
        X = SUM( Nodes % X(1:n) * Basis(1:n) )
        Y = SUM( Nodes % Y(1:n) * Basis(1:n) )
        Z = SUM( Nodes % Z(1:n) * Basis(1:n) )
        CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,X,Y,Z )
        s = s * SqrtMetric
      END IF
 
      C = SUM(Basis(1:n) * Compressibility(1:n))

      TotalVolume = TotalVolume + s
      TotalVolumeCompress = TotalVolumeCompress + s*C
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE CompressibilityIntegrate
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE CompressibilityScale
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE CompressibilityScale_init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver   
  TYPE(Model_t) :: Model     
  REAL(KIND=dp) :: dt        
  LOGICAL :: Transient

  CALL ListAddConstReal( Model % Simulation, &
      'res: Relative Volume Change',0.0_dp) 
  CALL ListAddConstReal( Model % Simulation, &
      'res: Mean Pressure on Surface',0.0_dp)
  CALL ListAddConstReal( Model % Simulation, &
      'res: Suggested Compressibility',0.0_dp)
  !CALL ListAddConstReal( Model % Simulation, &
  !    'res: Iterations',0.0_dp)

END SUBROUTINE CompressibilityScale_Init



!------------------------------------------------------------------------------
!> Subroutine for computing the artificial compressibility from the volume change 
!> of elements. The volume change is obtained by extending the displacement field
!> of a test load to the fluid domain.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE CompressibilitySolver( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver   !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model     !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt        !< Timestep size for time dependent simulations
  LOGICAL :: Transient !< Steady state or transient simulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Matrix_t), POINTER  :: StiffMatrix
  TYPE(Variable_t), POINTER :: DisplacementSol, PressureSol
  TYPE(Element_t), POINTER :: Element
  TYPE(Mesh_t), POINTER :: Mesh
  
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName

  REAL(KIND=dp), POINTER :: CompressibilityFunction(:), DisplacementSolValues(:), &
      PressureSolValues(:)
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:), ElemDisp(:,:), ElemPres(:)
  REAL(KIND=dp) :: Norm, ReferencePressure

  INTEGER, POINTER :: NodeIndexes(:), Perm(:), DisplacementSolPerm(:), PressureSolPerm(:)
  INTEGER :: k, t, i, j, n, istat, NSDOFs, dim

  TYPE(ValueList_t), POINTER :: Params
  TYPE(Nodes_t) :: Nodes0, Nodes1
 
  LOGICAL :: AllocationsDone = .FALSE., Found, DisplacedShape, PressureExists

  SAVE STIFF, FORCE, Nodes0, Nodes1, ElemPres, ElemDisp, AllocationsDone

  CALL Info('CompressibilitySolver',' ')
  CALL Info('CompressibilitySolver','----------------------------')
  CALL Info('CompressibilitySolver','Compressibility Solver')
  CALL Info('CompressibilitySolver','----------------------------')
  CALL Info('CompressibilitySolver',' ')

!------------------------------------------------------------------------------
! Get variables needed for solution
!------------------------------------------------------------------------------
  IF ( .NOT.ASSOCIATED( Solver % Matrix ) ) RETURN

  Mesh => Solver % Mesh

  CompressibilityFunction => Solver % Variable % Values
  Perm => Solver % Variable % Perm
  
  IF(Transient) THEN
    CALL Warn('CompressibilitySolver','Implemented only for steady state')
    PRINT *,'AC interval',MINVAL(CompressibilityFunction),MAXVAL(CompressibilityFunction)
    RETURN
  END IF

  dim = CoordinateSystemDimension()

  StiffMatrix => Solver % Matrix
!------------------------------------------------------------------------------
! Get initial values ( Default is DisplacementSolvers 'Displacement' )
!------------------------------------------------------------------------------
  Params => GetSolverParams()
  VarName = GetString( Params,'Displacement Variable Name', Found )
  IF ( .NOT. Found ) VarName = 'Displacement'
  DisplacementSol => VariableGet( Mesh % Variables, VarName )
  IF(.NOT. ASSOCIATED( DisplacementSol ) ) THEN
    CALL Fatal( 'CompressibilitySolver', 'No variable for displacement associated!' )
  END IF

  DisplacementSolPerm => DisplacementSol % Perm
  DisplacementSolValues => DisplacementSol % Values
  NSDOFs = DisplacementSol % DOFs

  VarName = GetString( Params,'Pressure Variable Name', Found )
  IF ( .NOT. Found ) VarName = 'Pressure'
  PressureSol => VariableGet( Mesh % Variables, VarName )
  PressureExists = ASSOCIATED(PressureSol)
  IF(PressureExists) THEN
    PressureSolPerm => PressureSol % Perm
    PressureSolValues => PressureSol % Values
    PressureExists = .FALSE.
  ELSE
    ReferencePressure = ListGetConstReal( Solver % Values,'Reference Pressure',Found)
    IF(.NOT. Found) THEN
      CALL Fatal( 'CompressibilitySolver', 'No variable for "pressure" or "Reference Pressure" available!' )
    END IF
  END IF

!------------------------------------------------------------------------------
! Get keyword values
!------------------------------------------------------------------------------

  DisplacedShape = ListGetLogical(Solver % Values,'Displaced Shape', Found)
  
!------------------------------------------------------------------------------
! Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
  IF ( .NOT. AllocationsDone ) THEN
    n = Mesh % MaxElementNodes ! just big enough for elemental arrays
    ALLOCATE( FORCE( n ), STIFF(n,n), ElemDisp(3,n), &
        Nodes0 % X(n), Nodes0 % Y(n), Nodes0 % Z(n), &
        Nodes1 % X(n), Nodes1 % Y(n), Nodes1 % Z(n), ElemPres( n ), &
        STAT=istat ) 
    IF ( istat /= 0 ) CALL Fatal( 'CompressibilitySolve', 'Memory allocation error.' )
    AllocationsDone = .TRUE.
  END IF

  ElemPres = 0.0_dp
  ElemDisp = 0.0_dp

!------------------------------------------------------------------------------
! Initialize the system and do the assembly
!------------------------------------------------------------------------------
  CALL DefaultInitialize()
  
  DO t=1,Solver % NumberOfActiveElements
    Element => GetActiveElement(t)
    n = GetElementNOFNodes()
    NodeIndexes => Element % NodeIndexes
!------------------------------------------------------------------------------
    
    IF(PressureExists) THEN
      ElemPres(1:n) = PressureSolValues( PressureSolPerm(NodeIndexes(1:n)) )
    ELSE      
      ElemPres(1:n) = ReferencePressure
    END IF

    DO i = 1,n
      k = DisplacementSolPerm(NodeIndexes(i))
      DO j= 1,NSDOFs
        ElemDisp(j,i) = DisplacementSolValues( NSDOFs * (k-1) + j)
      END DO
    END DO

!------------------------------------------------------------------------------
!     Get element local matrix and rhs vector
!------------------------------------------------------------------------------
     CALL LocalMatrix(  STIFF, FORCE, Element, ElemDisp, ElemPres, n )
!------------------------------------------------------------------------------
!     Update global matrix and rhs vector from local matrix & vector
!------------------------------------------------------------------------------
     CALL DefaultUpdateEquations( STIFF, FORCE )

   END DO! <- elements

!------------------------------------------------------------------------------
   CALL DefaultFinishAssembly()
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!  Solve the system:
!------------------------------------------------------------------------------
   Norm = DefaultSolve()

!------------------------------------------------------------------------------
!  Eliminate negative entries if needed
!------------------------------------------------------------------------------

   IF(ListGetLogical(Solver % Values,'Eliminate Negative',Found) ) THEN
     CompressibilityFunction = MAX(0.0_dp, CompressibilityFunction )
   END IF


CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, Element, ElemDisp, ElemPres, n )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), ElemDisp(:,:), ElemPres(:) 
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(1,1,1),  PresAtIp
    REAL(KIND=dp) :: DetJ1, DetJ0, U, V, W, S, x, dVolume, Volume0, Volume1
    LOGICAL :: Stat
    INTEGER :: t, p, q, dim, NBasis, k
    TYPE(GaussIntegrationPoints_t) :: IP
!------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()
    STIFF  = 0.0_dp
    FORCE  = 0.0_dp
!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------

    NBasis = n
    IP = GaussPoints( Element )

    IF(DisplacedShape) THEN
      Nodes1 % x(1:n) = Model % Nodes % x(Element % NodeIndexes) - ElemDisp(1,1:n)
      Nodes1 % y(1:n) = Model % Nodes % y(Element % NodeIndexes) - ElemDisp(2,1:n)
      IF(dim == 3) Nodes1 % z(1:n) = Model % Nodes % z(Element % NodeIndexes) - ElemDisp(3,1:n)
    ELSE      
      Nodes1 % x(1:n) = Model % Nodes % x(Element % NodeIndexes) + ElemDisp(1,1:n)
      Nodes1 % y(1:n) = Model % Nodes % y(Element % NodeIndexes) + ElemDisp(2,1:n)
      IF(dim == 3) Nodes1 % z(1:n) = Model % Nodes % z(Element % NodeIndexes) + ElemDisp(3,1:n)
    END IF

    Nodes0 % x(1:n) = Model % Nodes % x(Element % NodeIndexes)
    Nodes0 % y(1:n) = Model % Nodes % y(Element % NodeIndexes)
    IF(dim == 3) Nodes0 % z(1:n) = Model % Nodes % z(Element % NodeIndexes)

    DO t=1,IP % n
      U = IP % u(t)
      V = IP % v(t)
      W = IP % w(t)
      S = IP % s(t)
            
!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------

      stat = ElementInfo( Element, Nodes1, U, V, W, DetJ1, Basis, dBasisdx, ddBasisddx, .FALSE. )
      stat = ElementInfo( Element, Nodes0, U, V, W, DetJ0, Basis, dBasisdx, ddBasisddx, .FALSE. )

      Volume0 = DetJ0
      Volume1 = DetJ1
      S = S * DetJ0

      IF ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
          CurrentCoordinateSystem() == CylindricSymmetric ) THEN

        x = SUM( Basis(1:n) * Nodes1 % x(1:n) )
        Volume1 = Volume1 * 2 * PI * x

        x = SUM( Basis(1:n) * Nodes0 % x(1:n) )
        Volume0 = Volume0 * 2 * PI * x
       
        S = S * x       
      END IF
            
      IF(DisplacedShape) THEN
        dVolume = Volume0 - Volume1
      ELSE
        dVolume = Volume1 - Volume0
      END IF
     

!------------------------------------------------------------------------------
!      Load at the integration point
!------------------------------------------------------------------------------
      PresAtIp = SUM( Basis(1:n) * ElemPres(1:n) )

!------------------------------------------------------------------------------
!      Finally, the elemental matrix & vector
!------------------------------------------------------------------------------       

      DO p=1,NBasis
        DO q=1,NBasis             
          STIFF(p,q) = STIFF(p,q) + s * PresAtIp * Basis(p) * Basis(q)
        END DO
      END DO
      
      DO p = 1, NBasis
        FORCE(p) = FORCE(p) + s * Basis(p) * dVolume / Volume0
      END DO

    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE CompressibilitySolver
!------------------------------------------------------------------------------
