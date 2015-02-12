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
! *  Authors: Peter Råback
! *  Email:   Peter.Rabackcsc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 12.02.2002
! * 
! ******************************************************************************

!------------------------------------------------------------------------------
!> Subroutine that may be used to normalize the amplitude of artificial compressibility for
!> optimal fluid-structure coupling.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE CompressibilityScale( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver   !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model     !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt        !< Timestep size for time dependent simulations
  LOGICAL :: TransientSimulation !< Steady state or transient simulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  REAL(KIND=dp), ALLOCATABLE :: Pressure(:), &
      Displacement(:,:), Compressibility(:)
  TYPE(Solver_t), POINTER :: FlowSolver
  REAL(KIND=dp) :: SideVolume, SidePressure, SideArea, TotalVolume, &
      InitVolume, TotalVolumeCompress, CompressSuggest, &
      CompressScale, CompressScaleOld, Relax, Norm, &
      Time, PrevTime=0.0d0, MinimumSideVolume=1.0d10, &
      TransitionVolume, dVolume
  LOGICAL :: Stat, GotIt, SubroutineVisited=.False., &
      ScaleCompressibility, FsiBC
  INTEGER :: i,j,k,n,pn,t,DIM,mat_id,TimeStepVisited
  INTEGER, POINTER :: NodeIndexes(:)
  TYPE(Variable_t), POINTER :: Var,Dvar
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Material
  TYPE(Nodes_t) :: ElementNodes
  TYPE(Element_t), POINTER   :: CurrentElement
  CHARACTER(LEN=MAX_NAME_LEN) :: EquationName
!------------------------------------------------------------------------------

  SAVE SubroutineVisited, TimeStepVisited, PrevTime, InitVolume, &
      MinimumSideVolume

  Time = Solver % DoneTime
  IF(ABS(Time-PrevTime) > 1.0d-20) THEN 
    TimeStepVisited = 0
    PrevTime = Time
  END IF

  DO i=1,Model % NumberOfSolvers
    FlowSolver => Model % Solvers(i)
    IF ( ListGetString( FlowSolver % Values, 'Equation' ) == 'navier-stokes' ) EXIT
  END DO

  Mesh => Model % Meshes
  DO WHILE( ASSOCIATED(Mesh) )
    IF ( Mesh % OutputActive ) EXIT 
    Mesh => Mesh % Next
  END DO


  CALL SetCurrentMesh( Model, Mesh )
  Var => VariableGet( Mesh % Variables, 'Flow Solution', .TRUE. )
  Dvar => VariableGet( Mesh % Variables, 'Displacement', .TRUE.)

  ALLOCATE( ElementNodes % x(Mesh % MaxElementNodes) )
  ALLOCATE( ElementNodes % y(Mesh % MaxElementNodes) )
  ALLOCATE( ElementNodes % z(Mesh % MaxElementNodes) )
  ALLOCATE( Compressibility(   Mesh % MaxElementNodes ) )
  ALLOCATE( Pressure(   Mesh % MaxElementNodes ) )
  ALLOCATE( Displacement( 3,Mesh % MaxElementNodes ) )

  TotalVolume = 0.0d0
  TotalVolumeCompress = 0.0d0
  SidePressure = 0.0d0
  SideVolume = 0.0d0
  SideArea = 0.0d0
  
  DIM = CoordinateSystemDimension()
  EquationName = ListGetString( Solver % Values, 'Equation' )


  DO t=1,Solver % Mesh % NumberOfBulkElements

    CurrentElement => Solver % Mesh % Elements(t)
    
    IF ( .NOT. CheckElementEquation( Model, CurrentElement, EquationName ) &
        .AND. .NOT. CheckElementEquation( Model, CurrentElement, 'navier-stokes' ) ) &
        CYCLE
    
    n = CurrentElement % TYPE % NumberOfNodes
    NodeIndexes => CurrentElement % NodeIndexes

    ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
    ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
    ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)
    
    mat_id = ListGetInteger( Model % Bodies(CurrentElement % BodyId) % &
        Values, 'Material',minv=1,maxv=Model % NumberOfMaterials )
    
    Material => Model % Materials(mat_id) % Values

    Compressibility(1:n) = &
        ListGetReal(Material,'Artificial Compressibility',n,NodeIndexes,gotIt)
    IF(.NOT. gotIt) Compressibility(1:n) = 0.0d0

    CALL CompressibilityIntegrate(CurrentElement, n, ElementNodes, &
        Compressibility, TotalVolume, TotalVolumeCompress)
  END DO

  ! Compute the force acting on the boundary
  DO t = Mesh % NumberOfBulkElements + 1, &
            Mesh % NumberOfBulkElements + &
               Mesh % NumberOfBoundaryElements

!------------------------------------------------------------------------------
     CurrentElement => Mesh % Elements(t)
     IF ( CurrentElement % TYPE % ElementCode == 101 ) CYCLE

!------------------------------------------------------------------------------
!    Set the current element pointer in the model structure to 
!    reflect the element being processed
!------------------------------------------------------------------------------
     Model % CurrentElement => Mesh % Elements(t)
!------------------------------------------------------------------------------
     n = CurrentElement % TYPE % NumberOfNodes
     NodeIndexes => CurrentElement % NodeIndexes

     DO k=1, Model % NumberOfBCs
        IF ( Model % BCs(k) % Tag /= CurrentElement % BoundaryInfo % Constraint ) CYCLE

        FsiBC = ListGetLogical(Model % BCs(k) % Values,'Force BC',stat ) .OR. &
                ListGetLogical(Model % BCs(k) % Values,'Fsi BC',stat ) 
        IF(.NOT. FsiBC ) CYCLE

        ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
        ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
        ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)

        Pressure(1:n) = Var % Values(Var % DOFs * Var % Perm(NodeIndexes))

        Displacement = 0.0d0
        DO i=1,n
          DO j=1,DIM
            Displacement(j,i) = &
                DVar % Values(DVar % DOFs * (DVar % Perm(NodeIndexes(i))-1)+j)
           END DO
         END DO

        CALL PressureIntegrate(CurrentElement, n, ElementNodes)
     END DO
  END DO

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
    dVolume = TotalVolume-InitVolume
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
  CALL ListAddConstReal( Model % Simulation, &
      'res: Iterations',1.0d0*TimeStepVisited)

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
    
    INTEGER :: N_Integ, CoordSys
    
    LOGICAL :: stat
    INTEGER :: i,t
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    
    Ident = 0.0d0
    DO i=1,3
      Ident(i,i) = 1.0d0
    END DO
    
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
    IntegStuff = GaussPoints( Element )
    
    CoordSys = CurrentCoordinateSystem()
    

!------------------------------------------------------------------------------
    DO t=1,IntegStuff % n
!------------------------------------------------------------------------------
       u = IntegStuff % u(t)
       v = IntegStuff % v(t)
       w = IntegStuff % w(t)
!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, u, v, w, &
           SqrtElementMetric, Basis, dBasisdx )
       
       s = SqrtElementMetric * IntegStuff % s(t)

       IF ( CoordSys /= Cartesian ) THEN
         X = SUM( Nodes % X(1:n) * Basis(1:n) )
         Y = SUM( Nodes % Y(1:n) * Basis(1:n) )
         Z = SUM( Nodes % Z(1:n) * Basis(1:n) )
         CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,X,Y,Z )
         s = s * SqrtMetric
       END IF
       
       Normal = Normalvector( Element,Nodes, u,v, .TRUE. )
       
       NormalDisplacement = 0.0
       DO i=1,DIM
         NormalDisplacement = NormalDisplacement + &
             SUM(Basis(1:n) * Displacement(i,1:n)) * Normal(i)
       END DO
       
       Pres = SUM( Basis(1:n) * Pressure(1:n))

       SideVolume = SideVolume + s * ABS(NormalDisplacement)
       SidePressure = SidePressure + s * ABS(Pres)
       SideArea = SideArea + s
        
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
    INTEGER :: i,p,q,t,DIM, NBasis, CoordSys
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
    REAL(KIND=dp) :: X,Y,Z,Metric(3,3),SqrtMetric,Symb(3,3,3),dSymb(3,3,3,3)

!------------------------------------------------------------------------------
    DIM = CoordinateSystemDimension()
    CoordSys = CurrentCoordinateSystem()

    Metric = 0.0d0
    Metric(1,1) = 1.0d0
    Metric(2,2) = 1.0d0
    Metric(3,3) = 1.0d0

!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------

    NBasis = n
    IntegStuff = GaussPoints( Element )

!------------------------------------------------------------------------------
    DO t=1,IntegStuff % n
      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)

!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, U, V, W, SqrtElementMetric, &
          Basis, dBasisdx )
      s = IntegStuff % s(t) * SqrtElementMetric

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
!> Subroutine for computing the artificial compressibility from the volume change 
!> of elements. The volume change is obtained by extending the displacement field
!> of a test load to the fluid domain.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE CompressibilitySolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver   !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model     !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt        !< Timestep size for time dependent simulations
  LOGICAL :: TransientSimulation !< Steady state or transient simulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Matrix_t), POINTER  :: StiffMatrix
  TYPE(Variable_t), POINTER :: DisplacementSol, PressureSol
  TYPE(Element_t), POINTER :: CurrentElement

  CHARACTER(LEN=MAX_NAME_LEN) :: DisplacementVariableName, PressureVariableName

  REAL(KIND=dp), POINTER :: CompressibilityFunction(:), DisplacementSolValues(:), &
      PressureSolValues(:)
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:), ElemDisplacement(:,:), ElemPressure(:)
  REAL(KIND=dp) :: Norm, ReferencePressure

  INTEGER, POINTER :: NodeIndexes(:), Perm(:), DisplacementSolPerm(:), PressureSolPerm(:)
  INTEGER :: k, t, i, j, n, istat, NSDOFs, DIM

  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Nodes_t) :: Nodes0, Nodes1
 
  LOGICAL :: AllocationsDone = .FALSE., Found, DisplacedShape, PressureExists

  SAVE STIFF, FORCE, Nodes0, Nodes1, ElemPressure, ElemDisplacement, AllocationsDone

  CALL Info('CompressibilitySolver',' ')
  CALL Info('CompressibilitySolver','----------------------------')
  CALL Info('CompressibilitySolver','Compressibility Solver')
  CALL Info('CompressibilitySolver','----------------------------')
  CALL Info('CompressibilitySolver',' ')

!------------------------------------------------------------------------------
! Get variables needed for solution
!------------------------------------------------------------------------------
  IF ( .NOT.ASSOCIATED( Solver % Matrix ) ) RETURN

  CompressibilityFunction => Solver % Variable % Values
  Perm => Solver % Variable % Perm
  
  IF(TransientSimulation) THEN
    CALL Warn('CompressibilitySolver','Implemented only for steady state')
    PRINT *,'AC interval',MINVAL(CompressibilityFunction),MAXVAL(CompressibilityFunction)
    RETURN
  END IF

  DIM = CoordinateSystemDimension()

  StiffMatrix => Solver % Matrix
!------------------------------------------------------------------------------
! Get initial values ( Default is DisplacementSolvers 'Displacement' )
!------------------------------------------------------------------------------
  SolverParams => GetSolverParams()
  DisplacementVariableName = GetString( SolverParams, &
     'Displacement Variable Name', Found )
  IF ( .NOT. Found ) THEN
    DisplacementVariableName = 'Displacement'
  END IF

  DisplacementSol => VariableGet( Solver % Mesh % Variables, DisplacementVariableName )
  IF ( ASSOCIATED( DisplacementSol ) ) THEN
    DisplacementSolPerm => DisplacementSol % Perm
    DisplacementSolValues => DisplacementSol % Values
    NSDOFs = DisplacementSol % DOFs
  ELSE
    CALL Warn( 'CompressibilitySolver', 'No variable for displacement associated.' )
    CALL Warn( 'CompressibilitySolver', 'Quitting execution of CompressibilitySolver.' ) 
    RETURN
  END IF

  SolverParams => GetSolverParams()
  PressureVariableName = GetString( SolverParams, &
      'Pressure Variable Name', Found )
  IF ( .NOT. Found ) THEN
    PressureVariableName = 'Pressure'
  END IF

  PressureSol => VariableGet( Solver % Mesh % Variables, PressureVariableName )
  IF ( ASSOCIATED( PressureSol ) ) THEN
    PressureExists = .TRUE.
    PressureSolPerm => PressureSol % Perm
    PressureSolValues => PressureSol % Values
  ELSE
    PressureExists = .FALSE.
    CALL Info( 'CompressibilitySolver', 'No variable for pressure associated.' )
    ReferencePressure = ListGetConstReal( Solver % Values,'Reference Pressure',Found)
    IF(.NOT. Found) THEN
      CALL Warn( 'CompressibilitySolver', 'No pressure defined' )
      CALL Warn( 'CompressibilitySolver', 'Quitting execution of CompressibilitySolver.' ) 
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
     N = Solver % Mesh % MaxElementNodes ! just big enough for elemental arrays
     ALLOCATE( FORCE( N ), STIFF(N,N), ElemDisplacement(3,N), &
         Nodes0 % X(N), Nodes0 % Y(N), Nodes0 % Z(N), &
         Nodes1 % X(N), Nodes1 % Y(N), Nodes1 % Z(N), ElemPressure( N ), &
         STAT=istat ) 

     IF ( istat /= 0 ) CALL Fatal( 'CompressibilitySolve', 'Memory allocation error.' )
     AllocationsDone = .TRUE.
  END IF

  ElemPressure = 0.0d0
  ElemDisplacement = 0.0d0

!------------------------------------------------------------------------------
! Initialize the system and do the assembly
!------------------------------------------------------------------------------
  CALL DefaultInitialize()
  
  DO t=1,Solver % NumberOfActiveElements
    CurrentElement => GetActiveElement(t)
    n = GetElementNOFNodes()
    NodeIndexes => CurrentElement % NodeIndexes
!------------------------------------------------------------------------------
    
    IF(PressureExists) THEN
      ElemPressure(1:n) = PressureSolValues( PressureSolPerm(NodeIndexes(1:n)) )
    ELSE      
      ElemPressure = ReferencePressure
    END IF

    DO i = 1,n
      k = DisplacementSolPerm(NodeIndexes(i))
      DO j= 1,NSDOFs
        ElemDisplacement(j,i) = DisplacementSolValues( NSDOFs * (k-1) + j)
      END DO
    END DO

!------------------------------------------------------------------------------
!     Get element local matrix and rhs vector
!------------------------------------------------------------------------------
     CALL LocalMatrix(  STIFF, FORCE, CurrentElement, ElemDisplacement, ElemPressure, n )
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
     DO i = 1, SIZE(Perm)  
       j = Perm(i)
       IF(j == 0) CYCLE
       CompressibilityFunction(j) = MAX(0.0d0, CompressibilityFunction(j) )
     END DO
   END IF


CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, Element, ElemDisp, ElemPres, n )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), ElemDisp(:,:), ElemPres(:) 
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: StokesCompressibility
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(1,1,1),  PresAtIp
    REAL(KIND=dp) :: DetJ1, DetJ0, U, V, W, S, x, dVolume, Volume0, Volume1
    LOGICAL :: Stat
    INTEGER :: t, p, q, dim, NBasis, k
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
!------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()
    STIFF  = 0.0d0
    FORCE  = 0.0d0
!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------

    NBasis = n
    IntegStuff = GaussPoints( Element )

    IF(DisplacedShape) THEN
      Nodes1 % x(1:n) = Model % Nodes % x(Element % NodeIndexes) - ElemDisp(1,1:n)
      Nodes1 % y(1:n) = Model % Nodes % y(Element % NodeIndexes) - ElemDisp(2,1:n)
      IF(DIM == 3) Nodes1 % z(1:n) = Model % Nodes % z(Element % NodeIndexes) - ElemDisp(3,1:n)
    ELSE      
      Nodes1 % x(1:n) = Model % Nodes % x(Element % NodeIndexes) + ElemDisp(1,1:n)
      Nodes1 % y(1:n) = Model % Nodes % y(Element % NodeIndexes) + ElemDisp(2,1:n)
      IF(DIM == 3) Nodes1 % z(1:n) = Model % Nodes % z(Element % NodeIndexes) + ElemDisp(3,1:n)
    END IF

    Nodes0 % x(1:n) = Model % Nodes % x(Element % NodeIndexes)
    Nodes0 % y(1:n) = Model % Nodes % y(Element % NodeIndexes)
    IF(DIM == 3) Nodes0 % z(1:n) = Model % Nodes % z(Element % NodeIndexes)

    DO t=1,IntegStuff % n
      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)
      S = IntegStuff % s(t)
            
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
          STIFF(p,q) = STIFF(p,q) &
              + s * PresAtIp * Basis(p) * Basis(q)
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
