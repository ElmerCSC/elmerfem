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
! *
! ******************************************************************************
! *
! *  Authors: Peter Råback
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 13.04.2011
! *
! *****************************************************************************/


!------------------------------------------------------------------------------
!> Initilization for the primary solver: ShearrateSolver
!> \ingroup Solvers
!------------------------------------------------------------------------------
  SUBROUTINE ShearrateSolver_Init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE

    TYPE(Model_t)  :: Model
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: DT
    LOGICAL :: Transient
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: SolverParams
    INTEGER :: dim
    CHARACTER(LEN=MAX_NAME_LEN) :: VarName, FieldName
    LOGICAL :: GotIt
!------------------------------------------------------------------------------
    SolverParams => GetSolverParams()
    dim = CoordinateSystemDimension()

    VarName = GetString(SolverParams,'Target Variable',GotIt )
    IF(.NOT. gotIt) VarName = TRIM('Velocity')

    IF ( .NOT. ListCheckPresent( SolverParams,'Variable') ) THEN
      FieldName = 'Shearrate '//TRIM(VarName)
      CALL ListAddString( SolverParams, 'Variable', FieldName )
    END IF

    IF( GetLogical( SolverParams,'Calculate Viscosity',GotIt ) ) THEN
      FieldName = 'Viscosity '//TRIM(VarName)
      CALL ListAddString( SolverParams,&	
          NextFreeKeyword('Exported Variable',SolverParams),FieldName)
    END IF

    CALL ListAddInteger( SolverParams, 'Time derivative order', 0 )

    ! Add linear system defaults: cg+ILU0
    IF(.NOT. ListCheckPresent(SolverParams,'Linear System Solver')) &
      CALL ListAddString(SolverParams,'Linear System Solver','Iterative')
    IF(.NOT. ListCheckPresent(SolverParams,'Linear System Iterative Method')) &
      CALL ListAddString(SolverParams,'Linear System Iterative Method','cg')
    IF(.NOT. ListCheckPresent(SolverParams,'Linear System Preconditioning')) &
      CALL ListAddString(SolverParams,'Linear System Preconditioning','ILU0')
    IF(.NOT. ListCheckPresent(SolverParams,'Linear System Max Iterations')) &
      CALL ListAddInteger(SolverParams,'Linear System Max Iterations',500)
    IF(.NOT. ListCheckPresent(SolverParams,'Linear System Convergence Tolerance')) &
      CALL ListAddConstReal(SolverParams,'Linear System Convergence Tolerance',1.0e-8_dp)

!------------------------------------------------------------------------------
  END SUBROUTINE ShearrateSolver_Init
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Subroutine for computing the shearrate of a given velocity field
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE ShearrateSolver( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE CoordinateSystems
  USE DefUtils
  USE Materialmodels

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL ::  Transient
  REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(ValueList_t),POINTER :: SolverParams
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName
  INTEGER :: i,j,dim,DOFs
  LOGICAL :: ConstantBulkMatrix, ConstantBulkMatrixInUse
  LOGICAL :: CalculateViscosity, GotIt
  REAL(KIND=dp) :: Unorm
  REAL(KIND=dp), POINTER CONTIG :: ForceVector(:), ViscVector(:), SaveRHS(:)
  REAL(KIND=dp), POINTER :: ShearrateField(:), ViscField(:)
  TYPE(Variable_t), POINTER :: ShearrateSol, ViscSol
  LOGICAL :: AllocationsDone = .FALSE. 

  SAVE ViscVector, AllocationsDone 
 
  CALL Info( 'ShearrateSolver', '-------------------------------------',Level=4 )
  CALL Info( 'ShearrateSolver','Computing the shearrate field',Level=4 )
  CALL Info( 'ShearrateSolver', '-------------------------------------',Level=4 )

  dim = CoordinateSystemDimension()
!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
  IF ( COUNT( Solver % Variable % Perm > 0 ) <= 0 ) RETURN
  
  SolverParams => GetSolverParams()

  ShearrateSol => Solver % Variable
  Dofs = ShearrateSol % DOFs
  IF(Dofs /= 1) THEN
    CALL Fatal('ShearrateSolver','Shearrate should have just 1 component')
  END IF
  
  VarName = GetString(SolverParams,'Target Variable',GotIt )
  IF(.NOT. gotIt) VarName = TRIM('Velocity')

  CalculateViscosity = GetLogical( SolverParams,'Calculate Viscosity',GotIt)
  IF( CalculateViscosity ) THEN
    IF( .NOT. AllocationsDone ) THEN 
      ALLOCATE(ViscVector(SIZE(Solver % Matrix % RHS)))  
      AllocationsDone = .TRUE.
    END IF
    ViscVector = 0.0_dp
    ViscSol => VariableGet( Solver % Mesh % Variables, 'Viscosity '//TRIM(VarName) )
    ViscField => ViscSol % Values
  END IF

  ConstantBulkMatrix = GetLogical( SolverParams, 'Constant Bulk Matrix', GotIt )
  ConstantBulkMatrixInUse = ConstantBulkMatrix .AND. &
      ASSOCIATED(Solver % Matrix % BulkValues)
  ForceVector => Solver % Matrix % rhs 
  ShearrateField => Solver % Variable % Values

  IF ( ConstantBulkMatrixInUse ) THEN
    Solver % Matrix % Values = Solver % Matrix % BulkValues        
    ForceVector = 0.0_dp
  ELSE
    CALL DefaultInitialize()
  END IF

  CALL BulkAssembly()
  IF(.NOT. ConstantBulkMatrixInUse ) THEN
    CALL DefaultFinishBulkAssembly()
  END IF

  ! No Flux BCs 
  CALL DefaultFinishAssembly()
        
!------------------------------------------------------------------------------     
  IF( CalculateViscosity ) THEN
    CALL Info( 'ShearrateSolver','Solving for viscosity',Level=5 )
    Solver % Matrix % RHS => ViscVector
    Solver % Variable % Values => ViscField
    UNorm = DefaultSolve()
    Solver % Variable % Values => ShearrateField
    Solver % Matrix % RHS => ForceVector
  END IF
  
  CALL Info( 'ShearrateSolver','Solving for shearrare',Level=5 )
  UNorm = DefaultSolve()

!------------------------------------------------------------------------------     
  
CONTAINS


!------------------------------------------------------------------------------
  SUBROUTINE BulkAssembly()
!------------------------------------------------------------------------------
       
    INTEGER :: elem,t,i,j,p,q,n,nd, Rank
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:), FORCE2(:)
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3)
    REAL(KIND=dp) :: weight,grad(3),detJ,ShearRate,Velo(3),dVelodx(3,3),s,x,y,z, &
        u,v,w,ViscAtIp, RhoAtIP, mu, muder0
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
    REAL(KIND=dp), ALLOCATABLE :: Vx(:), Vy(:), Vz(:), NodalRho(:), NodalVisc(:)
    LOGICAL :: Found
    TYPE(ValueList_t), POINTER :: Material
    LOGICAL :: AllocationsDone = .FALSE.

    SAVE Nodes, STIFF, FORCE, FORCE2, Vx, Vy, Vz, Basis, dBasisdx, &
        NodalRho, NodalVisc, AllocationsDone
    
    IF(.NOT. AllocationsDone ) THEN
      n = Solver % Mesh % MaxElementNodes 
      ALLOCATE( STIFF(n,n), FORCE(n), FORCE2(n), Vx(n), Vy(n), Vz(n), &
          Basis(n), NodalVisc(n), NodalRho(n), dBasisdx(n,3) )
      Vx = 0.0_dp
      Vy = 0.0_dp
      Vz = 0.0_dp
      AllocationsDone = .TRUE.
    END IF


    DO elem = 1,Solver % NumberOFActiveElements
         
      ! Element information
      ! ---------------------
      Element => GetActiveElement(elem)
      CALL GetElementNodes( Nodes )
      nd = GetElementNOFDOFs()
      n  = GetElementNOFNodes()

      IF( CalculateViscosity ) THEN
        Material => GetMaterial()
        NodalVisc(1:n) = GetReal( Material,'Viscosity')
        NodalRho(1:n) = GetReal( Material,'Density')
      END IF

      CALL GetScalarLocalSolution( Vx, ComponentName(VarName,1) )
      CALL GetScalarLocalSolution( Vy, ComponentName(VarName,2) )      
      IF( dim > 2 ) THEN
        CALL GetScalarLocalSolution( Vz, ComponentName(VarName,3) )
      END IF

      ! Integrate local stresses:
      ! -------------------------
      IntegStuff = GaussPoints( Element )
      STIFF  = 0.0_dp
      FORCE  = 0.0_dp
      FORCE2 = 0.0_dp


      DO t=1,IntegStuff % n

        u = IntegStuff % u(t)
        v = IntegStuff % v(t)
        w = IntegStuff % w(t)

        Found = ElementInfo( Element, Nodes, u, v, w, detJ, Basis, dBasisdx )

        !------------------------------------------------------------------------------
        !   Coordinate system dependent information
        !------------------------------------------------------------------------------
        x = SUM( Nodes % x(1:n) * Basis(1:n) )
        y = SUM( Nodes % y(1:n) * Basis(1:n) )
        z = SUM( Nodes % z(1:n) * Basis(1:n) )
        CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )

        s = detJ * SqrtMetric * IntegStuff % s(t)

        IF( CalculateViscosity ) THEN
          RhoAtIP  = SUM( Nodalrho(1:n) * Basis(1:n) )
          ViscAtIP = SUM( NodalVisc(1:n) * Basis(1:n) )
          IF( ListCheckPresent( Material, 'Viscosity Model' ) ) THEN
            mu = EffectiveViscosity( ViscAtIP, RhoAtIp, Vx, Vy, Vz, &
                Element, Nodes, n, n, u, v, w,  muder0, LocalIP=t )
          ELSE
            mu = ViscAtIP
          END IF
        END IF

        IF ( .NOT. ConstantBulkMatrixInUse ) THEN
          DO p=1,nd
            DO q=1,nd
              STIFF(p,q) = STIFF(p,q) + s * Basis(q) * Basis(p)
            END DO
          END DO
        END IF
        
        DO j=1,3
          dVelodx(1,j) = SUM( Vx(1:nd)*dBasisdx(1:nd,j) )
          dVelodx(2,j) = SUM( Vy(1:nd)*dBasisdx(1:nd,j) )
          dVelodx(3,j) = SUM( Vz(1:nd)*dBasisdx(1:nd,j) )
        END DO
        
        Velo(1) = SUM( Basis(1:nd) * Vx(1:nd) )
        Velo(2) = SUM( Basis(1:nd) * Vy(1:nd) )
        Velo(3) = SUM( Basis(1:nd) * Vz(1:nd) )
        
        ShearRate = SQRT( SecondInvariant(Velo,dVelodx,Metric,Symb)/2 )

        FORCE(1:nd) = FORCE(1:nd) + Basis(1:nd) * s * ShearRate

        IF( CalculateViscosity ) THEN
          FORCE2(1:nd) = FORCE2(1:nd) + Basis(1:nd) * s * mu
        END IF
      END DO
      
!------------------------------------------------------------------------------
!      Update global matrices from local matrices 
!------------------------------------------------------------------------------

      IF ( .NOT. ConstantBulkMatrixInUse ) THEN
        CALL DefaultUpdateEquations( STIFF, FORCE(1:nd) )
      ELSE
        CALL DefaultUpdateForce( FORCE(1:nd) )        
      END IF

      IF( CalculateViscosity ) THEN
        Solver % Matrix % RHS => ViscVector
        CALL DefaultUpdateForce( FORCE2 )
        Solver % Matrix % RHS => ForceVector
      END IF


    END DO
    

!------------------------------------------------------------------------------
  END SUBROUTINE BulkAssembly
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
END SUBROUTINE ShearrateSolver
!------------------------------------------------------------------------------

