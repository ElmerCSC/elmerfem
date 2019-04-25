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
! *  Authors: Peter Råback
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 14.02.2008
! *
! *****************************************************************************/




!------------------------------------------------------------------------------
!>  Subroutine for computing vorticity of vector fields. May be used to compute 
!>  either the whole vorticity vector (3D) or just its z-component (2D).
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE VorticitySolver( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE CoordinateSystems
  USE DefUtils

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
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName, CondName
  INTEGER :: i,j,dim,DOFs,ActiveDir,VorDofs
  LOGICAL :: ConstantBulkMatrix, ConstantBulkMatrixInUse, CSymmetry
  LOGICAL :: GotIt, GotCoeff, Visited = .FALSE., DirMask(3)
  REAL(KIND=dp) :: Unorm, Totnorm
  REAL(KIND=dp), POINTER CONTIG :: SaveRHS(:)
  REAL(KIND=dp), POINTER CONTIG :: ForceVectors(:,:), ForceVector(:)
#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: at0,at1,at2
#else
  REAL(KIND=dp) :: at0,at1,at2,CPUTime,RealTime
#endif
  TYPE(Variable_t), POINTER :: VorticitySol
  
  SAVE Visited

  CALL Info( 'VorticitySolver', '-------------------------------------',Level=4 )
  CALL Info( 'VorticitySolver', 'Computing the vorticity field',Level=4 )
  CALL Info( 'VorticitySolver', '-------------------------------------',Level=4 )

  dim = CoordinateSystemDimension()
!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
  IF ( COUNT( Solver % Variable % Perm > 0 ) <= 0 ) RETURN
  
  SolverParams => GetSolverParams()

  IF( dim == 2 ) THEN
    VorticitySol => Solver % Variable
    Dofs = VorticitySol % DOFs
    IF(Dofs /= 1) THEN
      PRINT *,'DOFS',Dofs
      CALL Fatal('VorticitySolver','Vorticity should have 1 component in 2D')
    END IF
  ELSE IF( dim == 3) THEN    
    VarName = GetString(SolverParams,'Vorticity Result Variable',GotIt )
    IF(.NOT. gotIt) VarName = 'Vorticity'
    VorticitySol => VariableGet( Solver % Mesh % Variables,  VarName )
    IF( ASSOCIATED(VorticitySol) ) THEN
      Dofs = VorticitySol % DOFs
      IF(Dofs /= 3) THEN
        PRINT *,'DOFS',Dofs
        CALL Fatal('VorticitySolver','Vorticity should have 3 components in 3D')
      END IF
    ELSE
      CALL Fatal('VorticitySolver','Vorticity Result Variable is missing: '//TRIM(VarName))      
    END IF
  END IF

  CSymmetry = CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric
  
  VarName = GetString(SolverParams,'Vorticity Variable',GotIt )
  IF(.NOT. GotIt) VarName = GetString(SolverParams,'Target Variable',GotIt )
  IF(.NOT. GotIt) VarName = TRIM('Velocity')

  ! For future use
  CondName = ListGetString(SolverParams,'Vorticity Coefficient',GotCoeff )
  
  at0 = RealTime()
  
  ConstantBulkMatrix = GetLogical( SolverParams, 'Constant Bulk Matrix', GotIt )
  ConstantBulkMatrixInUse = ConstantBulkMatrix .AND. &
      ASSOCIATED(Solver % Matrix % BulkValues)
  
  IF ( ConstantBulkMatrixInUse ) THEN
    Solver % Matrix % Values = Solver % Matrix % BulkValues        
    Solver % Matrix % RHS = 0.0_dp
  ELSE
    CALL DefaultInitialize()
  END IF

  ! If vorticity has many components, compute them one-by-one
  IF(Dofs > 1) THEN
    ALLOCATE(ForceVectors(SIZE(Solver % Matrix % RHS),Dofs-1))  
    ForceVectors = 0.0_dp
    SaveRHS => Solver % Matrix % RHS
  END IF

  CALL BulkAssembly()
  IF(.NOT. ConstantBulkMatrixInUse ) THEN
    CALL DefaultFinishBulkAssembly()
  END IF

  ! No flux BCs
  CALL DefaultFinishAssembly()

  at1 = RealTime()
  WRITE(Message,* ) 'Assembly Time: ',at1-at0
  CALL Info( 'VorticitySolver', Message, Level=5 )
        
!------------------------------------------------------------------------------     

  IF(Dofs > 1) THEN
    TotNorm = 0._dp
    DO i=1,Dofs
      IF(i==1) THEN
        Solver % Matrix % RHS => SaveRHS
      ELSE
        Solver % Matrix % RHS => ForceVectors(:,i-1)
      END IF
      UNorm = DefaultSolve()
      TotNorm = TotNorm + Unorm ** 2
      DO j=1,Solver % Matrix % NumberOfRows
        VorticitySol % Values(DOFs*(j-1)+i) = Solver % Variable % Values(j)
      END DO
    END DO
    Solver % Matrix % RHS => SaveRHS
    TotNorm = SQRT(TotNorm)
 
    DEALLOCATE( ForceVectors )
    Solver % Matrix % RHS => SaveRHS
    Solver % Variable % Norm = Totnorm
  ELSE
    TotNorm = DefaultSolve()
  END IF
!------------------------------------------------------------------------------     


  at2 = RealTime()
  WRITE(Message,* ) 'Solution Time: ',at2-at1
  CALL Info( 'VorticitySolver', Message, Level=5 )
  
  WRITE( Message, * ) 'Result Norm: ',TotNorm
  CALL Info( 'VorticitySolver', Message, Level=4 )
  
CONTAINS


!------------------------------------------------------------------------------
  SUBROUTINE BulkAssembly()
!------------------------------------------------------------------------------
       
    INTEGER :: elem,t,i,j,p,q,n,nd, Rank
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:,:)
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: weight,grad(3),C(3,3),detJ
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
    REAL(KIND=dp), ALLOCATABLE :: Vx(:), Vy(:), Vz(:), GradVx(:), GradVy(:), GradVz(:), Coeff(:)
    LOGICAL :: Found
    TYPE(ValueList_t), POINTER :: Material
    
    SAVE Coeff, Nodes
    
    n = MAX( Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes )
    ALLOCATE( STIFF(n,n), FORCE(dim,n), Coeff(n) )
    ALLOCATE( Vx(n), Vy(n), Vz(n), GradVx(3), GradVy(3), GradVz(3), &
        Basis(n), dBasisdx(n,3) )

    DO elem = 1,Solver % NumberOFActiveElements
         
      ! Element information
      ! ---------------------
      Element => GetActiveElement(elem)
      CALL GetElementNodes( Nodes )
      nd = GetElementNOFDOFs()
      n  = GetElementNOFNodes()

!      IF(GotCoeff) Coeff(1:n) = GetReal( GetMaterial(), Coeff, CondName, Found )

      CALL GetScalarLocalSolution( Vx, ComponentName(VarName,1) )
      CALL GetScalarLocalSolution( Vy, ComponentName(VarName,2) )
      CALL GetScalarLocalSolution( Vz, ComponentName(VarName,3) )

      ! Integrate local stresses:
      ! -------------------------
      IntegStuff = GaussPoints( Element )
      STIFF  = 0.0_dp
      FORCE  = 0.0_dp

      C = 0.0_dp
      DO i=1,dim
        C(i,i) = 1.0_dp
      END DO
      
      DO t=1,IntegStuff % n
        Found = ElementInfo( Element, Nodes, IntegStuff % u(t), &
                IntegStuff % v(t), IntegStuff % w(t), detJ, Basis, dBasisdx )
        
        Weight = IntegStuff % s(t) * detJ
        IF ( CSymmetry ) Weight = Weight * SUM( Basis(1:n) * Nodes % x(1:n) )
        
        IF ( .NOT. ConstantBulkMatrixInUse ) THEN
          DO p=1,nd
            DO q=1,nd
              STIFF(p,q) = STIFF(p,q) + Weight * Basis(q) * Basis(p)
            END DO
          END DO
        END IF
        
        GradVx(1:dim) = MATMUL( Vx(1:nd), dBasisdx(1:nd,1:dim) )
        GradVy(1:dim) = MATMUL( Vy(1:nd), dBasisdx(1:nd,1:dim) )
        GradVz(1:dim) = MATMUL( Vz(1:nd), dBasisdx(1:nd,1:dim) )

        IF(Dofs == 1) THEN
          FORCE(1,1:nd) = FORCE(1,1:nd) + &
              Basis(1:nd) * Weight * ( GradVy(1) - GradVx(2) )
        ELSE 
          FORCE(1,1:nd) = FORCE(1,1:nd) + &
              Basis(1:nd) * Weight * ( GradVz(2) - GradVy(3) )
          FORCE(2,1:nd) = FORCE(2,1:nd) + &
              Basis(1:nd) * Weight * ( GradVx(3) - GradVz(1) )
          FORCE(3,1:nd) = FORCE(3,1:nd) + &
              Basis(1:nd) * Weight * ( GradVy(1) - GradVx(2) )
        END IF
      END DO
      
!------------------------------------------------------------------------------
!      Update global matrices from local matrices 
!------------------------------------------------------------------------------

      IF ( .NOT. ConstantBulkMatrixInUse ) THEN
        IF(Dofs > 1) Solver % Matrix % RHS => SaveRHS
        CALL DefaultUpdateEquations( STIFF, FORCE(1,1:nd) )
      ELSE
        CALL DefaultUpdateForce( FORCE(1,1:nd) )        
      END IF

      ! Assembly the 2nd and 3rd r.h.s. in 3D case
      DO i=2,Dofs
        Solver % Matrix % RHS => ForceVectors(:,i-1)
        CALL DefaultUpdateForce( FORCE(i,1:nd) )
      END DO
    END DO
    
    DEALLOCATE( STIFF, FORCE, Basis, dBasisdx, Coeff, Vx, Vy, Vz, &
        GradVx, GradVy, GradVz )
!------------------------------------------------------------------------------
  END SUBROUTINE BulkAssembly
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE VorticitySolver
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Initialization for the primary solver: VorticitySolver.
!> \ingroup Solvers
!------------------------------------------------------------------------------
  SUBROUTINE VorticitySolver_Init( Model,Solver,dt,Transient )
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
    CHARACTER(LEN=MAX_NAME_LEN) :: VarName, VorName, HelpName
    LOGICAL :: GotIt
!------------------------------------------------------------------------------
    SolverParams => GetSolverParams()
    dim = CoordinateSystemDimension()

    VarName = GetString(SolverParams,'Vorticity Variable',GotIt )
    IF(gotIt) THEN
      VorName = 'Curl '//TRIM(VarName)
    ELSE
      VorName = 'Vorticity'
    END IF

    IF ( .NOT. ListCheckPresent( SolverParams,'Variable') ) THEN
      CALL ListAddInteger( SolverParams, 'Variable DOFs', 1 )
      IF(dim == 2) THEN
        CALL ListAddString( SolverParams, 'Variable', VorName )
      ELSE IF(dim == 3) THEN
        HelpName = 'F['//TRIM(Vorname)//':3]'
        CALL ListAddString(SolverParams, 'Exported Variable 1',HelpName)
        CALL ListAddString(SolverParams, 'Vorticity Result Variable','F')
        CALL ListAddString( SolverParams, 'Variable','-nooutput vorticity_temp' )
      ELSE
        CALL Fatal('VortictySolver_init','Vorticity computation makes sense only in 2D and 3D')
      END IF
    END IF

    IF( GetLogical( SolverParams,'Calculate Abs',GotIt) ) THEN
      CALL ListAddString( SolverParams,&
          NextFreeKeyword('Exported Variable',SolverParams),TRIM(VorName)//'_abs')        
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
  END SUBROUTINE VorticitySolver_Init
!------------------------------------------------------------------------------

