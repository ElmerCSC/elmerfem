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
!> Subroutine computes the divergence of vector fields using Galerkin method. 
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE DivergenceSolver( Model,Solver,dt,Transient )
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
  INTEGER :: i,j,dim,DOFs
  LOGICAL :: ConstantBulkMatrix, ConstantBulkMatrixInUse, CSymmetry
  LOGICAL :: GotIt, GotCoeff, Visited = .FALSE.
  REAL(KIND=dp) :: Norm
#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: at0,at1,at2
#else
  REAL(KIND=dp) :: at0,at1,at2,CPUTime,RealTime
#endif
  TYPE(Variable_t), POINTER :: DivergenceSol
  
  SAVE Visited
 
  CALL Info( 'DivergenceSolver', '-------------------------------------',Level=4 )
  CALL Info( 'DivergenceSolver','Computing the divergence field',Level=4 )
  CALL Info( 'DivergenceSolver', '-------------------------------------',Level=4 )

  dim = CoordinateSystemDimension()
!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------

  SolverParams => GetSolverParams()
  DivergenceSol => Solver % Variable

  IF( ASSOCIATED(DivergenceSol)) THEN
    IF ( COUNT( Solver % Variable % Perm > 0 ) <= 0 ) THEN
      CALL Warn('DivergenceSolver','Size is zero, nothing to compute')
      RETURN
    END IF
    IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) THEN
      CALL Warn('DivergenceSolver','No matrix associated, exiting')
      RETURN
    END IF
    Dofs = DivergenceSol % DOFs
    IF(Dofs /= 1) CALL Fatal('DivergenceSolver','Divergence should have 1 component in 2D')
  ELSE
     CALL Fatal('DivergenceSolver','Variable does not exist or its size is zero!')      
  END IF

  CSymmetry = CurrentCoordinateSystem() == AxisSymmetric .OR. &
       CurrentCoordinateSystem() == CylindricSymmetric

  VarName = GetString(SolverParams,'Divergence Variable',GotIt )
  IF(.NOT. GotIt) VarName = GetString(SolverParams,'Target Variable',GotIt )
  IF(.NOT. gotIt) VarName = TRIM('Velocity')

  ! For future use
  CondName = ListGetString(SolverParams,'Divergence Coefficient',GotCoeff )
  
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

  CALL BulkAssembly()
  IF( ConstantBulkMatrix ) THEN
    IF(.NOT. ConstantBulkMatrixInUse ) THEN
      CALL DefaultFinishBulkAssembly( BulkUpdate = .TRUE.)
    END IF
  END IF

  CALL DefaultFinishAssembly()
  CALL DefaultDirichletBCs()

  at1 = RealTime()
  WRITE(Message,* ) 'Assembly Time: ',at1-at0
  CALL Info( 'DivergenceSolver', Message, Level=5 )
        
!------------------------------------------------------------------------------     

  Norm = DefaultSolve()

!------------------------------------------------------------------------------     
  
  at2 = RealTime()
  WRITE(Message,* ) 'Solution Time: ',at2-at1
  CALL Info( 'DivergenceSolver', Message, Level=5 )
  
  WRITE( Message, * ) 'Result Norm: ',Norm
  CALL Info( 'DivergenceSolver', Message, Level=4 )
  
CONTAINS


!------------------------------------------------------------------------------
  SUBROUTINE BulkAssembly()
!------------------------------------------------------------------------------
       
    INTEGER :: elem,t,i,j,p,q,n,nd, Rank
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:)
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: weight,grad(3),C(3,3),detJ,Source,x
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
    REAL(KIND=dp), ALLOCATABLE :: Vx(:), Vy(:), Vz(:), Coeff(:)
    LOGICAL :: Found
    TYPE(ValueList_t), POINTER :: Material
    
    SAVE Coeff, Nodes
    
    n = MAX( Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes )
    ALLOCATE( STIFF(n,n), FORCE(n), Coeff(n) )
    ALLOCATE( Vx(n), Vy(n), Vz(n), Basis(n), dBasisdx(n,3) )

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
      IF(dim == 3) CALL GetScalarLocalSolution( Vz, ComponentName(VarName,3) )

      IntegStuff = GaussPoints( Element )
      STIFF  = 0.0_dp
      FORCE  = 0.0_dp

!      C = 0.0_dp
!      DO i=1,dim
!        C(i,i) = 1.0_dp
!      END DO
      
      DO t=1,IntegStuff % n
        Found = ElementInfo( Element, Nodes, IntegStuff % u(t), &
                IntegStuff % v(t), IntegStuff % w(t), detJ, Basis, dBasisdx )
        
        Weight = IntegStuff % s(t) * detJ
        IF ( CSymmetry ) THEN
          x = SUM( Basis(1:n) * Nodes % x(1:n) ) 
          Weight = Weight * x
        END IF
        
        IF ( .NOT. ConstantBulkMatrixInUse ) THEN
          DO p=1,nd
            DO q=1,nd
              STIFF(p,q) = STIFF(p,q) + Weight * Basis(q) * Basis(p)
            END DO
          END DO
        END IF
        
        Source = SUM( dBasisdx(1:nd,1) * Vx(1:nd) )
        Source = Source + SUM( dBasisdx(1:nd,2) * Vy(1:nd) )
        IF(DIM == 3) Source = Source + SUM( dBasisdx(1:nd,3) * Vz(1:nd) )

	IF( CSymmetry ) THEN
          Source = Source + SUM( Basis(1:nd) * Vx(1:nd) ) / x
        END IF

        FORCE(1:nd) = FORCE(1:nd) + Basis(1:nd) * Weight * Source
      END DO
      
!------------------------------------------------------------------------------
!      Update global matrices from local matrices 
!------------------------------------------------------------------------------

      IF ( .NOT. ConstantBulkMatrixInUse ) THEN
        CALL DefaultUpdateEquations( STIFF, FORCE(1:nd) )
      ELSE
        CALL DefaultUpdateForce( FORCE(1:nd) )        
      END IF
    END DO
    
    DEALLOCATE( STIFF, FORCE, Basis, dBasisdx, Coeff, Vx, Vy, Vz )

!------------------------------------------------------------------------------
  END SUBROUTINE BulkAssembly
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
END SUBROUTINE DivergenceSolver
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Initialization for the primary solver, i.e. DivergenceSolver.
!> \ingroup Solvers
!------------------------------------------------------------------------------
  SUBROUTINE DivergenceSolver_Init( Model,Solver,dt,Transient )
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
    CHARACTER(LEN=MAX_NAME_LEN) :: VarName, DivName
    LOGICAL :: GotIt
!------------------------------------------------------------------------------
    SolverParams => GetSolverParams()
    dim = CoordinateSystemDimension()

    VarName = GetString(SolverParams,'Divergence Variable',GotIt )
    IF(gotIt) THEN
      DivName = 'Div '//TRIM(VarName)
    ELSE
      DivName = 'Divergence'
    END IF

    IF ( .NOT. ListCheckPresent( SolverParams,'Variable') ) THEN
      CALL ListAddInteger( SolverParams, 'Variable DOFs', 1 )
      CALL ListAddString( SolverParams, 'Variable', DivName )
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
  END SUBROUTINE DivergenceSolver_Init
!------------------------------------------------------------------------------

