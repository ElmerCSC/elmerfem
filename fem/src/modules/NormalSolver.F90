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
! *  Authors: Peter Råback, Juha Ruokolainen
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 20.06.2007
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> Solver for computing average normals for nodes using the Galerkin method.
!> For real meshes the direction of the normal is discrete for curved surfaces
!> and this solution should give a smoother normal direction. 
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE NormalSolver( Model,Solver,dt,Transient )
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
  TYPE(ValueList_t), POINTER :: SolverParams
  CHARACTER(LEN=MAX_NAME_LEN) :: Vname, VarName, CondName
  INTEGER :: i,j,k,dim,DOFs
  LOGICAL :: ConstantBulkMatrix, ConstantBulkMatrixInUse, CSymmetry
  LOGICAL :: GotIt, Visited = .FALSE., SetD
  REAL(KIND=dp) :: Unorm, Totnorm, nrm
  REAL(KIND=dp), ALLOCATABLE, TARGET :: ForceVector(:,:)
  REAL(KIND=dp), POINTER  CONTIG :: SaveRHS(:)
#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: at0,at1,at2
#else
  REAL(KIND=dp) :: at0,at1,at2,CPUTime,RealTime
#endif
  TYPE(Variable_t), POINTER :: NrmSol

  REAL(KIND=dp), ALLOCATABLE :: Values(:)
  
  SAVE Visited

 
  CALL Info( 'NormalSolver', '-------------------------------------',Level=4 )
  CALL Info( 'NormalSolver', 'Computing the normals',Level=4 )
  CALL Info( 'NormalSolver', '-------------------------------------',Level=4 )

  dim = CoordinateSystemDimension()
!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
  IF ( COUNT( Solver % Variable % Perm > 0 ) <= 0 ) RETURN


  SolverParams => GetSolverParams()

  VarName = GetString(SolverParams,'Normals Result Variable',GotIt )
  IF(.NOT. GotIt) THEN
    CALL Fatal('NormalSolver','> Normals Result Variable < not found!')
  END IF

  NrmSol => VariableGet( Solver % Mesh % Variables,  VarName )
  IF(ASSOCIATED(NrmSol)) THEN
     Dofs = NrmSol % DOFs
     IF(Dofs /= DIM) THEN
       CALL Fatal('NormalSolver','The normals should have DOFs equal to DIM')
     END IF
  ELSE
    CALL DefaultVariableAdd( VarName, dim, Var = NrmSol ) 
  END IF
  

  CSymmetry = CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric
  
  at0 = RealTime()
  
  ConstantBulkMatrix = GetLogical( SolverParams, 'Constant Bulk Matrix', GotIt )
  ConstantBulkMatrixInUse = ConstantBulkMatrix .AND. &
      ASSOCIATED(Solver % Matrix % BulkValues)
  
  IF ( ConstantBulkMatrixInUse ) THEN
    Solver % Matrix % Values = Solver % Matrix % BulkValues        
  ELSE
    CALL DefaultInitialize()
  END IF

  ALLOCATE(ForceVector(SIZE(Solver % Matrix % RHS),dim))  
  ForceVector = 0.0_dp
  SaveRHS => Solver % Matrix % RHS
!    
  CALL BulkAssembly()
  CALL DefaultFinishBulkAssembly()

  !No flux BCs for this solver
  CALL DefaultFinishAssembly()

  at1 = RealTime()
  WRITE(Message,* ) 'Assembly Time: ',at1-at0
  CALL Info( 'NormalSolver', Message, Level=5 )
!        
!------------------------------------------------------------------------------     
  SetD = GetLogical(GetSolverParams(), 'Set Dirichlet Conditions',GotIt)
  IF ( SetD ) THEN
    ALLOCATE(Values(SIZE(Solver % Matrix % Values)))
    Vname = Solver % Variable % Name
    Values = Solver % Matrix % Values
  END IF
  DO i=1,dim
    Solver % Matrix % RHS => ForceVector(:,i)

    IF ( SetD ) THEN
      Solver % Variable % Name = TRIM(VarName) // ' ' // TRIM(I2S(i))
      Solver % Matrix % Values = Values
      CALL DefaultDirichletBCs()
      Solver % Variable % Name = Vname
    ELSE
      CALL DefaultDirichletBCs()
    END IF

    Solver % Variable % Values = 0._dp
    UNorm = DefaultSolve()

    DO j=1,Solver % Matrix % NumberOfRows
      NrmSol % Values(DOFs*(j-1)+i) = Solver % Variable % Values(j)
    END DO
    Solver % Matrix % RHS => SaveRHS
  END DO
  IF ( SetD ) DEALLOCATE(Values)

  DO i=1,Solver % Matrix % NumberOFRows
    nrm = 0.0_dp
    DO j=1,dim
      k = DOFs*(i-1)+j
      nrm = nrm + NrmSol % Values(k)**2
    END DO
    IF ( nrm > 0 ) THEN
      nrm = 1.0_dp / SQRT(nrm)
      DO j=1,dim
        k = DOFs*(i-1)+j
        NrmSol % Values(k) = nrm*NrmSol % Values(k)
      END DO
    END IF
  END DO
!------------------------------------------------------------------------------     
 
  DEALLOCATE( ForceVector )
  Solver % Matrix % RHS => SaveRHS
  
  at2 = RealTime()
  WRITE(Message,* ) 'Solution Time: ',at2-at1
  CALL Info( 'NormalSolver', Message, Level=5 )
  
CONTAINS


!------------------------------------------------------------------------------
  SUBROUTINE BulkAssembly()
!------------------------------------------------------------------------------
       
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:,:), Basis(:)
    REAL(KIND=dp), POINTER :: ArrayPtr(:) => NULL()
    REAL(KIND=dp) :: Weight, Normal(3), detJ, Point(3), r(3)

    LOGICAL :: Found, CheckOrientation

    INTEGER :: elem,t,i,j,p,q,n,nd, Rank

    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    
    n = MAX( Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes )
    ALLOCATE( STIFF(n,n), FORCE(dim,n), Basis(n) )

    DO elem = 1,Solver % NumberOfActiveElements
      ! Element information
      ! ---------------------
      Element => GetActiveElement(elem)
      CALL GetElementNodes( Nodes )
      nd = GetElementNOFDOFs()
      n  = GetElementNOFNodes()

      ArrayPtr => ListGetConstRealArray1(GetBodyParams(Element), 'Point on Negative Side', CheckOrientation)
      IF (CheckOrientation) THEN
        Point = 0.0d0
        DO i=1,SIZE(ArrayPtr)
          Point(i) = ArrayPtr(i)
        END DO
      END IF

      ! Integrate local stresses:
      ! -------------------------
      IntegStuff = GaussPoints( Element )
      STIFF  = 0.0_dp
      FORCE  = 0.0_dp

      DO t=1,IntegStuff % n
        Found = ElementInfo( Element, Nodes, IntegStuff % u(t), &
           IntegStuff % v(t), IntegStuff % w(t), detJ, Basis )
        
        Weight = IntegStuff % s(t) * detJ
        IF ( CSymmetry ) Weight = Weight * &
              SUM( Basis(1:n) * Nodes % x(1:n) )
        
        IF ( .NOT. ConstantBulkMatrixInUse ) THEN
          DO p=1,nd
            DO q=1,nd
              STIFF(p,q) = STIFF(p,q) + Weight * Basis(q) * Basis(p)
            END DO
          END DO
        END IF
        
        Normal = NormalVector( Element, Nodes, &
           IntegStuff % u(t), IntegStuff % v(t), .TRUE. )
        IF (CheckOrientation) THEN
          r(1) = SUM(Basis(1:n) * Nodes % x(1:n)) - Point(1)
          r(2) = SUM(Basis(1:n) * Nodes % y(1:n)) - Point(2)
          r(3) = SUM(Basis(1:n) * Nodes % z(1:n)) - point(3)
          IF (SUM(Normal*r) < 0.0d0) Normal = -Normal
        END IF
        DO i=1,dim
          FORCE(i,1:nd) = FORCE(i,1:nd) + Weight*Normal(i)*Basis(1:nd)
        END DO
      END DO
      
!------------------------------------------------------------------------------
!      Update global matrices from local matrices 
!------------------------------------------------------------------------------
      IF ( .NOT. ConstantBulkMatrixInUse ) THEN
        Solver % Matrix % RHS => SaveRHS
        CALL DefaultUpdateEquations( STIFF, FORCE(1,1:nd) )
      END IF

      DO i=1,dim
        Solver % Matrix % RHS => ForceVector(:,i)
        CALL DefaultUpdateForce( FORCE(i,1:nd) )
      END DO
    END DO
    
    DEALLOCATE( STIFF, FORCE, Basis )
!------------------------------------------------------------------------------
  END SUBROUTINE BulkAssembly
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE NormalSolver
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Initialization for the primary solver: NormalSolver.
!> \ingroup Solvers
!------------------------------------------------------------------------------
  SUBROUTINE NormalSolver_Init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
    USE DefUtils

    TYPE(Model_t)  :: Model
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: DT
    LOGICAL :: Transient
!------------------------------------------------------------------------------
    INTEGER :: dim
    TYPE(ValueList_t), POINTER :: SolverParams
    CHARACTER(LEN=MAX_NAME_LEN) :: VarName
    LOGICAL :: Found
!------------------------------------------------------------------------------
    SolverParams => GetSolverParams()
    dim = CoordinateSystemDimension()

    IF ( .NOT. ListCheckPresent( SolverParams,'Variable') ) THEN
      CALL ListAddString( SolverParams, 'Variable','-nooutput nrm_temp' )
    END IF
    
    VarName = GetString(SolverParams,'Normals Result Variable', Found )
    IF( .NOT. Found ) THEN
      CALL ListAddString( SolverParams,'Normals Result Variable','Normals')
      CALL ListAddString( SolverParams,  'Exported Variable 1', 'Normals[Normals:2]' )
    END IF

    CALL ListAddInteger( SolverParams, 'Time derivative order', 0 )
!------------------------------------------------------------------------------
  END SUBROUTINE NormalSolver_Init
!------------------------------------------------------------------------------
