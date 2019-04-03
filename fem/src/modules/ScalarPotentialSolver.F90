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
!> Initialization for the primary solver: ScalarPotentialSolver
!> \ingroup Solvers
!------------------------------------------------------------------------------
  SUBROUTINE ScalarPotentialSolver_Init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
    USE DefUtils

    TYPE(Model_t)  :: Model
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: DT
    LOGICAL :: Transient
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: SolverParams
    INTEGER :: dim
!------------------------------------------------------------------------------
    SolverParams => GetSolverParams()
    dim = CoordinateSystemDimension()

    IF ( .NOT. ListCheckPresent( SolverParams,'Variable') ) THEN
      CALL ListAddInteger( SolverParams, 'Variable DOFs', 1 )
      CALL ListAddString( SolverParams,'Variable','Scalar Potential' )
    END IF
    CALL ListAddInteger( SolverParams, 'Time derivative order', 0 )

    ! Add linear system defaults: cg+ILU0
    IF(.NOT. ListCheckPresent(SolverParams,'Linear System Solver')) THEN
      CALL Info('ScalarPotentialSolver_init','Setting defaults for linear system solver')
      CALL ListAddString(SolverParams,'Linear System Solver','Iterative')
      IF(.NOT. ListCheckPresent(SolverParams,'Linear System Iterative Method')) &
          CALL ListAddString(SolverParams,'Linear System Iterative Method','cg')
      IF(.NOT. ListCheckPresent(SolverParams,'Linear System Preconditioning')) &
          CALL ListAddString(SolverParams,'Linear System Preconditioning','ILU0')
      IF(.NOT. ListCheckPresent(SolverParams,'Linear System Max Iterations')) &
          CALL ListAddInteger(SolverParams,'Linear System Max Iterations',500)
      IF(.NOT. ListCheckPresent(SolverParams,'Linear System Residual Output')) &
          CALL ListAddInteger(SolverParams,'Linear System Residual Output',10)
      IF(.NOT. ListCheckPresent(SolverParams,'Linear System Convergence Tolerance')) &
          CALL ListAddConstReal(SolverParams,'Linear System Convergence Tolerance',1.0e-10_dp)
    END IF
      
!------------------------------------------------------------------------------
  END SUBROUTINE ScalarPotentialSolver_Init
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Subroutine for computing the potential corresponding to a given flux.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE ScalarPotentialSolver( Model,Solver,dt,Transient )
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
  INTEGER :: i,j,dim,DOFs,ActiveDir
  LOGICAL :: ConstantBulkMatrix, ConstantBulkMatrixInUse, CSymmetry
  LOGICAL :: GotIt, Visited = .FALSE., DirMask(3)
  REAL(KIND=dp) :: Unorm, Totnorm
  REAL(KIND=dp), POINTER :: ForceVector(:,:), SaveRHS(:)
#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: at0,at1,at2
#else
  REAL(KIND=dp) :: at0,at1,at2,CPUTime,RealTime
#endif
  TYPE(Variable_t), POINTER :: ScalarPotentialSol
  
  SAVE Visited

 
  CALL Info( 'ScalarPotentialSolver', '-------------------------------------',Level=4 )
  CALL Info( 'ScalarPotentialSolver','Computing scalar potential of a vector field',Level=4 )
  CALL Info( 'ScalarPotentialSolver', '-------------------------------------',Level=4 )

  dim = CoordinateSystemDimension()
!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
  IF ( COUNT( Solver % Variable % Perm > 0 ) <= 0 ) RETURN
  
  SolverParams => GetSolverParams()
  
  CSymmetry = CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric
  
  VarName = GetString(SolverParams,'Flux Variable')
  CondName = ListGetString(SolverParams,'Flux Coefficient',GotIt)
  IF(.NOT. GotIt) CondName = 'none'

  at0 = RealTime()
  
  ConstantBulkMatrix = GetLogical( SolverParams, 'Constant Bulk Matrix', GotIt )
  ConstantBulkMatrixInUse = ConstantBulkMatrix .AND. &
      ASSOCIATED(Solver % Matrix % BulkValues)
  
  IF ( ConstantBulkMatrixInUse ) THEN
    Solver % Matrix % Values = Solver % Matrix % BulkValues        
    Solver % Matrix % rhs = 0.0_dp
  ELSE
    CALL DefaultInitialize()
  END IF
  
  CALL BulkAssembly()
  IF(.NOT. ConstantBulkMatrixInUse ) THEN
    CALL DefaultFinishBulkAssembly()
  END IF

  ! No flux BCs
  CALL DefaultFinishAssembly()

  CALL DefaultDirichletBCs()

  
  at1 = RealTime()
  WRITE(Message,* ) 'Assembly Time: ',at1-at0
  CALL Info( 'ScalarPotentialSolver', Message, Level=5 )
!        
!------------------------------------------------------------------------------     

  TotNorm = DefaultSolve()

!------------------------------------------------------------------------------     

  at2 = RealTime()
  WRITE(Message,* ) 'Solution Time: ',at2-at1
  CALL Info( 'ScalarPotentialSolver', Message, Level=5 )
  
  WRITE( Message, * ) 'Result Norm: ',TotNorm
  CALL Info( 'ScalarPotentialSolver', Message, Level=4 )
  
CONTAINS


!------------------------------------------------------------------------------
  SUBROUTINE BulkAssembly()
!------------------------------------------------------------------------------
       
    INTEGER :: elem,t,i,j,p,q,n,nd, Rank
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:)
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: weight,grad(3),C(3,3),coeff,detJ
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
    REAL(KIND=dp), ALLOCATABLE :: VectorField(:,:)
    LOGICAL :: Found, GotCoeff
    TYPE(ValueList_t), POINTER :: Material
    REAL(KIND=dp), POINTER :: Conductivity(:,:,:)=>NULL()
    
    SAVE Conductivity, Nodes
    
    n = MAX( Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes )
    ALLOCATE( STIFF(n,n), FORCE(n) )
    ALLOCATE( VectorField(3,n), Basis(n), dBasisdx(n,3) )

    DO elem = 1,Solver % NumberOFActiveElements
         
      ! Element information
      ! ---------------------
      Element => GetActiveElement(elem)
      CALL GetElementNodes( Nodes )
      nd = GetElementNOFDOFs()
      n  = GetElementNOFNodes()
      
      CALL GetRealArray( GetMaterial(), Conductivity, CondName, GotCoeff )
      IF ( GotCoeff ) THEN
        C = 0.0_dp
        Rank = GetTensorRank(Conductivity)
      END IF

      DO i=1,dim
        CALL GetScalarLocalSolution( VectorField(i,:), ComponentName(VarName,i) )
      END DO

      ! Integrate local stresses:
      ! -------------------------
      IntegStuff = GaussPoints( Element )
      STIFF  = 0.0_dp
      FORCE  = 0.0_dp
      
      
      DO t=1,IntegStuff % n
        Found = ElementInfo( Element, Nodes, IntegStuff % u(t), &
            IntegStuff % v(t), IntegStuff % w(t), detJ, Basis, dBasisdx )
        
        Weight = IntegStuff % s(t) * detJ
        IF ( CSymmetry ) Weight = Weight * SUM( Basis(1:n) * Nodes % x(1:n) )

       IF ( .NOT. ConstantBulkMatrixInUse ) THEN 
         IF(GotCoeff) THEN
           SELECT CASE(Rank)
           CASE(1)
             DO i=1,dim
               C(i,i) = SUM( Basis(1:n) * Conductivity(1,1,1:n) )
             END DO
           CASE(2)
             DO i=1,dim
               C(i,i) = SUM( Basis(1:n) * Conductivity(i,1,1:n) )
             END DO
           CASE DEFAULT
             DO i=1,dim
               DO j=1,dim
                 C(i,j) = SUM( Basis(1:n) * Conductivity(i,j,1:n) )
               END DO
             END DO
           END SELECT

           ! Use negative sign! 
           DO p=1,nd
             DO q=1,nd
               DO i=1,dim
                 DO j=1,dim
                   STIFF(p,q) = STIFF(p,q) - Weight * &
                       dBasisdx(q,i) * C(i,j) * dBasisdx(p,j)
                 END DO
               END DO
             END DO
           END DO
         ELSE
           DO p=1,nd           
             DO q=1,nd
               STIFF(p,q) = STIFF(p,q) + Weight * &
                   SUM( dBasisdx(q,1:dim) * dBasisdx(p,1:dim) )
             END DO
           END DO
         END IF
       END IF
        
       DO p=1,nd
         DO i=1,dim
           FORCE(p) = FORCE(p) + Weight * &
               dBasisdx(p,i) * SUM( Basis(1:nd) * VectorField(i,1:nd))
         END DO
       END DO
       
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
   
   DEALLOCATE( VectorField, STIFF, FORCE, Basis, dBasisdx )
!------------------------------------------------------------------------------
  END SUBROUTINE BulkAssembly
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  FUNCTION GetTensorRank( Tensor ) RESULT ( Rank )
!------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: Tensor(:,:,:)
    INTEGER :: Rank
    
    IF ( SIZE(Tensor,1) == 1 ) THEN
      Rank = 1
    ELSE IF ( SIZE(Tensor,2) == 1 ) THEN
      Rank = 2
    ELSE
      Rank = 3
    END IF
!-----------------------------------------------------------------------------
  END FUNCTION GetTensorRank
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
END SUBROUTINE ScalarPotentialSolver
!------------------------------------------------------------------------------

