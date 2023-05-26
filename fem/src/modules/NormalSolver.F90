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
!> and this solution should give a smoother normal direction. Solver may also be used
!> to compute the normal component of displacement or velocity field. Then this
!> will be the result of the computation. 
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
  REAL(KIND=dp) :: at0,at1,at2
  TYPE(Variable_t), POINTER :: NrmSol
  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp), ALLOCATABLE :: Values(:)
  CHARACTER(*), PARAMETER :: Caller = 'NormalSolver'

  
  SAVE Visited

 
  CALL Info( Caller, '-------------------------------------',Level=4 )
  CALL Info( Caller, 'Computing the normals',Level=4 )
  CALL Info( Caller, '-------------------------------------',Level=4 )

  dim = CoordinateSystemDimension()
!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
  IF ( COUNT( Solver % Variable % Perm > 0 ) <= 0 ) RETURN


  SolverParams => GetSolverParams()
  Mesh => Solver % Mesh
  
  VarName = GetString(SolverParams,'Normals Result Variable',GotIt )
  IF(.NOT. GotIt) THEN
    CALL Fatal(Caller,'> Normals Result Variable < not found!')
  END IF

  NrmSol => VariableGet( Mesh % Variables,  VarName )
  IF(ASSOCIATED(NrmSol)) THEN
     Dofs = NrmSol % DOFs
     IF(Dofs /= DIM) THEN
       CALL Fatal(Caller,'The normals should have DOFs equal to DIM')
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
  
  CALL DefaultInitialize(Solver, ConstantBulkMatrixInUse)

  ALLOCATE(ForceVector(SIZE(Solver % Matrix % RHS),dim))  
  ForceVector = 0.0_dp
  SaveRHS => Solver % Matrix % RHS

  CALL BulkAssembly()
  IF (ConstantBulkMatrix) THEN
    CALL DefaultFinishBulkAssembly(BulkUpdate = .NOT.ConstantBulkMatrixInUse, RHSUpdate = .FALSE.)
  ELSE
    CALL DefaultFinishBulkAssembly()
  END IF

  !No flux BCs for this solver
  CALL DefaultFinishAssembly()

  at1 = RealTime()
  WRITE(Message,* ) 'Assembly Time: ',at1-at0
  CALL Info( Caller, Message, Level=5 )
!        
!------------------------------------------------------------------------------     
  SetD = GetLogical(SolverParams, 'Set Dirichlet Conditions',GotIt)
  IF ( SetD ) THEN
    ALLOCATE(Values(SIZE(Solver % Matrix % Values)))
    Vname = Solver % Variable % Name
    Values = Solver % Matrix % Values
  END IF
  DO i=1,dim
    Solver % Matrix % RHS => ForceVector(:,i)

    IF ( SetD ) THEN
      ! We compute the normals component-wise hence the Dirichlet condition also must be given so.
      Solver % Variable % Name = TRIM(VarName) // ' ' // I2S(i)
      Solver % Matrix % Values = Values
      CALL DefaultDirichletBCs()
      Solver % Variable % Name = Vname
    ELSE
      CALL DefaultDirichletBCs()
    END IF

    Solver % Variable % Values = 0._dp
    UNorm = DefaultSolve()

    NrmSol % Values(i::DOFs) = Solver % Variable % Values
  END DO
  Solver % Matrix % RHS => SaveRHS

  IF ( SetD ) DEALLOCATE(Values)

  ! Normalize the compute normal so that the length is always one.
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

  ! Optionally use the solver to project a vector field to normal direction
  BLOCK
    TYPE(Variable_t), POINTER :: Var1,Var2,Var3
    REAL(KIND=dp) :: nvec(3),dvec(3),dvec0(3)
    REAL(KIND=dp) :: maxd, limd, coeff
    INTEGER :: id,in      
    LOGICAL :: ScaleDt
    
    Vname = ListGetString( SolverParams,'Vector Field to Project',GotIt )
    IF(GotIt ) THEN 
      ScaleDt = ListGetLogical(SolverParams,'Vector Field Timestep scaling',GotIt )
    ELSE
      IF( ListGetLogical( SolverParams,'Project Displacement',GotIt ) ) THEN
        Vname = 'Displacement'
        ScaleDt = .FALSE.
      ELSE IF( ListGetLogical( SolverParams,'Project Velocity',GotIt ) ) THEN
        Vname = 'Velocity'
        ScaleDt = .TRUE.
      END IF
    END IF
    
    IF( GotIt ) THEN       
      CALL Info(Caller,'Projecting field "'//TRIM(Vname)//'" to normal direction',Level=5)

      ! We deal with displacement or velocity field with component-wise since
      ! otherwise we should deal with the 'flow solution' vector separately.
      Var1 => VariableGet( Mesh % Variables,TRIM(Vname)//' 1')
      IF(.NOT. ASSOCIATED(Var1) ) CALL Fatal(Caller,'Vector field component 1 does not exist: '//TRIM(Vname))
      Var2 => VariableGet( Mesh % Variables,TRIM(Vname)//' 2')
      IF(.NOT. ASSOCIATED(Var2) ) CALL Fatal(Caller,'Vector field component 2 does not exist: '//TRIM(Vname))
      IF(dim==3) THEN
        Var3 => VariableGet( Mesh % Variables,TRIM(Vname)//' 3')
        IF(.NOT. ASSOCIATED(Var3) ) CALL Fatal(Caller,'Vector field component 3 does not exist: '//TRIM(Vname))
      END IF

      dvec0(1) = ListGetCReal( SolverParams,'Vector Field Offset 1',GotIt)
      dvec0(2) = ListGetCReal( SolverParams,'Vector Field Offset 2',GotIt)
      dvec0(3) = ListGetCReal( SolverParams,'Vector Field Offset 3',GotIt)      
      
      nvec = 0.0_dp; dvec = 0.0_dp
      maxd = 0.0_dp
      DO i=1, Mesh % NumberOfNodes
        in = NrmSol % Perm(i)
        IF(in==0) CYCLE
        id = Var1 % Perm(i)
        IF(id==0) CYCLE

        DO j=1,dim
          k = DOFs*(in-1)+j
          nvec(j) = NrmSol % Values(k)
        END DO

        dvec(1) = Var1 % Values(id)
        dvec(2) = Var2 % Values(id)
        IF(dim==3) dvec(3) = Var3 % Values(id)

        dvec = dvec - dvec0
        dvec = SUM(dvec*nvec)*nvec

        maxd = MAX(maxd,SQRT(SUM(dvec(1:dim)**2)))
        
        DO j=1,dim
          k = DOFs*(in-1)+j
          NrmSol % Values(k) = dvec(j)
        END DO
      END DO

      ! If we have velocity then the normal displacement is given by
      ! d = dt*(n.v)n
      IF( ScaleDt .AND. Transient ) THEN
        NrmSol % Values = dt *  NrmSol % Values
        maxd = maxd * dt
      END IF

      limd = ListGetCReal( SolverParams,'Maximum Displacement',GotIt)
      IF( GotIt ) THEN
        IF( maxd > limd ) THEN
          coeff = limd / maxd
          WRITE(Message,'(A,ES12.3)') 'Limiting displacements with factor: ',coeff
          CALL Info(Caller,Message)
          NrmSol % Values = coeff *  NrmSol % Values         
        END IF
      END IF            
    END IF
      
  END BLOCK


  IF( ListGetLogical(SolverParams,'Enforce Symmetry',GotIt) ) THEN 
    CALL EnforceSymmetry()    
  END IF

  
  at2 = RealTime()
  WRITE(Message,* ) 'Solution Time: ',at2-at1
  CALL Info( Caller, Message, Level=5 )
  
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
    
    n = MAX( Mesh % MaxElementDOFs, Mesh % MaxElementNodes )
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
  SUBROUTINE EnforceSymmetry()
!------------------------------------------------------------------------------       
    LOGICAL :: Found    
    INTEGER :: elem,t,i,j,n,nd,m
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t), SAVE :: Nodes    

    CALL Info(Caller,'Enforcing symmetry of normal field',Level=6)

    m = 0
    DO elem = 1,Solver % NumberOfActiveElements
      ! Element information
      ! ---------------------
      Element => GetActiveElement(elem)
      CALL GetElementNodes( Nodes )
      nd = GetElementNOFDOFs()
      n  = GetElementNOFNodes()

      Found = .FALSE.
      DO i=1,n
        IF( ABS(Nodes % x(i)) < 1.0d-8) THEN
          Found = .TRUE.
          EXIT
        END IF
      END DO
      
      IF( Found ) THEN
        j = 3-i

        !PRINT *,'Indeces:',i,j,Element % NodeIndexes
        j = Element % NodeIndexes(j)
        i = Element % NodeIndexes(i)
        
        j = NrmSol % Perm(j)
        i = NrmSol % Perm(i)

        DO k=1,dim
          NrmSol % Values(Dofs*(i-1)+k) = NrmSOl % Values(Dofs*(j-1)+k)
        END DO
        m = m + 1
      END IF
    END DO

    CALL Info(Caller,'Enforced axial symmetry in nodes: '//I2S(m),Level=5)

!------------------------------------------------------------------------------
  END SUBROUTINE EnforceSymmetry
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

    CALL ListAddNewString( SolverParams, 'Variable','-nooutput nrm_temp' )
    
    VarName = GetString(SolverParams,'Normals Result Variable', Found )
    IF( .NOT. Found ) THEN
      CALL ListAddString( SolverParams,'Normals Result Variable','Normals')
      CALL ListAddString( SolverParams,  'Exported Variable 1', 'Normals[Normals:'//I2S(dim)//']' )
    END IF

    CALL ListAddInteger( SolverParams, 'Time derivative order', 0 )
!------------------------------------------------------------------------------
  END SUBROUTINE NormalSolver_Init
!------------------------------------------------------------------------------
