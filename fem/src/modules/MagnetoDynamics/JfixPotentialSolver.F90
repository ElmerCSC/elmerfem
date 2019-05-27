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
! *  Authors: Juha Ruokolainen, Peter Råback
! *  Email:   elmeradm@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 08 Jun 1997
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> The current density given as a source must be divergence free to allow a 
!> hope for a solution for the magnetic vector potential. This solver may be used
!> to enforce a given current density to be divergence free by solving for an
!> equivalent potential a.k.a. Jfix such that the required fixing field is a
!> gradient of the potential.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE JfixPotentialSolver( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE MagnetoDynamicsUtils
  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL ::  Transient
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(ValueList_t),POINTER :: SolverParams, BF
  INTEGER :: i,j,k,n,m,dim,dofs
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Matrix_t), POINTER :: A => NULL(),B
  REAL(KIND=dp) :: Norm
  LOGICAL:: SingleNodeBC, EnsureBC, Found
  INTEGER, ALLOCATABLE :: Def_Dofs(:,:,:)
  CHARACTER(LEN=MAX_NAME_LEN):: Equation
  INTEGER, POINTER :: Perm(:)
  REAL(KIND=dp), POINTER :: fixpot(:),fixpotim(:),tmpsol(:),pJfixSurfaceVec(:)
  TYPE(Variable_t), POINTER :: jfixpot, jfixpotim, svar, IterV 
  LOGICAL :: ComplexSystem, Visited = .FALSE.

  SAVE :: A, Perm, Def_Dofs, jfixPot, jfixPotim
  
  CALL Info('JfixPotentialSolver','Computing fixing potential for given current density',Level=6)
  
  dim = CoordinateSystemDimension()
  Mesh => GetMesh()
  SolverParams => GetSolverParams()

  ! Take pointers to the master solver stuff
  B => GetMatrix()
  svar => Solver % Variable
  
  dofs = Solver % Variable % DOFs  
  ComplexSystem = ( dofs == 2 ) 
  
  jfixpot => VariableGet( Mesh % Variables, 'Jfix')

  IF( ComplexSystem ) THEN
    jfixpotim => VariableGet( Mesh % Variables, 'Jfix Im')
  END IF

  IF( .NOT. ASSOCIATED(jfixPot)) THEN    
    ALLOCATE(Perm(SIZE(Solver % Variable % Perm)))
    Perm = 0    
    Equation=GetString(SolverParams,'Equation',Found)
    
    ! Add default strategies for Jfix solver
    ! The AV equation typically prefers different ones. 
    CALL ListAddNewString(SolverParams,'Jfix: Linear System Solver', 'Iterative')
    CALL ListAddNewString(SolverParams,'Jfix: Linear System Iterative Method', 'BiCGStab')
    CALL ListAddNewLogical(SolverParams,'Jfix: Linear System Use HYPRE', .FALSE.)
    CALL ListAddNewLogical(SolverParams,'Jfix: Use Global Mass Matrix',.FALSE.)
    CALL ListAddNewString(SolverParams,'Jfix: Linear System Preconditioning', 'Ilu')
    CALL ListAddNewConstReal(SolverParams,'Jfix: Linear System Convergence Tolerance', &
        0.001_dp*GetCReal(SolverParams,'Linear System Convergence Tolerance', Found))
    CALL ListAddNewLogical(SolverParams,'Jfix: Skip Compute Nonlinear Change',.TRUE.)
    CALL ListAddNewLogical(SolverParams,'Jfix: Nonlinear System Consistent Norm',.TRUE.)
    CALL ListAddNewInteger(SolverParams,'Jfix: Linear System Residual Output',20)
    CALL ListAddNewString(SolverParams,'Jfix: Nonlinear System Convergence Measure','Norm')
    CALL ListAddNewLogical(SolverParams,'Jfix: Linear System Complex',.FALSE.)
    CALL ListAddNewLogical(SolverParams,'Jfix: Apply Conforming BCs',.FALSE.)
        
    n=SIZE(Solver % Def_Dofs,1)
    m=SIZE(Solver % Def_Dofs,2)
    k=SIZE(Solver % Def_Dofs,3)
    ALLOCATE(Def_Dofs(n,m,k))
    Def_Dofs = Solver % Def_Dofs

    Solver % Def_Dofs = 0
    Solver % Def_Dofs(:,:,1)=1

    ! Set namespace for matrix creation (to omit conforming BCs)
    CALL ListSetNameSpace('jfix:')    
    A => CreateMatrix( CurrentModel, Solver, Solver % Mesh, &
        Perm, 1, MATRIX_CRS, .TRUE., Equation, .FALSE., .FALSE.,&
        NodalDofsOnly = .TRUE.)          
    CALL ListSetNameSpace('')    
    n = A % NumberOfRows
    ALLOCATE(A % RHS(n))
    A % rhs = 0.0_dp
    A % Complex = .FALSE.
    
    ! Put pointers in the module so that these can be used externally also
    jfixRhs => A % Rhs 

    ! Create the variable for the (real part) fixing potential
    ALLOCATE(fixpot(n))
    fixpot = 0._dp
    CALL VariableAddVector( Mesh % Variables, Mesh, &
        Solver,'Jfix',1,fixpot,Perm)
    jfixpot => VariableGet(Mesh % Variables, 'Jfix')
    
    ! For complex cases create separately the variable for the imaginary component
    ! These are allocated component-wise so that we may easily solve them separately.
    ! Its the same equation for both with different load vector. 
    IF( ComplexSystem ) THEN
      ALLOCATE( jfixRhsC(n), fixpotim(n) )
      fixpotim = 0.0_dp
      CALL VariableAddVector( Mesh % Variables, Mesh, &
          Solver,'Jfix Im',1,fixpotim,Perm)
      jfixpotim => VariableGet(Mesh % Variables,'Jfix Im') 
    END IF

    ! Set potential only on a single node
    ! This uses the library functionality of DefaultDirichlet
    SingleNodeBC = GetLogical( SolverParams,'Single Node Projection BC',Found ) 
    IF( SingleNodeBC ) THEN
      DO i=1,Model % NumberOfBodyForces
        BF => Model % BodyForces(i) % Values
        IF( ListCheckPrefix( BF,'Current Density') ) THEN
          CALL ListAddConstReal( BF,'Jfix Single Node',0.0_dp )
        END IF
      END DO
    END IF
    
    CALL Info('JfixPotentialSolver','Finished creating matrix equation',Level=10)
  ELSE
    Solver % Def_Dofs = 0
    Solver % Def_Dofs(:,:,1)=1
  END IF
      
  Visited = .TRUE.
  
  Solver % Variable => jfixpot
  Solver % Matrix => A
  IF(ParEnv % PEs > 1) CALL ParallelInitMatrix(Solver,A)
  
  IF( JfixPhase == 1 ) THEN
    A % Values = 0.0_dp
    A % rhs = 0.0_dp
    IF( ComplexSystem ) THEN
      JfixRhsC = CMPLX( 0.0_dp, 0.0_dp )
    END IF

    !CALL DefaultInitialize()
    CALL JfixBulkAssembly()
  
    IF(.NOT. ASSOCIATED( JfixSurfacePerm ) ) THEN
      EnsureBC = .NOT. ListGetLogical( SolverParams,'Jfix without boundaries',Found)
      CALL MarkOuterNodes(Mesh,Perm,n,JfixSurfacePerm,EnsureBC)
      IF( ComplexSystem ) THEN
        ALLOCATE( JfixSurfaceVecC(3*n) )    
      ELSE
        ALLOCATE( JfixSurfaceVec(3*n) )
        pJfixSurfaceVec => JfixSurfaceVec
        IF( ListGetLogical( SolverParams,'Jfix Surface Source Save') ) THEN      
          CALL VariableAddVector( Mesh % Variables,Mesh,Solver,'Jfix Surface Source',&
              3,pJfixSurfaceVec,JfixSurfacePerm)
        END IF      
        IF( ListGetLogical( SolverParams,'Jfix rhs Save') ) THEN      
          CALL VariableAddVector( Mesh % Variables,Mesh,Solver,'Jfix rhs',&
              1,A % rhs,Perm)
        END IF      
      END IF
    END IF
    IF( ComplexSystem ) THEN
      JfixSurfaceVecC = CMPLX( 0.0_dp, 0.0_dp ) 
    ELSE
      JfixSurfaceVec = 0.0_dp
    END IF

    
  ELSE IF( JfixPhase == 2 ) THEN
    CALL ListSetNameSpace('jfix:')    

    ! Set load before applying Dircihlet condistions
    IF( ComplexSystem ) A % rhs = REAL( JfixRhsC )      
        
    CALL JfixBCs()

    ! We can use the same BCs for real and complex currents i.e. always "Jfix"
    CALL DefaultDirichletBCs()

    n = 0
    IF( ALLOCATED( A % ConstrainedDOF ) ) THEN
      n = COUNT( A % ConstrainedDOF )
      n = NINT( ParallelReduction( 1.0_dp * n ) )
    END IF
    IF( n == 0 ) THEN
      CALL Warn('JfixPotentialSolver','No Dirichlet conditions used to define Jfix level!')
    ELSE
      CALL Info('JfixPotentialSolver','Number of dirichlet nodes: '//TRIM(I2S(n)),Level=7)
    END IF
    
    CALL Info('JfixPotentialSolver','Solving for Jfix',Level=10)
    ! Set load before applying Dircihlet condistions
    IF( ComplexSystem ) THEN
      A % rhs = REAL( JfixRhsC )      
      IF( ALLOCATED( A % ConstrainedDOF ) ) THEN
        WHERE( A % ConstrainedDOF ) A % rhs = 0.0_dp
      END IF
    END IF
    CALL SolveSystem(A,ParMatrix,A % rhs,jfixpot % Values,jfixpot % Norm,1,Solver)

    WRITE(Message,'(A,ES12.3)') 'Norm for Jfix computation: ',SUM( ABS( jfixpot % Values ) )
    CALL Info('JfixPotentialSolver',Message,Level=8)

    IF( ComplexSystem ) THEN
      CALL Info('JfixPotentialSolver','Solving for imaginary component of Jfix',Level=10)      
      A % rhs = AIMAG( JfixRhsC )
      IF( ALLOCATED( A % ConstrainedDOF ) ) THEN
        WHERE( A % ConstrainedDOF ) A % rhs = 0.0_dp
      END IF
      CALL SolveSystem(A,ParMatrix,A % rhs,jfixpotim % values,jfixpotim % Norm,1,Solver)
      
      ! This is temporal norm for debugging
      WRITE(Message,'(A,ES12.3)') 'Norm for Jfix Im computation: ',SUM( ABS( jfixpotim % Values ) )
      CALL Info('JfixPotentialSolver',Message,Level=8)
    END IF   
    
    CALL ListSetNameSpace('')

    IterV => VariableGet( Solver % Mesh % Variables, 'nonlin iter' )
    IterV % Values(1) = 1
  ELSE
    CALL Fatal('JfixPotentialSolver','Invalid JfixPhase')
  END IF
  
  Solver % Variable => svar
  Solver % Matrix => B
  Solver % Def_Dofs = Def_Dofs    

    
CONTAINS

  
!------------------------------------------------------------------------------
! Assemble nodal Poisson equation (matrix part only) for the solution of the
! Jfix field. The complex system is solved component-wise so this is always real.
!------------------------------------------------------------------------------
  SUBROUTINE JfixBulkAssembly()
!------------------------------------------------------------------------------
    IMPLICIT NONE       
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:)
    INTEGER :: elem,t,p,q,n
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    TYPE(Nodes_t) :: Nodes
    LOGICAL :: Found
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: weight,detJ
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)

    SAVE Nodes

    n = MAX(Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes)
    ALLOCATE( STIFF(n,n), FORCE(n),Basis(n), dBasisdx(n,3) )
    FORCE = 0.0_dp    
    
    DO elem = 1,GetNOFActive()
      ! Element information
      ! ---------------------
      Element => GetActiveElement(elem)
      
      IF( Element % PartIndex /= ParEnv % MyPe ) CYCLE      
      
      CALL GetElementNodes( Nodes )
      n  = GetElementNOFNodes()

      STIFF = 0._dp

      IntegStuff = GaussPoints( Element )
      
      DO t=1,IntegStuff % n
        Found = ElementInfo( Element, Nodes, IntegStuff % u(t), &
            IntegStuff % v(t), IntegStuff % w(t), detJ, Basis, dBasisdx )

        Weight = IntegStuff % s(t) * detJ
        DO p=1,n
          DO q=1,n
            STIFF(p,q) = STIFF(p,q) + Weight * SUM(dBasisdx(q,:)*dBasisdx(p,:))
          END DO
        END DO
      END DO
      
      CALL DefaultUpdateEquations(STIFF, FORCE)
    END DO
      
!------------------------------------------------------------------------------
  END SUBROUTINE JfixBulkAssembly
!------------------------------------------------------------------------------

  
!------------------------------------------------------------------------------
! This subroutine fixes the potential to zero where there is significant current density
! component in the normal direction. Alternatively, the user may set mixed conditions
! such that positive flux yields Dirichlet conditions and negative ones Neumann conditions. 
!------------------------------------------------------------------------------
  SUBROUTINE JfixBCs()
!------------------------------------------------------------------------------
    IMPLICIT NONE       
    INTEGER :: i,j,k1,k2,t,n,meshdim,ActParents,ParParents
    TYPE(Nodes_t) :: Nodes
    LOGICAL :: Found, JfixHybrid, JfixNeu, JfixDir, JfixAuto, JfixStatCurr
    TYPE(Element_t), POINTER :: Element, LElement, P1, P2
    REAL(KIND=dp) :: NrmEps,Jlen,JVec(3),Nrm(3),NrmProj,MaxProj,&
        MaxJVec,JEps,Jrel,Jabs,Jnrm
    TYPE(ValueList_t), POINTER :: BC
    
    SAVE Nodes

    meshdim = Solver % Mesh % MeshDim    
    
    JfixNeu = GetLogical( SolverParams,'Jfix Neumann BCs',Found ) 
    JfixDir = GetLogical( SolverParams,'Jfix Dirichlet BCs',Found )
    JfixHybrid = GetLogical( SolverParams,'Jfix Hybrid BCs',Found )

    JfixAuto = GetLogical( SolverParams,'Jfix Automated BCs',Found )
    IF(.NOT. Found ) JfixAuto = GetLogical( SolverParams, &
        'Automated Source Projection BCs', Found )

    JfixStatCurr = .FALSE.
    IF( Found ) THEN
      ! Old automated flag set false => use StatCurrentSolver mode
      JfixStatCurr = .NOT. JfixAuto
    ELSE IF(.NOT. (JfixNeu .OR. JfixDir .OR. JfixHybrid )) THEN
      JfixAuto = .TRUE.
      JfixStatCurr = .FALSE.
    END IF
      
      
    ! This keyword fixes all outer boundaries with same strategy.
    ! If we set Dirichlet conditions it is needless to set anything else.
    ! If we set Neumann conditions with the hubrid strategy we may still combine
    ! with selected Dirichlet nodes.
    IF( JfixDir .OR. JfixNeu .OR. JfixHybrid ) THEN
      DO j=1,Mesh % NumberOfNodes
        k1 = Perm(j)
        k2 = JfixSurfacePerm(j)
        IF( k1 == 0 .OR. k2 == 0 ) CYCLE
        IF( JfixDir ) THEN
          CALL UpdateDirichletDof(A, k1, 0._dp)
        ELSE
          IF( ComplexSystem ) THEN
            JfixRhsC(k1) = 0.0_dp
          ELSE
            JfixRhs(k1) = 0.0_dp
          END IF
        END IF
      END DO
      IF( JfixDir .OR. JfixNeu) RETURN
    END IF
    
    
    ! Find the tolerances for detecting automatically the Jfix BCs.
    IF( JfixAuto .OR. JfixHybrid ) THEN
      NrmEps = GetCReal(SolverParams, 'Jfix norm eps', Found)
      IF (.NOT. Found) NrmEps = 0.5_dp

      Jrel = GetCReal(SolverParams, 'Jfix relative eps', Found)
      IF (.NOT. Found) Jrel = 1.0e-6

      Jabs = GetCReal(SolverParams, 'Jfix absolute eps', Found)
      IF (.NOT. Found) Jabs = EPSILON( Jabs ) 

      MaxProj = 0.0_dp
      IF( ComplexSystem ) THEN
        MaxJVec = MAXVAL( ABS( JfixSurfaceVecC ) )
      ELSE      
        MaxJVec = MAXVAL( ABS( JfixSurfaceVec ) )
      END IF
      MaxJVec = ParallelReduction( MaxJVec, 2 ) 
      WRITE( Message,'(A,ES12.3)') 'Maximum source term on boundaries:',MaxJVec
      CALL Info('JfixBCs',Message,Level=8)

      JEps = MAX( Jabs, Jrel * MaxJVec )
      WRITE( Message,'(A,ES12.3)') 'Using jfix epsilon for flux:',Jeps
      CALL Info('JfixBCs',Message,Level=8)   
    END IF

      
    DO t=1,GetNOFBoundaryElements()
      Element => GetBoundaryElement(t)
      n = GetElementNOFNodes()

      IF (Element % TYPE % DIMENSION < meshdim-1 ) CYCLE

      IF( ANY( Perm(Element % NodeIndexes) == 0 ) ) CYCLE
      IF( ANY( JfixSurfacePerm(Element % NodeIndexes) == 0 ) ) CYCLE

      ! We may choose to set some BC to be Neumann overriding the automated strategy
      BC => GetBC( Element )
      IF( ASSOCIATED( BC ) ) THEN
        JfixNeu = ListGetLogical( BC,'Jfix Neumann BC',Found )

        IF(.NOT. Found .AND. JfixStatCurr ) THEN
          ! This is here to have backward compatibility with old non-automated strategy
          ! intended to be used with StatCurrentSolver.
          JfixNeu = ListCheckPrefix( BC,'Current Density')
        END IF
        IF( JfixNeu ) THEN
          IF( ComplexSystem ) THEN
            JfixRhsC(JfixSurfacePerm(Element % NodeIndexes)) = 0.0_dp
          ELSE
            JfixRhs(JfixSurfacePerm(Element % NodeIndexes)) = 0.0_dp
          END IF
          CYCLE
        END IF        
      END IF

      IF( .NOT. ( JfixAuto .OR. JfixHybrid ) ) CYCLE
      CALL GetElementNodes(Nodes)
      Nrm = NormalVector(Element,Nodes,Check=.TRUE.)
      
      DO i=1,n
        j = Element % NodeIndexes(i) 
        k1 = Perm(j)
        k2 = JfixSurfacePerm(j)

        IF( ComplexSystem ) THEN
          JVec = ABS( JfixSurfaceVecC(3*k2-2:3*k2) )
        ELSE
          JVec = JfixSurfaceVec(3*k2-2:3*k2)
        END IF

        JLen = SQRT( SUM( JVec**2 ) )

        IF( Jlen < Jeps ) CYCLE          
        Jnrm = SUM( Nrm * Jvec ) 

        ! This is to test whether to only fix positive flux BCs and 
        IF( JfixHybrid .AND. Jnrm < 0.0 ) CYCLE

        NrmProj = ABS( Jnrm ) / Jlen
        MaxProj = MAX( MaxProj, NrmProj ) 

        IF( NrmProj > Nrmeps ) THEN
          CALL UpdateDirichletDof(A, k1, 0._dp)
        END IF
      END DO
    END DO

    IF( JfixAuto .OR. JfixHybrid ) THEN
      WRITE( Message,'(A,ES12.3)') 'Maximum norm projection:',MaxProj
      CALL Info('JfixBCs',Message,Level=15)           
    END IF
      
!------------------------------------------------------------------------------
  END SUBROUTINE JfixBCs
!------------------------------------------------------------------------------
  

!------------------------------------------------------------------------------
END SUBROUTINE JfixPotentialSolver
!------------------------------------------------------------------------------
