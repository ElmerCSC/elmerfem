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
!


!-----------------------------------------------------------------------------
!>  
!>  , and upper and lower limiters.
!> \ingroup Solvers
!-----------------------------------------------------------------------------
SUBROUTINE ThicknessSolver( Model,Solver,dt,TransientSimulation )
  USE DefUtils
  USE Differentials
  USE MaterialModels
  IMPLICIT NONE

  !------------------------------------------------------------------------------
  !    external variables
  !------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t):: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
  !    Local variables
  !------------------------------------------------------------------------------

  LOGICAL ::&
       firstTime=.TRUE., Found, AllocationsDone = .FALSE., stat, &
       LimitDisp,  Bubbles = .False.,&
       SubstantialSurface = .TRUE.,&
       UseBodyForce = .TRUE., ApplyDirichlet=.FALSE.,  ALEFormulation=.FALSE. , &
       ConvectionVar,Compute_dhdt,UnFoundFatal=.TRUE.
  LOGICAL :: Passive
  LOGICAL, ALLOCATABLE ::  LimitedSolution(:,:), ActiveNode(:,:)

  INTEGER :: & 
       i,j,K,L, p, q, R, t,N,NMAX,MMAX,nfamily, deg, Nmatrix,&
       edge, bf_id,DIM,istat,LocalNodes,nocorr,&
       NSDOFs,NonlinearIter,iter, numberofsurfacenodes
  INTEGER, POINTER ::&
       ThickPerm(:), DHDTPrem(:),FlowPerm(:), NodeIndexes(:), EdgeMap(:,:)

#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: &
       at,st,totat,totst,Norm,PrevNorm,LocalBottom, cv, &
       Relax, MaxDisp, maxdh,LinearTol,NonlinearTol,RelativeChange,&
       smallestpossiblenumber, rr, ss
#else
  REAL(KIND=dp) :: &
       at,st,totat,totst,CPUTime,Norm,PrevNorm,LocalBottom, cv, &
       Relax, MaxDisp, maxdh,LinearTol,NonlinearTol,RelativeChange,&
       smallestpossiblenumber, rr, ss
#endif

  REAL(KIND=dp), POINTER :: ForceVector(:), Thick(:),DHDT(:),PreH(:,:), &
       FlowSolution(:),  PointerToResidualVector(:)

  REAL(KIND=dp), ALLOCATABLE :: ResidualVector(:), &
       STIFF(:,:),FORCE(:), TimeForce(:), LOAD(:),&
       MASS(:,:), Velo(:,:),  LowerLimit(:), UpperLimit(:), &
       OldValues(:), OldRHS(:),StiffVector(:),MeshVelocity(:,:)

  CHARACTER(LEN=MAX_NAME_LEN)  :: SolverName, VariableName, EquationName, FlowSolName, StabilizeFlag

  TYPE(Nodes_t)   :: ElementNodes
  TYPE(Element_t),POINTER :: CurrentElement
  TYPE(Variable_t), POINTER :: FlowSol, VarThickResidual,DHDTSol
  TYPE(ValueList_t), POINTER :: BodyForce, SolverParams, Material, Equation
  TYPE(Matrix_t), POINTER :: Systemmatrix
  !-----------------------------------------------------------------------------
  !      remember these variables
  !----------------------------------------------------------------------------- 
  SAVE STIFF, MASS, FORCE, &
       LOAD, &
       ElementNodes, AllocationsDone, Velo,  TimeForce, &
       UseBodyForce, LimitedSolution, LowerLimit, UpperLimit, ActiveNode, OldValues, OldRHS, &
       ResidualVector, StiffVector, MeshVelocity

  !------------------------------------------------------------------------------
  !    Get variabel/solver name
  !------------------------------------------------------------------------------
  VariableName = TRIM(Solver % Variable % Name)
  SolverName = 'ThicknessSolver ('// TRIM(Solver % Variable % Name) // ')'

  !    Get variables for the solution
  !------------------------------------------------------------------------------
  Thick     => Solver % Variable % Values     ! Nodal values for free surface displacement
  IF (.NOT.ASSOCIATED(Thick)) CALL Fatal(SolverName,'Variable values not associated')
  ThickPerm => Solver % Variable % Perm       ! Permutations for free surface displacement
  PreH => Solver % Variable % PrevValues

  !------------------------------------------------------------------------------
  !    Get constants and solver params
  !------------------------------------------------------------------------------
  DIM = CoordinateSystemDimension()
  smallestpossiblenumber = TINY(smallestpossiblenumber)
  SolverParams => GetSolverParams()

  SystemMatrix => Solver % Matrix
  ForceVector => Solver % Matrix % RHS


  LinearTol = GetConstReal( SolverParams, &
       'Linear System Convergence Tolerance',    Found )
  IF ( .NOT.Found ) THEN
     CALL Fatal(SolverName, 'No >Linear System Convergence Tolerance< found')
  END IF
  NonlinearTol  = GetConstReal( SolverParams, &
       'Nonlinear System Convergence Tolerance',    Found )
  IF ( .NOT.Found ) NonlinearTol = 1.0e-6
  NonlinearIter = GetInteger(   SolverParams, &
       'Nonlinear System Max Iterations', Found )
  IF ( .NOT.Found ) NonlinearIter = 1

  Compute_dhdt = GetLogical( SolverParams, &
         'Compute DHDT', Found)
  IF ( .NOT.Found ) Compute_dhdt=.False.

  ApplyDirichlet = GetLogical( SolverParams, &
       'Apply Dirichlet', Found)
  IF ( .NOT.Found ) THEN
     ApplyDirichlet = .FALSE.
     CALL Info(SolverName, 'No keyword > Apply Dirichlet < found. No limitation of solution',Level=6 )
  ELSE
     IF (ApplyDirichlet) THEN
        CALL Info(SolverName, 'Using Dirichlet method for limitation',Level=6 )
        IF (NonlinearIter < 2) THEN
           CALL Warn(SolverName, 'Keyword > Apply Dirichlet < set, but > Nonlinear System Max Iterations < set to lower than 2')
        END IF
     ELSE
        CALL Info(SolverName, 'No limitation of solution',Level=6 )
     END IF
  END IF

  ALEFormulation = GetLogical( SolverParams, &
       'ALE Formulation', Found)
  IF ( .NOT.Found ) THEN
     ALEFormulation = .FALSE.
  END IF
  IF (ALEFormulation) THEN 
     CALL Info(SolverName, 'Using horizontal ALE Formulation',Level=6 )
  ELSE
     CALL Info(SolverName, 'Using horizontal Eulerian Formulation',Level=6 )
  END IF

  StabilizeFlag = GetString( SolverParams, &
       'Stabilization Method',Found )
  SELECT CASE(StabilizeFlag)
     CASE('stabilized')
        Bubbles = .FALSE.
     CASE('bubbles')
        Bubbles = .TRUE.
     CASE DEFAULT
        Bubbles = .FALSE.
  END SELECT
  IF (Bubbles) THEN
     CALL Info(SolverName, 'Using residual free bubble stabilization',Level=6 )
  ELSE
     CALL Info(SolverName, &
          'Using residual squared-stabilized formulation.',Level=6 )
  END IF


  WRITE(Message,'(A,I0)') 'Mesh dimension: ', DIM
  CALL Info( SolverName, Message, Level=8 )

  !------------------------------------------------------------------------------
  !    Allocate some permanent storage, this is done first time only
  !------------------------------------------------------------------------------

  IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed) THEN
     NMAX = Model % MaxElementNodes
     MMAX = Model % Mesh % NumberOfNodes 
     K = SIZE( SystemMatrix % Values )
     L = SIZE( SystemMatrix % RHS )

     IF ( AllocationsDone ) THEN
        DEALLOCATE( ElementNodes % x,    &
             ElementNodes % y,    &
             ElementNodes % z,    &
             TimeForce,        &
             FORCE,    &
             STIFF, &
             MASS,  &
             LOAD,&
             Velo,  &
             MeshVelocity, &
             LowerLimit,                      &
             UpperLimit, &
             LimitedSolution,  &
             ActiveNode,                      & 
             ResidualVector, &
             OldValues,&
             OldRHS,&
             StiffVector)

     END IF


     IF (Bubbles) THEN
        Nmatrix = 2*NMAX
     ELSE
        Nmatrix = NMAX
     END IF

     ALLOCATE( ElementNodes % x( NMAX ),    &
          ElementNodes % y( NMAX ),    &
          ElementNodes % z( NMAX ),    &
          TimeForce( Nmatrix ),        &
          FORCE( Nmatrix ),    &
          STIFF( Nmatrix, Nmatrix ), &
          MASS( Nmatrix, Nmatrix ),  &
          LOAD(NMAX) , &
          Velo( 3, NMAX ), &
          MeshVelocity( 3,NMAX ), &
          LowerLimit( MMAX ), &
          UpperLimit( MMAX ), &
          LimitedSolution( MMAX, 2 ),  &
          ActiveNode( MMAX, 2 ),       &  
          ResidualVector( L ),         &
          OldValues(K), &
          OldRHS(L), &
          StiffVector( L ), &
          STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal(SolverName,'Memory allocation error 1, Aborting.')
     END IF

     CALL Info(SolverName,'Memory allocations done' )
     AllocationsDone = .TRUE.
     ActiveNode = .FALSE.
     ResidualVector = 0.0_dp
  END IF
  LimitedSolution=.FALSE.


  !------------------------------------------------------------------------------
  !    Get variables for the residual
  !------------------------------------------------------------------------------
  VarThickResidual => VariableGet( Model % Mesh % Variables, TRIM(VariableName) // ' Residual',UnFoundFatal=UnFoundFatal)

  PointerToResidualVector => VarThickResidual % Values

  !------------------------------------------------------------------------------
  !    Get Flow solution
  !------------------------------------------------------------------------------
   ConvectionVar=.True.
   FlowSolName =  GetString(SolverParams ,'Flow Solution Name', Found)
   IF(.NOT.Found) THEN        
        WRITE(Message,'(A)') &
           '<Flow Solution Name> Not Found; will look for <convection velocity> in body forces'
        CALL Info(SolverName,Message,level=10)
        ConvectionVar=.False.
        NSDOFS=GetInteger(SolverParams ,'Convection Dimension',Found)
        IF(.NOT.Found) &
            CALL Fatal(SolverName,'if <Flow Solution Name> not given prescribe <Convection Dimension>')
   ELSE
        FlowSol => VariableGet( Solver % Mesh % Variables, FlowSolName,UnFoundFatal=UnFoundFatal)
        FlowPerm     => FlowSol % Perm
        NSDOFs     =  FlowSol % DOFs
        FlowSolution => FlowSol % Values
   END IF

  !------------------------------------------------------------------------------
  ! Non-linear iteration loop
  !------------------------------------------------------------------------------
  DO iter=1,NonlinearIter
     !------------------------------------------------------------------------------
     !    assign matrices
     !------------------------------------------------------------------------------
     LocalNodes = Model % NumberOfNodes
     !Norm = Solver % Variable % Norm     
     WRITE(Message,'(a,I0,a,I0)') 'Non-linear Iteration ', iter,' out of max. ',NonlinearIter
     CALL Info( SolverName, Message, Level=4)
     !------------------------------------------------------------------------------
     !    Do some additional initialization, and go for it
     !------------------------------------------------------------------------------
     totat = 0.0_dp
     totst = 0.0_dp
     at = CPUTime()
     CALL Info( SolverName, 'start assembly',Level=6 )
     CALL DefaultInitialize()
     !------------------------------------------------------------------------------
     !    Do the assembly
     !------------------------------------------------------------------------------
     DO t=1,Solver % NumberOfActiveElements
        CurrentElement => GetActiveElement(t)
        Passive=CheckPassiveElement(CurrentElement)
        n = GetElementNOFNodes()
        NodeIndexes => CurrentElement % NodeIndexes

        ! set coords of highest occurring dimension to zero (to get correct path element)
        !-------------------------------------------------------------------------------
        ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
        IF (NSDOFs == 1) THEN
           ElementNodes % y(1:n) = 0.0
           ElementNodes % z(1:n) = 0.0
        ELSE IF (NSDOFs == 2) THEN
           ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
           ElementNodes % z(1:n) = 0.0_dp
        ELSE
           WRITE(Message,'(a,i0,a)')&
                'It is not possible to compute Thickness evolution if Flow Sol DOFs=',&
                NSDOFs, ' . Aborting'
           CALL Fatal( SolverName, Message) 
           STOP   
        END IF


        ! get pointers on Equation, Material and body-Force section input
        !----------------------------------------------------------------
        Equation => GetEquation()
        Material => GetMaterial()
        BodyForce => GetBodyForce()


        ! get lower limit for solution 
        !-----------------------------
        LowerLimit(CurrentElement % Nodeindexes(1:N)) = &
             ListGetReal(Material,'Min ' // TRIM(VariableName),n,CurrentElement % NodeIndexes, Found) 
        IF (.NOT.Passive) LimitedSolution(CurrentElement % Nodeindexes(1:N), 1) = Found
        ! get upper limit for solution 
        !-----------------------------
        UpperLimit(CurrentElement % Nodeindexes(1:N)) = &
             ListGetReal(Material,'Max ' // TRIM(VariableName),n,CurrentElement % NodeIndexes, Found)              
        IF (.NOT.Passive) LimitedSolution(CurrentElement % Nodeindexes(1:N), 2) = Found

        ! get flow soulution and velocity field from it
        !----------------------------------------------
        Velo = 0.0_dp
        !----------------------------------------------------

        ! get velocity profile
        IF (ConvectionVar) Then
          DO i=1,n
             j = NSDOFs*FlowPerm(NodeIndexes(i))
              !2D problem - 1D Thickness evolution
              IF((DIM == 2) .AND. (NSDOFs == 1)) THEN 
                 Velo(1,i) = FlowSolution( j ) 
                 Velo(2,i) = 0.0_dp
              !2D problem - 2D Thickness evolution (plane view pb)
              ELSE IF ((DIM == 2) .AND. (NSDOFs == 2)) THEN
                 Velo(1,i) = FlowSolution( j-1 ) 
                 Velo(2,i) = FlowSolution( j ) 
              !3D problem - 2D Thickness evolution 
              ELSE IF ((DIM == 3) .AND. (NSDOFs == 2)) THEN
                 Velo(1,i) = FlowSolution( j-1 ) 
                 Velo(2,i) = FlowSolution( j ) 
              ELSE
                 WRITE(Message,'(a,i0,a,i0,a)')&
                      'DIM=', DIM, ' NSDOFs=', NSDOFs, ' does not combine. Aborting'
                 CALL Fatal( SolverName, Message)
              END IF
           END DO
       ELSE
          IF (ASSOCIATED( BodyForce ) ) THEN
               Velo(1,1:n) = GetReal( BodyForce, 'Convection Velocity 1',Found )
               if (NSDOFs.eq.2) Velo(2,1:n) = GetReal( BodyForce, 'Convection Velocity 2',Found )
          END IF
        END IF
        !------------------------------------------------------------------------------
        ! Get mesh velocity
        !------------------------------------------------------------------------------
        MeshVelocity = 0.0_dp
        CALL GetVectorLocalSolution( MeshVelocity, 'Mesh Velocity',CurrentElement)
        !

        !------------------------------------------------------------------------------
        !      get the accumulation/ablation rate (i.e. normal surface flux)
        !      from the body force section
        !------------------------------------------------------------------------------
        LOAD=0.0_dp
        IF (ASSOCIATED( BodyForce ) ) THEN
              LOAD(1:n) = LOAD(1:n) +   &
                      GetReal( BodyForce, 'Top Surface Accumulation', Found )
              LOAD(1:n) = LOAD(1:n) +   &
                      GetReal( BodyForce, 'Bottom Surface Accumulation', Found )
        END IF


        !------------------------------------------------------------------------------
        !      Get element local matrix, and rhs vector
        !------------------------------------------------------------------------------
        CALL LocalMatrix( STIFF, MASS, FORCE,&
             LOAD,  Velo, NSDOFs, MeshVelocity, &
             CurrentElement, n, ElementNodes, NodeIndexes, &
             TransientSimulation,&
              ALEFormulation)

        !------------------------------------------------------------------------------
        !      If time dependent simulation add mass matrix to stiff matrix
        !------------------------------------------------------------------------------
        TimeForce = 0.0_dp
        IF ( TransientSimulation ) THEN
           !------------------------------------------------------------------------------
           !        NOTE: This will replace STIFF and LocalForce with the
           !              combined information...
           !------------------------------------------------------------------------------
           CALL Default1stOrderTime( MASS, STIFF, FORCE )
        END IF
        !------------------------------------------------------------------------------
        !      Update global matrices from local matrices
        !------------------------------------------------------------------------------
        IF (Bubbles) CALL Condensate( N, STIFF, FORCE, TimeForce )
        !------------------------------------------------------------------------------
        !      Update global matrix and rhs vector from local matrix & vector
        !------------------------------------------------------------------------------
        CALL DefaultUpdateEquations( STIFF, FORCE )
        !------------------------------------------------------------------------------
     END DO ! End loop bulk elements
     CALL DefaultFinishBulkAssembly()

     !------------------------------------------------------------------------------
     !     Neumann & Newton boundary conditions
     !------------------------------------------------------------------------------
     !
     ! MIND: In weak formulation it is not possible to prescribe a contact angle on
     !       a boundary in this solver. This has to be taken care of in the boundary
     !       condition for the stress tensor in the Navier-Stokes Solver. Thus, in
     !       generally it does not make sense to prescribe a Neumann type of
     !       condition here.

     !------------------------------------------------------------------------------
     !    FinishAssemebly must be called after all other assembly steps, but before
     !    Dirichlet boundary settings. Actually no need to call it except for
     !    transient simulations.
     !------------------------------------------------------------------------------
     CALL DefaultFinishAssembly()
     CALL DefaultDirichletBCs()

     !------------------------------------------------------------------------------
     !    Manipulation of the assembled matrix due to limits
     !------------------------------------------------------------------------------
     OldValues = SystemMatrix % Values
     OldRHS = ForceVector

     IF (ApplyDirichlet) THEN
        ! manipulation of the matrix
        !---------------------------
        DO i=1,Model % Mesh % NumberOfNodes
           k = ThickPerm(i)           
           IF ((ActiveNode(i,1) .OR. ActiveNode(i,2)) .AND. (k > 0)) THEN
              CALL ZeroRow( SystemMatrix, k ) 
              CALL SetMatrixElement( SystemMatrix, k, k, 1.0_dp ) 
              IF(ActiveNode(i,1)) THEN
                 SystemMatrix % RHS(k) = LowerLimit(i)
              ELSE
                 SystemMatrix % RHS(k) = UpperLimit(i)
              END IF
           END IF
        END DO
     END IF

        CALL Info( SolverName, 'Assembly done', Level=6 )
        !------------------------------------------------------------------------------
        !    Solve System  and check for convergence
        !------------------------------------------------------------------------------
        at = CPUTime() - at
        st = CPUTime() 

        PrevNorm = Solver % Variable % Norm

        Norm = DefaultSolve()

       if (TransientSimulation.and.Compute_dhdt) then
          DHDTSol => VariableGet( Model % Mesh % Variables, 'DHDT',UnFoundFatal=UnFoundFatal)
          DHDT => DHDTSol % Values

          Do i=1,Solver % Mesh % NumberOfNodes
             IF ( DHDTSol % Perm(i) .ge. 1 ) THEN
                DHDT(DHDTSol % Perm(i))=(Solver % Variable % Values(ThickPerm(i))-PreH(ThickPerm(i),1))/dt
             END IF
          End Do
       Endif

        IF ( PrevNorm + Norm /= 0.0_dp ) THEN
           RelativeChange = 2.0_dp * ABS( PrevNorm-Norm ) / (PrevNorm + Norm)
        ELSE
           RelativeChange = 0.0_dp
        END IF

        WRITE( Message, * ) 'Result Norm   : ',Norm
        CALL Info( SolverName, Message, Level=4 )
        WRITE( Message, * ) 'Relative Change : ',RelativeChange
        CALL Info( SolverName, Message, Level=4 )

        !------------------------------------------------------------------------------
        ! compute residual
        !------------------------------------------------------------------------------ 
        SystemMatrix % Values = OldValues
        ForceVector = OldRHS

        IF ( ParEnv % PEs > 1 ) THEN !!!!!!!!!!!!!!!!!!!!!! we have a parallel run
           CALL ParallelInitSolve( SystemMatrix, Thick, ForceVector, ResidualVector )
           CALL ParallelMatrixVector( SystemMatrix, Thick, StiffVector, .TRUE. )
           ResidualVector =  StiffVector - ForceVector
           CALL ParallelSumVector( SystemMatrix, ResidualVector )
        ELSE !!!!!!!!!!!!!!!!!!!!!! serial run 
           CALL CRS_MatrixVectorMultiply( SystemMatrix, Thick, StiffVector)
           ResidualVector =  StiffVector - ForceVector
        END IF
        !-----------------------------
        ! determine "active" nodes set
        !-----------------------------
        IF (ApplyDirichlet) THEN
           numberofsurfacenodes = 0
           DO i=1,Model % NumberOfNodes
              l= ThickPerm(i)  
              IF (l<1) CYCLE
              numberofsurfacenodes = numberofsurfacenodes + 1
              !---------------------------------------------------------
              ! if upper limit is exceeded, manipulate matrix in any case
              !----------------------------------------------------------
              IF ((LimitedSolution(i,1)).AND.(Thick(l)-LowerLimit(i)<0.0_dp )) THEN
                 ActiveNode(i,1) = .TRUE.
              END IF
              IF ((LimitedSolution(i,2)).AND.(Thick(l)-UpperLimit(i)>0.0_dp )) THEN
                 ActiveNode(i,2) = .TRUE.
              END IF

              IF ( LimitedSolution(i,1) .AND. ResidualVector(l) < -LinearTol & 
                       .AND. iter>1 ) ActiveNode(i,1) = .FALSE.
              IF ( LimitedSolution(i,2) .AND. ResidualVector(l) >  LinearTol & 
                       .AND. iter>1 ) ActiveNode(i,2) = .FALSE.

              IF( .NOT.ActiveNode(i,1) .AND. .NOT.ActiveNode(i,2) ) THEN
                 PointerToResidualVector(VarThickResidual % Perm(i)) = 0.0_dp
              ELSE
                 PointerToResidualVector(VarThickResidual % Perm(i)) = ResidualVector(l)
              END IF
            END DO
        END IF
        !------------------------------------------
        ! special treatment for periodic boundaries
        !------------------------------------------

        !------------------------------------------------------------------------------
        ! Relaxation
        !------------------------------------------------------------------------------

        st = CPUTIme()-st
        totat = totat + at
        totst = totst + st

        WRITE(Message,'(a,F8.2,F8.2)') 'Assembly: (s)', at, totat
        CALL Info( SolverName, Message, Level=4 )
        WRITE(Message,'(a,F8.2,F8.2)') ' Solve:    (s)', st, totst
        CALL Info( SolverName, Message, Level=4 )
        !------------------------------------------------------------------------------
        ! write some info on max/min values
        !------------------------------------------------------------------------------
        WRITE(Message,'(a,e13.6,a,e13.6)') &
             'Max/min values Thickness:', MAXVAL(Thick(:)),'/',MINVAL( Thick(:))
        CALL Info(SolverName,Message,Level=4)
        IF (ApplyDirichlet) THEN
           !           WRITE(Message,'(a,i10)') 'Deactivated Periodic BC nodes:', k
           !          CALL Info(SolverName,Message,Level=1)
           WRITE(Message,'(a,i0)') 'Number of surface nodes: ', numberofsurfacenodes
           CALL Info(SolverName,Message,Level=4)
           WRITE(Message,'(a,i0)') 'Number of constrained points (lower limit): ', COUNT(ActiveNode(:,1))
           CALL Info(SolverName,Message,Level=4)
           WRITE(Message,'(a,i0)') 'Number of constrained points (upper limit): ', COUNT(ActiveNode(:,2))
           CALL Info(SolverName,Message,Level=4)
        END IF
        !----------------------
        ! check for convergence
        !----------------------
        IF ( RelativeChange < NonlinearTol ) THEN
           WRITE(Message,'(a,i0,a)') 'Converged after', iter, ' iterations'
           CALL Info(SolverName,Message,Level=4)
           EXIT
        ELSE

        END IF
     END DO ! End loop non-linear iterations
     !------------------------------------------------------------------------------
   CONTAINS

     !------------------------------------------------------------------------------
     !==============================================================================
     SUBROUTINE LocalMatrix( STIFF, MASS, FORCE,&
          LOAD,  Velo, NSDOFs, MeshVelo, &
          Element, nCoord, Nodes, NodeIndexes, &
          TransientSimulation,&
          ALEFormulation)
       !------------------------------------------------------------------------------
       !    INPUT:  LOAD(:)   nodal values of the accumulation/ablation function
       !            
       !            Element         current element
       !            n               number of nodes
       !            Nodes           current node points
       !
       !    OUTPUT: STIFF(:,:)
       !            MASS(:,:)
       !            FORCE(:)
       !------------------------------------------------------------------------------
       !      external variables:
       !      ------------------------------------------------------------------------
       REAL(KIND=dp) ::&
            STIFF(:,:), MASS(:,:), FORCE(:), LOAD(:), &
            Velo(:,:), MeshVelo(:,:)

       INTEGER :: nCoord, NodeIndexes(:), NSDOFs
       TYPE(Nodes_t) :: Nodes
       TYPE(Element_t), POINTER :: Element
       LOGICAL :: TransientSimulation,ALEFormulation

       !------------------------------------------------------------------------------
       !      internal variables:
       !      ------------------------------------------------------------------------
       REAL(KIND=dp) ::&
            Basis(2*nCoord),dBasisdx(2*nCoord,3), &
            Vgauss(3),  Source, &
            X,Y,Z,U,V,W,S,SqrtElementMetric, SU(2*nCoord),SW(2*nCoord),hK,UNorm,divu
       REAL(KIND=dp) :: Tau,Tau1,Tau2
       REAL(KIND=dp),PARAMETER :: r_switch=2.0
       LOGICAL :: TransientStab
       REAL(KIND=dp) :: Tau2_factor
       TYPE(ElementType_t), POINTER :: SaveElementType
       INTEGER :: LinType(2:4) = [202,303,404]

       LOGICAL :: Stat, UseLinear
       INTEGER :: i,j,t,p,q, n
       TYPE(GaussIntegrationPoints_t) :: IntegStuff
       !------------------------------------------------------------------------------

       FORCE = 0.0_dp
       STIFF = 0.0_dp
       MASS  = 0.0_dp

       IF (Bubbles) THEN
          n = nCoord * 2
       ELSE
          n = nCoord
       END IF

       UseLinear = GetLogical( GetSolverParams(), 'Use linear elements', Stat )
       UseLinear = UseLinear .OR. ANY(ActiveNode(NodeIndexes,:))
       UseLinear = UseLinear .AND. Element % TYPE % BasisFunctionDegree==2

       IF ( UseLinear ) THEN
         SaveElementType => Element % TYPE
         Element % TYPE => GetElementType(LinType(GetElementFamily()))
       END IF

       hK = ElementDiameter( Element, Nodes )

       !
       !      Numerical integration:
       !      ----------------------
       IF (Bubbles) THEN
          IntegStuff = GaussPoints( Element, Element % TYPE % gausspoints2)
       ELSE
          IntegStuff = GaussPoints( Element )
       END IF

       SU = 0.0_dp
       SW = 0.0_dp

       DO t = 1,IntegStuff % n
          U = IntegStuff % u(t)
          V = IntegStuff % v(t)
          W = IntegStuff % w(t)
          S = IntegStuff % s(t)
          !
          !        Basis function values & derivatives at the integration point:
          !        -------------------------------------------------------------
          stat = ElementInfo( Element,Nodes,U,V,W,SqrtElementMetric, &
               Basis,dBasisdx, Bubbles=Bubbles )

          !        Correction from metric
          !        ----------------------
          S = S * SqrtElementMetric

          IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
             X = SUM( Nodes % x(1:nCoord) * Basis(1:nCoord) )
             Y = SUM( Nodes % y(1:nCoord) * Basis(1:nCoord) )
             Z = SUM( Nodes % z(1:nCoord) * Basis(1:nCoord) )
             S = S * X
          END IF
          !
          !        Velocities and (norm of) gradient of free surface and source function 
          !        at Gauss point
          !        ---------------------------------------------------------------------

          Vgauss=0.0_dp

          IF (.NOT.ALEFormulation) THEN
             DO i=1,NSDOFs
                Vgauss(i) = SUM( Basis(1:nCoord)*(Velo(i,1:nCoord)))
             END DO
          ELSE
             DO i=1,NSDOFs
                Vgauss(i) = SUM( Basis(1:nCoord)*(Velo(i,1:nCoord) - MeshVelo(i,1:nCoord)))
             END DO
          END IF

          divu = 0.0_dp
          DO i=1,NSDOFs
             divu = divu +  SUM( dBasisdx(1:nCoord,i)*(Velo(i,1:nCoord)))
          END DO

          UNorm = SQRT( SUM( Vgauss(1:NSDOFs)**2 ) )
          Tau=0._dp
          TransientStab=ListGetLogical(GetSolverParams(), 'Transient Stabilisation',Stat )
          IF (.NOT.Stat) TransientStab=.FALSE.
          IF (.NOT.TransientStab) THEN 
            IF (UNorm .NE. 0.0_dp) Tau = hK / ( 2*Unorm )
          ELSE
          ! STAB PARAMETER (SEE Akin and Tezduyar, Calculation of the advective
        ! limit ..., Comput. Methods Appl. Mech. Engrg. 193 (2004)
            !1/Tau1
            Tau1=0._dp
            Do i=1,n
               Tau1=Tau1+ABS(SUM(dBasisdx(i,1:NSDOFS)*Vgauss(1:NSDOFS)))
            End do
            !1/Tau2
            Tau2=0._dp
            Tau2_factor=ListGetConstReal(GetSolverParams(), 'Tau2 factor', Stat )
            IF (.NOT.Stat) Tau2_factor=1.0e-3
            IF (TransientSimulation) Tau2=Tau2_factor*2.0/dt
            !Tau=(1/Tau1^r+1/tau2^r)^(-1/r)
            Tau=Tau1**r_switch+Tau2**r_switch
            IF (Tau.NE.0._dp) Tau=(Tau)**(-1.0/r_switch)
          END IF

          IF ( .NOT. Bubbles ) THEN
             DO p=1,n
                SU(p) = 0.0_dp
                DO i=1,NSDOFs
                   SU(p) = SU(p) + Vgauss(i) * dBasisdx(p,i)
                END DO

                SW(p) = 0.0_dp
                DO i=1,NSDOFs
                   SW(p) = SW(p) + Vgauss(i) * dBasisdx(p,i)
                END DO
             END DO
          END IF

          !        Stiffness matrix:
          !        -----------------
          DO p=1,n
             DO q=1,n
                DO i=1,NSDOFs
                   STIFF(p,q) = STIFF(p,q) + &
                        s * Vgauss(i) * dBasisdx(q,i) * Basis(p)
                END DO
                STIFF(p,q) =  STIFF(p,q) + s * Tau * SU(q) * SW(p)
                STIFF(p,q) =  STIFF(p,q) + s * divu * Basis(q) * (Basis(p) + Tau*SW(p))
             END DO
          END DO


          !        Mass Matrix:
          !        ------------
          IF ( TransientSimulation ) THEN
             DO p=1,n
                DO q=1,n
                   MASS(p,q) = MASS(p,q) +  &
                        S * Basis(q) * (Basis(p) + Tau*SW(p))
                END DO
             END DO
          END IF

          !        Get accumulation/ablation function 
          !        --------------------------------------------------------- 
          Source = 0.0_dp
          Source=SUM(Basis(1:nCoord)*LOAD(1:nCoord))

          !        Assemble force vector:
          !        ---------------------
          FORCE(1:n) = FORCE(1:n) &
               + Source * (Basis(1:n) + Tau*SW(1:n)) * s
       END DO

       IF (UseLinear) THEN
         EdgeMap => GetEdgeMap(GetElementFamily())
         n = ELement % TYPE % NumberOfNodes
         DO i=n+1,n+SIZE(EdgeMap,1)
           j=EdgeMap(i-n,1)
           k=EdgeMap(i-n,2)
           STIFF(i,:) =  0._dp
           STIFF(:,i) =  0._dp
           MASS(i,:)  =  0._dp
           MASS(:,i)  =  0._dp
           STIFF(i,i) =  1._dp
           STIFF(i,j) = -0.5_dp
           STIFF(i,k) = -0.5_dp
           FORCE(i) = 0._dp
           Element % TYPE => SaveElementType
         END DO
       END IF

       !------------------------------------------------------------------------------
     END SUBROUTINE LocalMatrix

     !------------------------------------------------------------------------------
   END SUBROUTINE ThicknessSolver
!------------------------------------------------------------------------------
