!/*****************************************************************************
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
! ******************************************************************************
! *
! *                    Author:       Juha Ruokolainen
! *
! *                    Address: CSC - IT Center for Science Ltd.
! *                                Keilaranta 14, P.O. BOX 405
! *                                  02101 Espoo, Finland
! *                                  Tel. +358 0 457 2723
! *                                Telefax: +358 0 457 2302
! *                              EMail: Juha.Ruokolainen@csc.fi
! *
! *                       Date: 08 Jun 1997
! *
! *                Modified by: Thomas, Mikko, Juha
! *
! *       Date of modification: 19 Dec 2005
! *
! *****************************************************************************/
 
!------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE AdvReactSolverLtd( Model,Solver,dt,TransientSimulation )
!DEC$ATTRIBUTES DLLEXPORT :: AdvReactSolver
!------------------------------------------------------------------------------
!******************************************************************************
!
! Solves the advection-reaction equation with discontinous Galerkin method!
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh, materials, BCs, etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear equation solver options
!
!  DOUBLE PRECISION :: dt,
!     INPUT: Timestep size for time dependent simulations
!
!  LOGICAL :: TransientSimulation
!     INPUT: Steady state or transient simulation
!
!******************************************************************************
     USE DefUtils

     IMPLICIT NONE
!------------------------------------------------------------------------------
! external variables
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     TYPE(Solver_t), TARGET :: Solver
     REAL(KIND=dp) :: dt
     LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER ::&
          BC, BodyForce, Material, Constants, SolverParams
     TYPE(Element_t), POINTER ::&
          Element, Face, &
       ParentElement, LeftParent, RightParent
     LOGICAL ::&
          AllocationsDone = .FALSE., FirstTime = .TRUE., Found, Stat, LimitSolution
     LOGICAL, ALLOCATABLE ::  LimitedSolution(:,:)
     INTEGER ::&
          n1,n2, DIM, k, n, t, istat, i, j, NumberOfFAces, Indexes(128), NonlinearIter, &
          LocalNodes, iter, NumberOfModelNodes,NumberOfMeshNodes, NumberOfBoundaryNodes
     INTEGER, POINTER, DIMENSION(:) :: FlowPerm, NodeIndexes, Perm
     REAL(KIND=dp) :: &
          Norm, PrevNorm, RelativeChange, Relax, RelaxIncrease, RelaxMax, &
          UzawaParameter, UzawaParameterDegrease, UzawaParameterMin, UzawaParameterOrig, &
          NonlinearTol, NonLinearTolMin, NonLinearTolDegrease, OldNonlinearTol, &
          at,at0,totat,st,totst,t1,CPUTime,RealTime
     REAL(KIND=dp), ALLOCATABLE :: &
          MASS(:,:), STIFF(:,:), LOAD(:), &
          FORCE(:), Velocity(:,:), Gamma(:), Ref(:), PrevValues(:), &
          UpperLimit(:),LowerLimit(:),LagrangeMultiplier(:,:)
     REAL(KIND=dp), POINTER, DIMENSION(:) ::&
          FlowValues, Values
     TYPE(Variable_t), POINTER ::&
          Var, ExportVar, FlowVariable
     CHARACTER(LEN=MAX_NAME_LEN)  ::&
          SolverName, VariableName, FlowSolverName, ExportVariable
!*******************************************************************************

     TYPE( Element_t ), POINTER :: Faces(:)
     TYPE(Nodes_t) :: ElementNodes

     INTEGER :: DOFs

     SAVE MASS, STIFF, LOAD, FORCE, Velocity, Gamma, &
          SolverName, VariableName, ExportVariable,  &
          Relax, RelaxIncrease, RelaxMax, &
          UzawaParameter, UzawaParameterDegrease, UzawaParameterMin, UzawaParameterOrig, &
          OldNonlinearTol, NonLinearTolMin,NonLinearTolDegrease, &
          FirstTime, AllocationsDone, LimitSolution, DIM, &
          PrevValues, NumberOfModelNodes, NumberOfMeshNodes, &
          UpperLimit, LowerLimit, LagrangeMultiplier,LimitedSolution
!*******************************************************************************
     !-------------------------------
     ! is there something to be done?
     !-------------------------------
     IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
     
     !---------------------------
     ! stuff done the first time
     !---------------------------
     IF(FirstTime) THEN
        DIM = CoordinateSystemDimension()
        Constants => GetConstants()
        SolverParams => GetSolverParams()

        SolverName = 'AdvReactDGLtd ('// TRIM(Solver % Variable % Name) // ')'
        VariableName = TRIM(Solver % Variable % Name)

        NonlinearIter = GetInteger(   SolverParams, &
             'Nonlinear System Max Iterations', Found )
        IF ( .NOT.Found ) NonlinearIter = 1

        NonlinearTol  = GetConstReal( SolverParams, &
             'Nonlinear System Convergence Tolerance',    Found )

        OldNonlinearTol = NonlinearTol  
        NonlinearTolDegrease = GetConstReal( SolverParams, &
             'Nonlinear System Convergence Tolerance Degrease',    Found )
        IF (.NOT.Found) NonlinearTolDegrease = 1.0D00

        NonlinearTolMin  = GetConstReal( SolverParams, &
             'Nonlinear System Convergence Tolerance Min',    Found )

        Relax = GetConstReal( SolverParams, &
             'Nonlinear System Relaxation Factor',Found )

        IF ( .NOT.Found ) Relax = 1.0D00

        RelaxIncrease = GetConstReal( SolverParams, &
             'Nonlinear System Relaxation Factor Increase',Found )
        IF (.NOT.Found) RelaxIncrease = 1.0D00

        RelaxMax = GetConstReal( SolverParams, &
             'Nonlinear System Relaxation Factor Max',Found )
        IF (.NOT.Found) RelaxMax = Relax     

        UzawaParameter = GetConstReal( SolverParams, &
             'Uzawa Parameter', Found)
        IF (.NOT.Found) THEN                
           CALL INFO(SolverName, 'No Uzawa Parameter found. No limitting of solution', level=1)
           LimitSolution = .FALSE.
           CALL INFO(SolverName, 'Setting nonlinear iterations to 1.')
           NonlinearIter = 1
        ELSE
           WRITE(Message,'(a,e13.4)') 'Uzawa Parameter:',  UzawaParameter
           CALL INFO(SolverName, Message, level=1)
           LimitSolution = .TRUE.
           UzawaParameterOrig = UzawaParameter
        END IF
        IF (LimitSolution) THEN
           UzawaParameterDegrease = GetConstReal( SolverParams, &
                'Uzawa Parameter Degrease', Found)
           IF (.NOT.Found) THEN                
              CALL INFO(SolverName, 'Keyword >Uzawa Parameter Degrease< not found', Level=1)
              CALL INFO(SolverName, 'Setting to 1.',Level=1)
              UzawaParameterDegrease = 1.0
           ELSE IF(UzawaParameterDegrease > 1.0D00) THEN
              CALL INFO(SolverName, '0.0 < Uzawa Parameter Degrease < 1.0!',Level=1)
              CALL INFO(SolverName, 'Setting to 1.',Level=1)
              UzawaParameterDegrease = 1.0
           ELSE
              WRITE(Message,'(a,e13.4)') 'Uzawa Parameter Degrease:',  UzawaParameterDegrease
              CALL INFO(SolverName, Message, level=1)
           END IF
           UzawaParameterMin = GetConstReal( SolverParams, &
                'Uzawa Parameter Min', Found)
           IF (.NOT.Found) THEN                
              CALL INFO(SolverName, 'Keyword >Uzawa Parameter Min< not found', Level=9)
              CALL INFO(SolverName, 'Setting to value of Uzawa Parameter.',Level=1)
              UzawaParameterMin = UzawaParameter
           ELSE IF(UzawaParameterMin .GE. UzawaParameter) THEN
              CALL INFO(SolverName, 'Uzawa Parameter Min < Uzawa Paramter!',Level=1)
              CALL INFO(SolverName, 'Setting to value of Uzawa Parameter.',Level=1)
              UzawaParameterMin = UzawaParameter
           ELSE
              WRITE(Message,'(a,e13.4)') 'Uzawa Parameter Min:',  UzawaParameterMin
              CALL INFO(SolverName, Message, level=1)
           END IF
        END IF
     ELSE
        UzawaParameter = UzawaParameterOrig 
     END IF
     !------------------------------------------------------------------------------
     !    Get variables needed for solution
     !------------------------------------------------------------------------------  
     Var   => VariableGet( Solver % Mesh % Variables, TRIM(VariableName) )
     IF (ASSOCIATED(Var)) THEN
        Perm   => Var % Perm
        Values => Var % Values
        LocalNodes = COUNT( Perm > 0 )
     ELSE
        WRITE(Message,'(A,A)') VariableName, ' not associated.'
        CALL FATAL(Solvername,Message)
     END IF
     !----------------------------------------------
     ! Initialize & allocate some permanent storage
     !----------------------------------------------
     IF ( .NOT.AllocationsDone .OR. Solver % Mesh % Changed) THEN
        NumberOfModelNodes = SIZE(Solver % Variable % Values)
        NumberOfMeshNodes = Model % Mesh % NumberOfNodes
        N = 2 * MAX(Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes )
        IF ( AllocationsDone ) &
             DEALLOCATE (FORCE, MASS, STIFF, LOAD,  &
             Velocity, Gamma, PrevValues,&
             LagrangeMultiplier,   &
             LowerLimit,            &
             UpperLimit,            &
             LimitedSolution)
        ALLOCATE( FORCE(N), MASS(n,n), STIFF(N,N), LOAD(N),         &
             Velocity(3,N), Gamma(n), PrevValues(NumberOfModelNodes),       &
             LagrangeMultiplier(2,NumberOfMeshNodes),             &
             LowerLimit(NumberOfMeshNodes),                        &
             UpperLimit(NumberOfMeshNodes),                        &
             LimitedSolution(2,NumberOfMeshNodes),                 & 
             STAT = istat )        
        IF ( istat /= 0 ) CALL FATAL(SolverName,'Memory allocation error.' )
        AllocationsDone = .TRUE.
        CALL INFO(SolverName,'Allocations done',Level=1)
        ExportVariable = GetString( Solver % Values, 'Exported Variable 1', Found )  
        IF (.NOT.Found) CALL FATAL(SolverName,'Mandatory keyword >Variable Variable 1< not found')
        WRITE(Message,'(A,A)') 'Exported Variable:', ExportVariable
        CALL INFO(SolverName,Message,Level=1)
        LagrangeMultiplier = 0.0d0
     END IF

    !--------------------------------------
    ! determine if interfaces are 1d or 2d
    !--------------------------------------   
     IF ( DIM  == 2 ) THEN
        Faces => Solver % Mesh % Edges
        NumberOfFaces = Solver % Mesh % NumberOfEdges
     ELSE
        Faces => Solver % Mesh % Faces
        NumberOfFaces = Solver % Mesh % NumberOfFaces
     END IF


     !---------------------
     ! Get velocity field
     !---------------------
     FlowSolverName = GetString( Solver % Values, 'Flow Solver Name', Found )    
     IF (.NOT.Found) FlowSolverName = 'flow solution'
     FlowVariable => VariableGet( Solver % Mesh % Variables, FlowSolverName )


     IF ( ASSOCIATED( FlowVariable ) ) THEN
       FlowPerm    => FlowVariable % Perm
       FlowValues  => FlowVariable % Values
     ELSE
       CALL FATAL(SolverName, 'No variable for velocity associated.')
     END IF

     !------------------------------------------------------------
     !  The first time around this has been done by the caller...
     !------------------------------------------------------------
     IF ( TransientSimulation .AND. .NOT.FirstTime ) THEN
        CALL InitializeTimestep( Solver )
     END IF

     !-----------------------
     ! Save current solution
     !-----------------------
 
     PrevValues = Values(1:LocalNodes)

     !----------------------------------
     ! Determine new Nonlinear Tolerance
     !----------------------------------
     IF (.NOT.FirstTime) THEN
        WRITE( Message, * ) 'Nonlinear Tolerance: old:',OldNonLinearTol,&
             ' new:', MAX(NonLinearTolMin,OldNonLinearTol*NonLinearTolDegrease)
        CALL Info( SolverName, Message, Level=4 )
        NonLinearTol = MAX(NonLinearTolMin,OldNonLinearTol*NonLinearTolDegrease)
        OldNonLinearTol = NonLinearTol
     ELSE
        WRITE( Message, * ) 'Nonlinear Tolerance at start:',NonLinearTol
        CALL Info( SolverName, Message, Level=4 )
        OldNonLinearTol = NonLinearTol
     END IF

     !------------------------------------------------------------------------------
     ! non-linear system iteration loop
     !------------------------------------------------------------------------------
     DO iter=1,NonlinearIter
        IF (.NOT.FirstTime) THEN
           WRITE( Message, * ) 'Relaxation Factor: old:',Relax,&
                ' new:', MIN(RelaxMax,Relax * Relaxincrease) 
           CALL Info( SolverName, Message, Level=4 )
           Relax = MIN(RelaxMax,Relax * Relaxincrease)
        ELSE
           WRITE( Message, * ) 'Relaxation Factor:',Relax
           CALL Info( SolverName, Message, Level=4 )
        END IF

        IF (UzawaParameterDegrease < 1.0D00) THEN
           WRITE( Message, * ) 'Uzawa Parameter: ', UzawaParameter
           CALL Info( SolverName, Message, Level=4 )
        END IF

        FirstTime = .FALSE.

        at  = CPUTime()
        at0 = RealTime()


        CALL DefaultInitialize()

        !----------------------------
        ! Update Lagrange multiplier:
        !----------------------------
        IF (LimitSolution) THEN
           LowerLimit = 0.0d00
           UpperLimit = 0.0d00

           LimitedSolution = .TRUE.
           !--------------------------
           ! Lagrange multipliers
           !--------------------------
           DO t=1,Solver % NumberOfActiveElements
              Element => GetActiveElement(t)
              NodeIndexes => Element % Nodeindexes
              !--------------------------
              ! get Indizi for DG-Element
              !--------------------------
              N = GetElementDOFs( Indexes )
              !--------------------------
              ! get number of mesh nodes
              !--------------------------
              N = GetElementNOFNodes(Element)
              CALL GetElementNodes( ElementNodes )
              Material => GetMaterial()
              !----------------
              ! get lower limit
              !----------------
              LowerLimit(Nodeindexes(1:N)) = ListGetReal(Material,TRIM(Solver % Variable % Name) // & 
                   ' Lower Limit',n,NodeIndexes, Found)
              IF (.NOT.Found) THEN
                 LimitedSolution(1,Nodeindexes(1:N)) = .FALSE.
                 WRITE(Message,'(a,i)') 'No lower limit of solution for element no. ', t
                 CALL INFO(SolverName, Message, level=9)
                 LagrangeMultiplier(1,Nodeindexes(1:N)) = 0.0D00
              END IF
              !----------------
              ! get upper limit
              !----------------
              UpperLimit(Nodeindexes(1:N)) = ListGetReal(Material,TRIM(Solver % Variable % Name) // & 
                   ' Upper Limit',n,NodeIndexes, Found)
              IF (.NOT. Found) THEN
                 LimitedSolution(2,Nodeindexes(1:N)) = .FALSE.
                 WRITE(Message,'(a,i)') 'No upper limit of solution for element no. ', t
                 CALL INFO(SolverName, Message, level=9)
                 LagrangeMultiplier(2,Nodeindexes(1:N)) = 0.0D00
              END IF 
              ! --------------------------------
              ! get Lagrange Multipliers for ...
              ! --------------------------------
              DO i=1,N
                 ! ... lower limit
                 !----------------
                 IF (LimitedSolution(1,Nodeindexes(i))) &
                    LagrangeMultiplier(1,Nodeindexes(i)) = MIN( 0.0d0, LagrangeMultiplier(1,Nodeindexes(i)) + UzawaParameter * &
                      (Values(Perm(Indexes(i))) - LowerLimit(NodeIndexes(i))))
                 LimitedSolution(1,Nodeindexes(i)) = .FALSE. 
                 ! ... upper limit
                 !----------------
                 IF (LimitedSolution(2,Nodeindexes(i))) &
                    LagrangeMultiplier(2,Nodeindexes(i)) = MAX( 0.0d0, LagrangeMultiplier(2,Nodeindexes(i)) + UzawaParameter * &
                      (Values(Perm(Indexes(i))) - UpperLimit(Nodeindexes(i))))
                 LimitedSolution(2,Nodeindexes(i)) = .FALSE. 
              END DO
           END DO
        END IF
        CALL INFO(SolverName,'Max/min values Lagrange multiplier ...',Level=4)
        WRITE(Message,'(a,e13.6,a,e13.6)') &
             '... for lower constraint:', MAXVAL( LagrangeMultiplier(1,1:NumberOfMeshNodes) ),&
             '/',MINVAL( LagrangeMultiplier(1,1:NumberOfMeshNodes) )
        CALL INFO(SolverName,Message,Level=4)
        WRITE(Message,'(a,e13.6,a,e13.6)') &
             '... or upper constraint:', MAXVAL( LagrangeMultiplier(2,1:NumberOfMeshNodes) ),&
             '/',MINVAL( LagrangeMultiplier(2,1:NumberOfMeshNodes) )
        CALL INFO(SolverName,Message,Level=4)
        WRITE(Message,'(a,e13.6,a,e13.6)') &
             'Max/min values:', MAXVAL( Values(:)),'/',MINVAL( Values(:))
        CALL INFO(SolverName,Message,Level=4)

        !---------------------------------------------------------------------------
        ! Assembly of the bulk elements:
        !---------------------------------------------------------------------------
        DO t = 1, Solver % NumberOfActiveElements
           Element => GetActiveElement( t )
           NodeIndexes => Element % NodeIndexes 

           n = GetElementNOfNodes( Element )
           N = GetElementDOFs( Indexes )
           
           LOAD = 0.0D00
           !-----------
           ! Body Force
           !-----------
           BodyForce => GetBodyForce( Element )
           LOAD(1:n) = GetReal( BodyForce, TRIM(Solver % Variable % Name) // ' Source', Found )

           IF (.NOT.Found) THEN
              LOAD(1:n) = 0.0D00
              CALL INFO(SolverName,'No Source term found',Level=9)
           END IF

           !---------------------
           ! Lagrange multipliers
           !---------------------
           DO j=1,N
              LOAD(j) = LOAD(j) &
                   - LagrangeMultiplier(1,NodeIndexes(j)) &
                   - LagrangeMultiplier(2,NodeIndexes(j))
           END DO
!           PRINT *,LOAD
        
           k = FlowVariable % DOFs
           Velocity = 0.0d0
           DO i=1,k-1
              DO j=1,N
                 Velocity(i,j) = FlowValues(k*(FlowPerm(NodeIndexes(j))-1)+i)
              END DO
           END DO

           Material => GetMaterial()
           Gamma(1:n)  = GetReal( Material, TRIM(Solver % Variable % Name) // ' Gamma', Found )
           IF (.NOT.Found) Gamma(1:n) = 0.0D00

           CALL LocalMatrix( MASS, STIFF, FORCE, LOAD, Velocity, Gamma, Element, n ) 
           IF ( TransientSimulation ) CALL Default1stOrderTime( MASS, STIFF, FORCE )
           CALL DefaultUpdateEquations( STIFF, FORCE )
        END DO
        !----------------------------
        ! Assembly of the face terms:
        !----------------------------
        FORCE = 0.0d0
        DO t=1,NumberOfFaces
           Face => Faces(t)           
           IF ( .NOT. ActiveBoundaryElement(Face) ) CYCLE
           NodeIndexes => Face % NodeIndexes

           LeftParent  => Face % BoundaryInfo % Left
           RightParent => Face % BoundaryInfo % Right
           IF ( ASSOCIATED(RightParent) .AND. ASSOCIATED(LeftParent) ) THEN
              n  = GetElementNOFNodes( Face )
              n1 = GetElementNOFNodes( LeftParent )
              n2 = GetElementNOFNodes( RightParent )
              
              k = FlowVariable % DOFs
              Velocity = 0.0d0
              DO i=1,k-1  
                 DO j=1,N
                    Velocity(i,j) = FlowValues(k*(FlowPerm(NodeIndexes(j))-1)+i)
                 END DO
              END DO

              CALL LocalJumps( STIFF,Face,n,LeftParent,n1,RightParent,n2,Velocity )
              CALL DefaultUpdateEquations( STIFF, FORCE, Face )
           END IF
        END DO
        !---------------------------------
        ! Loop over the boundary elements:
        !---------------------------------
        DO t=1,Solver % Mesh % NumberOfBoundaryElements
           Element => GetBoundaryElement(t)
           NodeIndexes => Element % NodeIndexes

           IF( .NOT. ActiveBoundaryElement() )  CYCLE
           IF( GetElementFamily(Element) == 1 ) CYCLE

           ParentElement => Element % BoundaryInfo % Left
           IF ( .NOT. ASSOCIATED( ParentElement ) ) &
                ParentElement => Element % BoundaryInfo % Right
        
!           N = GetElementDOFs( Indexes )
           N  = GetElementNOFNodes( Element )
           n1 = GetElementNOFnodes( ParentElement )

           k = FlowVariable % DOFs
           Velocity = 0.0d0
           DO i=1,k-1 
              DO j=1,N
                 Velocity(i,j) = FlowValues(k*(FlowPerm(NodeIndexes(j))-1)+i)
              END DO
           END DO

           BC => GetBC()
           LOAD = 0.0d0
           Found = .FALSE.
           IF ( ASSOCIATED(BC) ) THEN
              LOAD(1:n) = GetReal( BC, Solver % Variable % Name, Found )
              DO j=1,N
                 LOAD(j) = LOAD(j) &
                      - LagrangeMultiplier(1,NodeIndexes(j)) &
                      - LagrangeMultiplier(2,NodeIndexes(j))
              END DO
           END IF

           MASS = 0.0d0
           CALL LocalMatrixBoundary(  STIFF, FORCE, LOAD, &
                Element, n, ParentElement, n1, Velocity, Found )
           
           IF ( TransientSimulation ) CALL Default1stOrderTime( MASS, STIFF, FORCE )
           CALL DefaultUpdateEquations( STIFF, FORCE )
        END DO
        
        CALL DefaultFinishAssembly()
        !---------------------
        ! Solve the system ...
        !---------------------
        at = CPUTime() - at
        st = CPUTime()

        PrevNorm = Solver % Variable % Norm
        Norm = DefaultSolve()
        
        st = CPUTIme()-st
        totat = totat + at
        totst = totst + st
        WRITE(Message,'(a,i4,a,i4,a,F8.2,F8.2)') 'iter: ',iter,'/max. ', NonlinearIter,' Assembly: (s)', at, totat
        CALL Info( SolverName, Message, Level=4 )
        WRITE(Message,'(a,i4,a,i4,a,F8.2,F8.2)') 'iter: ',iter,'/max. ', NonlinearIter,' Solve:    (s)', st, totst
        CALL Info( SolverName, Message, Level=4 )


        !-------------------------------
        ! ... and check for convergence
        !-------------------------------
        IF ( PrevNorm + Norm /= 0.0d0 ) THEN
           RelativeChange = 2.0d0 * ABS( PrevNorm-Norm ) / (PrevNorm + Norm)
        ELSE
           RelativeChange = 0.0d0
        END IF
        
        WRITE( Message, * ) 'Result Norm   : ',Norm
        CALL Info( SolverName, Message, Level=4 )
        WRITE( Message, * ) 'Relative Change : ',RelativeChange
        CALL Info( SolverName, Message, Level=4 )
        
        IF ( RelativeChange < NonlinearTol ) EXIT
        
        UzawaParameter = MAX(UzawaParameter*UzawaParameterDegrease,UzawaParameterMin)
        
        CALL ListAddConstReal(Solver % Values,  &
             'Nonlinear System Relaxation Factor', Relax )
!------------------------------------------------------------------------------
     END DO ! of the nonlinear iteration
!------------------------------------------------------------------------------

     CALL ListAddConstReal(Solver % Values,  &
          'Nonlinear System Tolerance', NonlinearTol)
     !-----------------------------------------------
     ! Average the elemental results to nodal values:
     !-----------------------------------------------
     ExportVar => VariableGet( Solver % Mesh % Variables, TRIM(ExportVariable) )
     IF ( ASSOCIATED( ExportVar ) ) THEN
        n1 = Solver % Mesh % NumberOfNodes
        ALLOCATE( Ref(n1) )
        Ref = 0
           
        IF ( ASSOCIATED( ExportVar % Perm, Solver % Variable % Perm ) ) THEN
           ALLOCATE( ExportVar % Perm(SIZE(Solver % Variable % Perm))  )
           ExportVar % Perm = 0
           DO i = 1,n1
              ExportVar % Perm(i) = i
           END DO
        END IF
           
        ExportVar % Values = 0.0d0
        DO t=1,Solver % NumberOfActiveElements
           Element => GetActiveElement(t) 
           n = GetElementDOFs( Indexes )
           n = GetElementNOFNodes()
           DO i=1,n
              k = Element % NodeIndexes(i)
              ExportVar % Values(k) = ExportVar % Values(k) + &
                   Solver % Variable % Values( Solver % Variable % Perm(Indexes(i)) )
              Ref(k) = Ref(k) + 1
           END DO
        END DO
           
        WHERE( Ref > 0 )
           ExportVar % Values(1:n1) = ExportVar % Values(1:n1) / Ref
        END WHERE
        DEALLOCATE( Ref )
     ELSE
        WRITE(Message,'(A,A,A)') 'Variable >',ExportVariable,'< not found'
        CALL WARN(SolverName,Message)
     END IF

        
   CONTAINS

!------------------------------------------------------------------------------      
     SUBROUTINE LocalMatrix(MASS, STIFF, FORCE, LOAD, Velo, Gamma, Element, n)
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: MASS(:,:), STIFF(:,:), FORCE(:), &
                  LOAD(:), Velo(:,:), Gamma(:)
       INTEGER :: n
       TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: Basis(n),dBasisdx(n,3)
       REAL(KIND=dp) :: detJ,U,V,W,S,A,L,cu(3),g
       LOGICAL :: Stat
       INTEGER :: i,p,q,t,dim
       TYPE(GaussIntegrationPoints_t) :: IntegStuff
       TYPE(Nodes_t) :: Nodes
       SAVE Nodes
!------------------------------------------------------------------------------
       dim = CoordinateSystemDimension()
       FORCE = 0.0d0
       STIFF = 0.0d0
       MASS  = 0.0d0
       CALL GetElementNodes( Nodes, Element )
!------------------------------------------------------------------------------
!      Numerical integration
!------------------------------------------------------------------------------
       IntegStuff = GaussPoints( Element )
 
       DO t=1,IntegStuff % n
         U = IntegStuff % u(t)
         V = IntegStuff % v(t)
         W = IntegStuff % w(t)
         S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
         stat = ElementInfo( Element, Nodes, U, V, W, detJ, &
                 Basis, dBasisdx )

         S = S * detJ
         L = SUM( LOAD(1:n) *  Basis(1:n) )
         g = SUM( Basis(1:n) * Gamma(1:n) )

         cu = 0.0d0
         DO i=1,dim
            cu(i) = SUM( Basis(1:n) * Velo(i,1:n) )
         END DO
!------------------------------------------------------------------------------
!        The advection-reaction equation
!------------------------------------------------------------------------------
         DO p=1,n
            DO q=1,n
              MASS(p,q)  = MASS(p,q)  + s * Basis(q) * Basis(p)
              STIFF(p,q) = STIFF(p,q) + s * g * Basis(q) * Basis(p)
              DO i=1,dim
                STIFF(p,q) = STIFF(p,q) - s * cu(i) * Basis(q) * dBasisdx(p,i)
              END DO
            END DO
         END DO
         FORCE(1:n) = FORCE(1:n) + s*L*Basis(1:n)
!------------------------------------------------------------------------------
      END DO
!------------------------------------------------------------------------------
     END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE FindParentUVW(Face, nFace, Parent, nParent, U, V, W, Basis)
!------------------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(Element_t) :: Face, Parent
      INTEGER :: nFace, nParent
      REAL( KIND=dp ) :: U, V, W, Basis(:)
!------------------------------------------------------------------------------
      INTEGER :: i,j
      REAL(KIND=dp) :: ParentU(nFace), ParentV(nFace), ParentW(nFace)
!------------------------------------------------------------------------------
      DO i = 1,nFace
        DO j = 1,nParent
          IF ( Face % NodeIndexes(i) == Parent % NodeIndexes(j) ) THEN
            ParentU(i) = Parent % Type % NodeU(j)
            ParentV(i) = Parent % Type % NodeV(j)
            ParentW(i) = Parent % Type % NodeW(j)
            EXIT
          END IF
        END DO
      END DO
      U = SUM( Basis(1:nFace) * ParentU(1:nFace) )
      V = SUM( Basis(1:nFace) * ParentV(1:nFace) )
      W = SUM( Basis(1:nFace) * ParentW(1:nFace) )
!------------------------------------------------------------------------------      
    END SUBROUTINE FindParentUVW
!------------------------------------------------------------------------------      


!------------------------------------------------------------------------------
    SUBROUTINE LocalJumps( STIFF,Face,n,LeftParent,n1,RightParent,n2,Velo )
!------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=dp) :: STIFF(:,:), Velo(:,:)
      INTEGER :: n,n1,n2
      TYPE(Element_t), POINTER :: Face, LeftParent, RightParent
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: FaceBasis(n), FacedBasisdx(n,3)
      REAL(KIND=dp) :: LeftBasis(n1), LeftdBasisdx(n1,3)
      REAL(KIND=dp) :: RightBasis(n2), RightdBasisdx(n2,3)
      REAL(KIND=dp) :: Jump(n1+n2), Average(n1+n2)
      REAL(KIND=dp) :: detJ, U, V, W, S, Udotn, xx, yy
      LOGICAL :: Stat
      INTEGER :: i, j, p, q, dim, t, nFace, nParent
      TYPE(GaussIntegrationPoints_t) :: IntegStuff
      REAL(KIND=dp) :: hE, Normal(3), cu(3), LeftOut(3)

      TYPE(Nodes_t) :: FaceNodes, LeftParentNodes, RightParentNodes
      SAVE FaceNodes, LeftParentNodes, RightParentNodes
!------------------------------------------------------------------------------
      dim = CoordinateSystemDimension()
      STIFF = 0.0d0

      CALL GetElementNodes( FaceNodes, Face )
      CALL GetElementNodes( LeftParentNodes,  LeftParent )
      CALL GetElementNodes( RightParentNodes, RightParent )
!------------------------------------------------------------------------------
!     Numerical integration over the edge
!------------------------------------------------------------------------------
      IntegStuff = GaussPoints( Face )

      LeftOut(1) = SUM(LeftParentNodes % x(1:n1)) / n1
      LeftOut(2) = SUM(LeftParentNodes % y(1:n1)) / n1
      LeftOut(3) = SUM(LeftParentNodes % z(1:n1)) / n1
      LeftOut(1) = SUM(FaceNodes % x(1:n)) / n - LeftOut(1)
      LeftOut(2) = SUM(FaceNodes % y(1:n)) / n - LeftOut(2)
      LeftOut(3) = SUM(FaceNodes % z(1:n)) / n - LeftOut(3)

      DO t=1,IntegStuff % n
        U = IntegStuff % u(t)
        V = IntegStuff % v(t)
        W = IntegStuff % w(t)
        S = IntegStuff % s(t)

        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo( Face, FaceNodes, U, V, W, detJ, &
             FaceBasis, FacedBasisdx )

        S = S * detJ

        Normal = NormalVector( Face, FaceNodes, U, V, .FALSE. )
        IF ( SUM( LeftOut*Normal ) < 0 ) Normal = -Normal

        ! Find basis functions for the parent elements:
        ! ---------------------------------------------
        CALL FindParentUVW( Face,n,LeftParent,n1,U,V,W,FaceBasis )
        stat = ElementInfo( LeftParent, LeftParentNodes, U, V, W, detJ, &
                LeftBasis, LeftdBasisdx )

        CALL FindParentUVW( Face,n,RightParent,n2,U,V,W,FaceBasis )
        stat = ElementInfo( RightParent, RightParentNodes, U, V, W, detJ, &
              RightBasis, RightdBasisdx )

        ! Integrate jump terms:
        ! ---------------------
        Jump(1:n1) = LeftBasis(1:n1)
        Jump(n1+1:n1+n2) = -RightBasis(1:n2)

        Average(1:n1) = LeftBasis(1:n1) / 2
        Average(n1+1:n1+n2) = RightBasis(1:n2) / 2

        cu = 0.0d0
        DO i=1,dim
          cu(i) = SUM( Velo(i,1:n) * FaceBasis(1:n) )
        END DO
        Udotn = SUM( Normal * cu )

        DO p=1,n1+n2
          DO q=1,n1+n2
            STIFF(p,q) = STIFF(p,q) + s * Udotn * Average(q) * Jump(p)
            STIFF(p,q) = STIFF(p,q) + s * ABS(Udotn)/2 * Jump(q) * Jump(p)
          END DO
        END DO
      END DO
!------------------------------------------------------------------------------
    END SUBROUTINE LocalJumps
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
    SUBROUTINE LocalMatrixBoundary( STIFF, FORCE, LOAD, &
        Element, n, ParentElement, np, Velo, InFlowBC )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: STIFF(:,:),  FORCE(:), LOAD(:), Velo(:,:)
     INTEGER :: n, np
     LOGICAL :: InFlowBC
     TYPE(Element_t), POINTER :: Element, ParentElement
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(n), dBasisdx(n,3)
     REAL(KIND=dp) :: ParentBasis(np), ParentdBasisdx(np,3)
     INTEGER :: i,j,p,q,t,dim

     REAL(KIND=dp) :: Normal(3), g, L, Udotn, cu(3), cu1(3), detJ,U,V,W,S
     LOGICAL :: Stat, Inflow
     TYPE(GaussIntegrationPoints_t) :: IntegStuff

     TYPE(Nodes_t) :: Nodes, ParentNodes
     SAVE Nodes, ParentNodes
!------------------------------------------------------------------------------
     dim = CoordinateSystemDimension()
     FORCE = 0.0d0
     STIFF = 0.0d0
 
     CALL GetElementNodes( Nodes, Element )
     CALL GetElementNodes( ParentNodes, ParentElement )

     Normal = NormalVector( Element, Nodes, 0.0d0, 0.0d0, .TRUE. ) 
     DO i=1,3
       cu(i) = SUM( Velo(i,1:n) ) / n
     END DO
     Inflow = InFlowBC .AND. SUM( Normal * cu ) < 0.0d0

     ! Numerical integration:
     !-----------------------
     IntegStuff = GaussPoints( Element )

     DO t=1,IntegStuff % n
       U = IntegStuff % u(t)
       V = IntegStuff % v(t)
       W = IntegStuff % w(t)
       S = IntegStuff % s(t)

       Normal = NormalVector( Element, Nodes, U, V, .TRUE. ) 

       ! Basis function values & derivatives at the integration point:
       ! -------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, U, V, W, detJ, &
               Basis, dBasisdx )
       S = S * detJ

       CALL FindParentUVW( Element, n, ParentElement, np, U, V, W, Basis )
       stat = ElementInfo( ParentElement, ParentNodes, U, V, W, &
            detJ, ParentBasis, ParentdBasisdx )

       L = SUM( LOAD(1:n) * Basis(1:n) )
       cu = 0.0d0
       DO i=1,dim
          cu(i)  = SUM( Velo(i,1:n) * Basis(1:n) )
       END DO
       Udotn = SUM( Normal * cu )

       DO p = 1,np
         IF ( Inflow ) THEN
            FORCE(p) = FORCE(p) - s * Udotn*L*ParentBasis(p)
         ELSE
           DO q=1,np
             STIFF(p,q) = STIFF(p,q) + s*Udotn*ParentBasis(q)*ParentBasis(p)
           END DO
         END IF
       END DO
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE LocalMatrixBoundary
!------------------------------------------------------------------------------
 END SUBROUTINE AdvReactSolverLtd
!------------------------------------------------------------------------------
