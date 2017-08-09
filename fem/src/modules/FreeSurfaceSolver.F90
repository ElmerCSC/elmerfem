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
! *  Authors: Thomas Zwinger, Peter R�back, Juha Ruokolainen, Mikko Lyly
! *  Email:   Thomas.Zwinger@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 17 May 2002
! *  Added limiters: 09 Jan 2008
! *
! ****************************************************************************/

!-----------------------------------------------------------------------------
!>  Initialisation routine to set initial free surface height 
!>  when the surface normal is not aligned to Z axis (rotated FS)
!> \ingroup Solvers
!-----------------------------------------------------------------------------
SUBROUTINE FreeSurfaceSolver_RotInit( Model,Solver,dt,TransientSimulation )
  USE DefUtils  
  USE ElementDescription
  IMPLICIT NONE

  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL ::  TransientSimulation
  REAL(KIND=dp) :: dt

  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Variable_t), POINTER :: VarFS, OrientVar

  INTEGER :: DIM, i
  INTEGER, POINTER ::  RefFSPerm(:), FreeSurfPerm(:)

  REAL(KIND=dp), POINTER :: RefFS(:), FreeSurf(:), OrientVarPointer(:), PArray(:,:) => NULL()
  REAL(KIND=dp) :: NodeHolder(3), Orientation(3),  RotationMatrix(3,3)

  LOGICAL :: RotateFS, Found
  CHARACTER(LEN=MAX_NAME_LEN)  :: VariableName, SolverName, OrientVarName

  SolverParams => GetSolverParams()
  DIM = CoordinateSystemDimension()
  SolverName = "FreeSurfaceSolver_RotInit"

  VariableName = ListGetString( SolverParams, 'Free Surface Variable Name', Found)
  IF ( .NOT.Found ) THEN
     CALL Fatal(SolverName, "Can't find 'Free Surface Variable Name'")
  END IF
  CALL ListAddString(Model % Constants, 'Free Surface Variable Name',TRIM(VariableName))

  PArray => ListGetConstRealArray( SolverParams,TRIM(VariableName)//' Orientation', Found)
  IF ( .NOT.Found ) THEN
     WRITE(Message, '(a,a,a)') "Solver to Initialize Rotated Free Surface is present but &
          &can't find ",VariableName," Orientation"
     CALL Fatal(SolverName, Message)
  END IF
  DO i=1,3
     Orientation(i) = PArray(i,1)
  END DO

  !Push to globals
  WRITE(OrientVarName, '(a,a)') TRIM(VariableName),' Orientation' 
  ALLOCATE(OrientVarPointer(3))

  CALL VariableAdd(Model % Mesh % Variables, Model % Mesh, Solver, OrientVarName,3,OrientVarPointer)
  OrientVar => VariableGet(Model %  Mesh % Variables, OrientVarName, .TRUE.)
  IF(.NOT. ASSOCIATED(OrientVar)) &
       CALL FATAL(SolverName, "Internal problem set/getting Orientation Variable")
  OrientVar % Values = Orientation

  RotationMatrix = ComputeRotationMatrix(Orientation)

  !Get Variables
  !--------------------------------------------
  RefFS     => Solver % Variable % Values
  IF (.NOT.ASSOCIATED(RefFS)) CALL Fatal(SolverName,'Variable values not associated')
  RefFSPerm => Solver % Variable % Perm

  VarFS => VariableGet( Model % Mesh % Variables, TRIM(VariableName) )
  IF (.NOT.ASSOCIATED(VarFS)) THEN
     WRITE (Message, '(a,a)') "Can't get variable: ",VariableName
     CALL Fatal( SolverName, Message)
  END IF
  FreeSurf => VarFS % Values
  FreeSurfPerm => VarFS % Perm


  DO i=1,Model % Mesh % NumberOfNodes
     IF(RefFSPerm(i) == 0) CYCLE

     NodeHolder(1) = Model % Mesh % Nodes % x(i)
     NodeHolder(2) = Model % Mesh % Nodes % y(i)
     NodeHolder(3) = Model % Mesh % Nodes % z(i)

     NodeHolder = MATMUL(RotationMatrix, NodeHolder)

     RefFS(RefFSPerm(i)) = NodeHolder(3)
     FreeSurf(FreeSurfPerm(i)) = NodeHolder(3)
  END DO

END SUBROUTINE FreeSurfaceSolver_RotInit

!-----------------------------------------------------------------------------
!>  Functions to compute Mesh Update components from scalar free surface
!>  in cases when free surface is not aligned to z axis (Rotated FS)
!> \ingroup Solvers
!-----------------------------------------------------------------------------
FUNCTION FreeSurfaceToMeshUpdate1( Model, nodenumber,inarray ) RESULT(mu)
  USE Types
  USE Defutils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER :: NodeNumber
  REAL(KIND=dp) ::inarray(*), y, mu, FreeSurfaceToMeshUpdate

  mu = FreeSurfaceToMeshUpdate(Model, nodenumber,inarray, 1)
END FUNCTION

FUNCTION FreeSurfaceToMeshUpdate2( Model, nodenumber,inarray ) RESULT(mu)
  USE Types
  USE Defutils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER :: NodeNumber
  REAL(KIND=dp) ::inarray(*), mu, FreeSurfaceToMeshUpdate

  mu = FreeSurfaceToMeshUpdate(Model, nodenumber,inarray, 2)
END FUNCTION

FUNCTION FreeSurfaceToMeshUpdate3( Model, nodenumber,inarray ) RESULT(mu)
  USE Types
  USE Defutils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER :: NodeNumber
  REAL(KIND=dp) ::inarray(*), mu, FreeSurfaceToMeshUpdate

  mu = FreeSurfaceToMeshUpdate(Model, nodenumber,inarray, 3)
END FUNCTION

FUNCTION FreeSurfaceToMeshUpdate( Model, nodenumber, inarray, axis ) RESULT(mu)
  USE Types
  USE Defutils
  USE ElementDescription
  IMPLICIT NONE

  TYPE(Model_t) :: Model
  TYPE(Variable_t), POINTER :: OrientVar, TimeVar
  INTEGER :: NodeNumber, axis, i
  REAL(KIND=dp) :: t0, FS, RefFS, Diff, inarray(*), mu, UnRotationMatrix(3,3), &
       Orientation(3), DiffVect(3), MuVect(3)
  CHARACTER(LEN=MAX_NAME_LEN)  :: USF_Name, VariableName
  LOGICAL :: Found, FirstTime

  SAVE :: FirstTime, t0
  
  FS = inarray(1)
  RefFS = inarray(2)
  USF_Name = "FreeSurfaceToMeshUpdate"

  TimeVar => VariableGet(Model % Mesh % Variables, 'time' )

  IF(FirstTime) THEN
     IF(.NOT. Found) CALL FATAL(USF_Name, "Can't get time variable")
     t0 = TimeVar % Values(1)
     FirstTime = .FALSE.
  END IF
  
  VariableName = ListGetString( Model % Constants, 'Free Surface Variable Name', Found)
  IF(.NOT. Found) THEN
     IF(t0 == TimeVar % Values(1)) THEN !Isn't set when first called before simulation
        mu = 0
        RETURN
     ELSE
        CALL FATAL(USF_Name, "Couldn't find Free Surface Variable Name")
     END IF
  END IF

  OrientVar => VariableGet(Model % Mesh % Variables, TRIM(VariableName)//' Orientation', Found)
  IF ( .NOT.Found ) THEN
     WRITE(Message, '(a,a,a)') "Can't find ",VariableName," Orientation"
     CALL Fatal(USF_Name, Message)
  END IF

  Orientation = OrientVar % Values
  UnRotationMatrix = TRANSPOSE(ComputeRotationMatrix(Orientation))

  Diff = FS - RefFS
  DiffVect = [0.0_dp, 0.0_dp, Diff]

  MuVect = MATMUL(UnRotationMatrix, DiffVect)
  mu = MuVect(axis)

END FUNCTION FreeSurfaceToMeshUpdate

!-----------------------------------------------------------------------------
!>  Solver for free surface evolution in 2d and 3d flows
!>  with or without surface flux, and upper and lower limiters.
!> \ingroup Solvers
!-----------------------------------------------------------------------------
SUBROUTINE FreeSurfaceSolver( Model,Solver,dt,TransientSimulation )
  USE DefUtils
  USE ElementDescription
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
       NeedOldValues, LimitDisp,  Bubbles = .TRUE.,&
       NormalFlux = .TRUE., SubstantialSurface = .TRUE.,&
       UseBodyForce = .TRUE., ApplyDirichlet=.FALSE.,  ALEFormulation=.FALSE.,&
       RotateFS, ReAllocate=.TRUE., ResetLimiters=.FALSE.
  LOGICAL, ALLOCATABLE ::  LimitedSolution(:,:), ActiveNode(:,:)

  INTEGER :: & 
       i,j,K,L, p, q, R, t,N,NMAX,MMAX,nfamily, deg, Nmatrix,&
       edge, bf_id,DIM,istat,LocalNodes,nocorr,&
       NSDOFs,NonlinearIter,iter, numberofsurfacenodes
  INTEGER, POINTER ::&
       FreeSurfPerm(:), FlowPerm(:), NodeIndexes(:), EdgeMap(:,:)

#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: &
       at,st,totat,totst,Norm,PrevNorm,LocalBottom, cv, &
       Relax, MaxDisp, maxdh, maxdh_comm, LinearTol,NonlinearTol,RelativeChange,&
       smallestpossiblenumber, rr, ss, Orientation(3), RotationMatrix(3,3),&
       NodeHolder(3)
#else
  REAL(KIND=dp) :: &
       at,st,totat,totst,CPUTime,Norm,PrevNorm,LocalBottom, cv, &
       Relax, MaxDisp, maxdh, maxdh_comm, LinearTol,NonlinearTol,RelativeChange,&
       smallestpossiblenumber, rr, ss, Orientation(3), RotationMatrix(3,3),&
       NodeHolder(3)
#endif

  REAL(KIND=dp), POINTER :: ForceVector(:), FreeSurf(:), PreFreeSurf(:,:), &
       FlowSolution(:), PrevFlowSol(:,:), PointerToResidualVector(:)

  REAL(KIND=dp), ALLOCATABLE :: ResidualVector(:), &
       STIFF(:,:),SourceFunc(:),FORCE(:), TimeForce(:), &
       MASS(:,:), Velo(:,:), Flux(:,:), LowerLimit(:), UpperLimit(:), &
       OldValues(:), OldRHS(:),StiffVector(:),MeshVelocity(:,:), ElemFreeSurf(:)

  CHARACTER(LEN=MAX_NAME_LEN)  :: SolverName, VariableName, FlowSolName,ConvectionFlag, StabilizeFlag

  TYPE(Nodes_t)   :: ElementNodes
  TYPE(Element_t),POINTER :: CurrentElement
  TYPE(Variable_t), POINTER :: FlowSol, VarSurfResidual, RotVar
  TYPE(ValueList_t), POINTER :: BodyForce, SolverParams, Material, Equation
  TYPE(Matrix_t), POINTER :: Systemmatrix
  !-----------------------------------------------------------------------------
  !      remember these variables
  !----------------------------------------------------------------------------- 
  SAVE STIFF, MASS, SourceFunc, FORCE, &
       ElementNodes, AllocationsDone, ReAllocate, Velo, TimeForce, &
       ElemFreeSurf, Flux, SubstantialSurface, NormalFlux,&
       UseBodyForce, LimitedSolution, LowerLimit, &
       UpperLimit, ActiveNode, ResetLimiters, OldValues, OldRHS, &
       ResidualVector, StiffVector, MeshVelocity
  !------------------------------------------------------------------------------
  !    Get variables for the solution
  !------------------------------------------------------------------------------
  FreeSurf     => Solver % Variable % Values     ! Nodal values for free surface displacement
  IF (.NOT.ASSOCIATED(FreeSurf)) CALL Fatal(SolverName,'Variable values not associated')
  FreeSurfPerm => Solver % Variable % Perm       ! Permutations for free surface displacement
  PreFreeSurf  => Solver % Variable % PrevValues ! Nodal values for free surface displacement
  !------------------------------------------------------------------------------
  !    Get variabel/solver name
  !------------------------------------------------------------------------------
  IF (VariableName .NE. TRIM(Solver % Variable % Name)) THEN
    VariableName = TRIM(Solver % Variable % Name)
    ReAllocate = .TRUE.
  ELSE
    Reallocate = .FALSE.
  END IF
    
  SolverName = 'FreeSurfaceSolver ('// TRIM(Solver % Variable % Name) // ')'
  
  !------------------------------------------------------------------------------
  !    if this partition (or the serial problem) has no free surface,
  !    then nothing to be doneGet variabel/solver name
  !------------------------------------------------------------------------------
  IF ( COUNT(FreeSurfPerm/=0)==0) THEN
     IF (ParEnv % PEs > 1) THEN
        WRITE(Message,'(A,I0,A)')  'Partition ', ParEnv % myPE, ' has no free surface'
        CALL Warn(SolverName,Message)
        CALL Warn(SolverName,'This is not good for load balance!')
     ELSE
        CALL Warn(SolverName,'A serial run without a free surface, but the solver switched in - weird!')
     END IF
     RETURN
  END IF
  !------------------------------------------------------------------------------
  !    Get constants and solver params
  !------------------------------------------------------------------------------
  DIM = CoordinateSystemDimension()
  smallestpossiblenumber = TINY(smallestpossiblenumber)
  SolverParams => GetSolverParams()

  SystemMatrix => Solver % Matrix
  ForceVector => Solver % Matrix % RHS

  cv = GetConstReal( SolverParams, 'Velocity Implicity', Found)
  IF(.NOT. Found) cv = 1.0_dp 
  WRITE(Message,'(a,F9.2)') 'Velocity implicity (1=fully implicit)=',cv
  CALL Info(SolverName, Message, Level=6 )

  LinearTol = GetConstReal( SolverParams, &
       'Linear System Convergence Tolerance',    Found )
  IF ( .NOT.Found ) THEN
     CALL Fatal(SolverName, 'No > Linear System Convergence Tolerance < found')
  END IF
  NonlinearTol  = GetConstReal( SolverParams, &
       'Nonlinear System Convergence Tolerance',    Found )
  NonlinearIter = GetInteger(   SolverParams, &
       'Nonlinear System Max Iterations', Found )
  IF ( .NOT.Found ) NonlinearIter = 1

  MaxDisp = GetConstReal( SolverParams, 'Maximum Displacement', LimitDisp)

  Relax = GetCReal( SolverParams, 'Relaxation Factor', Found)
  IF(.NOT. Found) Relax = 1.0_dp
  NeedOldValues = (Found .AND. (Relax < 1.0_dp)) .OR. LimitDisp 

  ApplyDirichlet = GetLogical( SolverParams, &
       'Apply Dirichlet', Found)
  IF ( .NOT.Found ) THEN
    ApplyDirichlet = .FALSE.
    CALL Info(SolverName, 'No keyword > Apply Dirichlet < found. No limitation of solution',Level=6 )
  ELSE
    IF (ApplyDirichlet) THEN
      CALL Info(SolverName, 'Using Dirichlet method for limitation',Level=6 )
      ResetLimiters = GetLogical( SolverParams, &
           'Reset Limiter', Found)
      IF (.NOT.Found) THEN
        ResetLimiters = .FALSE.
      ELSE
        CALL INFO(SolverName,"Limiters will be reset for each nonlinear iteration",Level=3)
      END IF
      IF (NonlinearIter < 2) THEN
        CALL Warn(SolverName, 'Keyword > Apply Dirichlet < set, but > Nonlinear System Max Iterations < set to lower than 2')
      END IF
    ELSE
      CALL Info(SolverName, 'No limitation of solution',Level=6 )
    END IF
  END IF

  ALEFormulation = GetLogical( SolverParams, &
       'ALE Formulation', Found)
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
     CALL Info(SolverName, &
          'Using residual squared-stabilized formulation.',Level=6 )
  CASE('bubbles')
     Bubbles = .TRUE.
     CALL Info(SolverName, 'Using residual free bubble stabilization',Level=6 )
  END SELECT

  RotVar => VariableGet( Model % Mesh % Variables,TRIM(VariableName)//' Orientation')
  IF(ASSOCIATED(RotVar)) THEN
     RotateFS = .TRUE.
     DO i=1,3
        Orientation(i) = RotVar % Values(i)
     END DO

     RotationMatrix = ComputeRotationMatrix(Orientation)

     WRITE(Message,'(A,f8.2,f8.2,f8.2)') 'Rotated Free Surface defined using vector: ',Orientation
     CALL Info(SolverName, Message,Level=6 )
  ELSE 
     RotateFS = .FALSE.
     CALL Info(SolverName, 'No Free Surface Orientation Vector found, assuming normal to z-axis',&
          Level=6 )
  END IF

  WRITE(Message,'(A,I0)') 'Mesh dimension: ', DIM
  CALL Info( SolverName, Message, Level=8 )

  !------------------------------------------------------------------------------
  !    Allocate some permanent storage, this is done first time only
  !------------------------------------------------------------------------------

  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed .OR. ReAllocate) THEN
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
           Velo,  &
           MeshVelocity, &
           Flux, &
           ElemFreeSurf,&
           SourceFunc )

      IF( ApplyDirichlet ) THEN
        DEALLOCATE( LowerLimit,                      &
             UpperLimit, &
             LimitedSolution,  &
             ActiveNode,                      & 
             ResidualVector, &
             StiffVector,  &
             OldValues, &
             OldRHS)
      END IF

      !IF( NeedOldValues ) DEALLOCATE( OldFreeSurf ) 
    END IF


    IF (Bubbles) THEN
      Nmatrix = 2*NMAX
    ELSE
      Nmatrix = NMAX
    END IF

    ALLOCATE(  ElementNodes % x( NMAX ),    &
         ElementNodes % y( NMAX ),    &
         ElementNodes % z( NMAX ),    &
         TimeForce( Nmatrix ),        &
         FORCE( Nmatrix ),    &
         STIFF( Nmatrix, Nmatrix ), &
         MASS( Nmatrix, Nmatrix ),  &
         Velo( 3, NMAX ), &
         MeshVelocity( 3,NMAX ), &
         Flux( 3, NMAX), &
         ElemFreeSurf( NMAX ),&
         SourceFunc( NMAX ), &
         STAT=istat )
    IF ( istat /= 0 ) THEN
      CALL Fatal(SolverName,'Memory allocation error 1, Aborting.')
    END IF

    ElemFreeSurf = 0._dp

    IF(NeedOldValues) THEN
      !ALLOCATE(OldFreeSurf(SIZE(FreeSurf)), STAT=istat)

      IF ( istat /= 0 ) THEN
        CALL Fatal(SolverName,'Memory allocation error 2, Aborting.')
      END IF
    END IF

    IF( ApplyDirichlet ) THEN
      ALLOCATE( LowerLimit( MMAX ), &
           UpperLimit( MMAX ), &
           LimitedSolution( MMAX, 2 ),  &
           ActiveNode( MMAX, 2 ),                      &  
           ResidualVector( L ),                    &
           StiffVector( L ), &
           OldValues( K ), &
           OldRHS( L ), &
           STAT=istat )
      IF ( istat /= 0 ) THEN
        CALL Fatal(SolverName,'Memory allocation error 3, Aborting.')
      END IF
      ActiveNode = .FALSE.
      ResidualVector = 0.0_dp
    END IF

    CALL Info(SolverName,'Memory allocations done' )
    AllocationsDone = .TRUE.
  END IF


  !   from previous timestep
  !IF( NeedOldValues) THEN
  !   OldFreeSurf = FreeSurf
  !END IF

  !------------------------------------------------------------------------------
  !    Get variables for the residual
  !------------------------------------------------------------------------------
  IF( ApplyDirichlet ) THEN
    VarSurfResidual => VariableGet( Model % Mesh % Variables, TRIM(VariableName) // ' Residual' )
    IF (.NOT.ASSOCIATED(VarSurfResidual)) THEN
      WRITE(Message,'(A)') '>' // TRIM(VariableName) // ' Residual < not associated'
      CALL Fatal( SolverName, Message)
    END IF
    PointerToResidualVector => VarSurfResidual % Values
  END IF

  IF (ResetLimiters)  ActiveNode = .FALSE.
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
        n = GetElementNOFNodes()
        IF(GetElementFamily() == 1) CYCLE
        NodeIndexes => CurrentElement % NodeIndexes

        ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
        ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
        ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)

        ! rotate element nodes round to reference direction
        !------------------------------------------------------------------------------
        IF(RotateFS) THEN
           DO i=1,n
              NodeHolder(1) = ElementNodes % x(i)
              NodeHolder(2) = ElementNodes % y(i)
              NodeHolder(3) = ElementNodes % z(i)

              NodeHolder = MATMUL(RotationMatrix,NodeHolder)
              
              ElementNodes % x(i) = NodeHolder(1)
              ElementNodes % y(i) = NodeHolder(2)
              ElementNodes % z(i) = NodeHolder(3)
           END DO
        END IF

        ! set coords of highest occurring dimension to zero (to get correct path element)
        !-------------------------------------------------------------------------------
        ElementNodes % z(1:n) = 0.0_dp
        IF (DIM == 2) THEN
           ElementNodes % y(1:n) = 0.0
        ELSE IF(DIM .NE. 3) THEN
           WRITE(Message,'(a,i0,a)')&
                'It is not possible to compute free-surface problems in DIM=',&
                DIM, ' dimensions. Aborting'
           CALL Fatal( SolverName, Message) 
           STOP   
        END IF

        ! get pointers on Equation, Material and body-Force section input
        !----------------------------------------------------------------
        Equation => GetEquation()
        Material => GetMaterial()
        BodyForce => GetBodyForce()

        IF( ApplyDirichlet ) THEN
          ! get lower limit for solution 
          !-----------------------------
          LowerLimit(CurrentElement % Nodeindexes(1:N)) = &
              ListGetReal(Material,'Min ' // TRIM(VariableName),n,CurrentElement % NodeIndexes, Found)
          LimitedSolution(CurrentElement % Nodeindexes(1:N), 1) = Found
          ! get upper limit for solution 
          !-----------------------------
          UpperLimit(CurrentElement % Nodeindexes(1:N)) = &
              ListGetReal(Material,'Max ' // TRIM(VariableName),n,CurrentElement % NodeIndexes, Found)              
          LimitedSolution(CurrentElement % Nodeindexes(1:N), 2) = Found
        END IF

        ! get flow soulution and velocity field from it
        !----------------------------------------------
        ConvectionFlag = GetString( Equation, 'Convection', Found )
        IF (.NOT. Found) &
             CALL Fatal(SolverName, 'No string for keyword > Convection < found in Equation')
        Velo = 0.0_dp
        ! constant (i.e., in section Material given) velocity
        !----------------------------------------------------
        IF ( ConvectionFlag == 'constant' ) THEN
           Velo(1,1:N) = GetReal( Material, 'Convection Velocity 1', Found )
           IF ( .NOT.Found ) &
                Velo(1,1:N) = GetReal( Equation, 'Convection Velocity 1', Found )

           Velo(2,1:N) = GetReal( Material, 'Convection Velocity 2', Found )
           IF ( .NOT.Found ) &
                Velo(2,1:N) = GetReal( Equation, 'Convection Velocity 2', Found )

           Velo(3,1:N) = GetReal( Material, 'Convection Velocity 3', Found )
           IF ( .NOT.Found ) &
                Velo(3,1:N) = GetReal( Equation, 'Convection Velocity 3', Found )
           ! computed velocity
           !------------------
        ELSE IF (ConvectionFlag == 'computed' ) THEN
           FlowSolName =  GetString( Equation,'Flow Solution Name', Found)
           IF(.NOT.Found) THEN        
              CALL Warn(SolverName,'Keyword > Flow Solution Name < not found in section >Equation<')
              CALL Warn(SolverName,'Taking default value > Flow Solution <')
           END IF
           FlowSol => VariableGet( Solver % Mesh % Variables, FlowSolName )
           IF ( ASSOCIATED( FlowSol ) ) THEN
              FlowPerm     => FlowSol % Perm
              NSDOFs     =  FlowSol % DOFs
              FlowSolution => FlowSol % Values
              PrevFlowSol => FlowSol % PrevValues
           ELSE
              WRITE(Message,'(A,A,A)') &
                   'Convection flag set to > computed <, but no variable >',FlowSolName,'< found'
              CALL Fatal(SolverName,Message)              
           END IF
           ! get velocity profile
           IF ( ASSOCIATED( FlowSol ) ) THEN
              DO i=1,n
                 j = NSDOFs*FlowPerm(NodeIndexes(i))

                 IF(TransientSimulation .AND. ABS(cv-1.0) > 0.001) THEN
                    IF((DIM == 2) .AND. (NSDOFs == 3)) THEN
                       Velo(1,i) = cv * FlowSolution( j-2 ) + (1-cv) * PrevFlowSol(j-2,1)
                       Velo(2,i) = cv * FlowSolution( j-1 ) + (1-cv) * PrevFlowSol(j-1,1)
                       Velo(3,i) = 0.0_dp
                    ELSE IF ((DIM == 3) .AND. (NSDOFs == 4)) THEN
                       Velo(1,i) = cv * FlowSolution( j-3 ) + (1-cv) * PrevFlowSol(j-3,1)
                       Velo(2,i) = cv * FlowSolution( j-2 ) + (1-cv) * PrevFlowSol(j-2,1)
                       Velo(3,i) = cv * FlowSolution( j-1 ) + (1-cv) * PrevFlowSol(j-1,1)
                    ELSE IF ((CurrentCoordinateSystem() == CylindricSymmetric) &
                         .AND. (DIM == 2) .AND. (NSDOFs == 4)) THEN  
                       Velo(1,i) = cv * FlowSolution( j-3 ) + (1-cv) * PrevFlowSol(j-3,1)
                       Velo(2,i) = cv * FlowSolution( j-2 ) + (1-cv) * PrevFlowSol(j-2,1)
                       Velo(3,i) = cv * FlowSolution( j-1 ) + (1-cv) * PrevFlowSol(j-1,1)
                    ELSE
                       WRITE(Message,'(a,i0,a,i0,a)')&
                            'DIM=', DIM, ' NSDOFs=', NSDOFs, ' does not combine. Aborting'
                       CALL Fatal( SolverName, Message)               
                    END IF
                 ELSE
                    IF((DIM == 2) .AND. (NSDOFs == 3)) THEN
                       Velo(1,i) = FlowSolution( j-2 ) 
                       Velo(2,i) = FlowSolution( j-1 ) 
                       Velo(3,i) = 0.0_dp
                    ELSE IF ((DIM == 3) .AND. (NSDOFs == 4)) THEN
                       Velo(1,i) = FlowSolution( j-3 ) 
                       Velo(2,i) = FlowSolution( j-2 ) 
                       Velo(3,i) = FlowSolution( j-1 ) 
                    ELSE IF ((CurrentCoordinateSystem() == CylindricSymmetric) &
                         .AND. (DIM == 2) .AND. (NSDOFs == 4)) THEN
                       Velo(1,i) = FlowSolution( j-3 ) 
                       Velo(2,i) = FlowSolution( j-2 ) 
                       Velo(3,i) = FlowSolution( j-1 ) 
                    ELSE
                       WRITE(Message,'(a,i0,a,i0,a)')&
                            'DIM=', DIM, ' NSDOFs=', NSDOFs, ' does not combine. Aborting'
                       CALL Fatal( SolverName, Message)
                    END IF
                 END IF
              END DO
           ELSE
              Velo=0.0_dp          
           END IF
        ELSE
           Velo=0.0_dp  
        END IF

        IF(RotateFS) THEN
           DO i=1,n
              Velo(:,i) = MATMUL(RotationMatrix,Velo(:,i))
           END DO
        END IF
        !------------------------------------------------------------------------------
        ! Get mesh velocity
        !------------------------------------------------------------------------------
        MeshVelocity = 0.0_dp
        CALL GetVectorLocalSolution( MeshVelocity, 'Mesh Velocity',CurrentElement)
        !------------------------------------------------------------------------------
        !      get the accumulation/ablation rate (i.e. normal surface flux)
        !      from the body force section
        !------------------------------------------------------------------------------
        SourceFunc = 0.0_dp
        Flux  = 0.0_dp
        SubstantialSurface = .TRUE.

        IF (ASSOCIATED( BodyForce ) ) THEN
           SubstantialSurface = .FALSE.
           ! Accumulation/ablation is given in normal direction of surface:
           !---------------------------------------------------------------
           SourceFunc(1:n) = GetReal( BodyForce, &
                TRIM(VariableName) // ' Accumulation', NormalFlux ) 
           ! Accumulation/ablation has to be computed from given flux:
           !----------------------------------------------------------
           IF (.NOT.NormalFlux) THEN
              Flux(1,1:n) = GetReal( BodyForce, TRIM(VariableName) // ' Accumulation Flux 1',Found)
              IF (.NOT.Found) Flux(1,1:n) = 0.0_dp
              IF (DIM >= 2) THEN
                 Flux(2,1:n) = GetReal( BodyForce, TRIM(VariableName) // ' Accumulation Flux 2',Found )
                 IF (.NOT.Found) Flux(2,1:n) = 0.0_dp
              ELSE
                 Flux(2,1:n) = 0.0_dp
              END IF
              IF (DIM == 3) THEN
                 Flux(3,1:n) = GetReal( BodyForce, TRIM(VariableName) // ' Accumulation Flux 3',Found )
                 IF (.NOT.Found) Flux(3,1:n) = 0.0_dp
              ELSE
                 Flux(3,1:n) = 0.0_dp
              END IF
              SourceFunc = 0.0_dp

              IF(RotateFS) THEN
                 DO i=1,n
                    Flux(:,i) = MATMUL(RotationMatrix, Flux(:,i))
                 END DO
              END IF
           END IF
        END IF

        IF( TransientSimulation) THEN
           ElemFreeSurf(1:n) = PreFreeSurf(FreeSurfPerm(NodeIndexes),1)
        END IF

        !------------------------------------------------------------------------------
        !      Get element local matrix, and rhs vector
        !------------------------------------------------------------------------------
        CALL LocalMatrix( STIFF, MASS, FORCE,&
             SourceFunc, ElemFreeSurf, Velo, MeshVelocity, CurrentElement,&
             n, ElementNodes, NodeIndexes, TransientSimulation,&
             Flux, NormalFlux, SubstantialSurface, ALEFormulation)

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

     IF (ApplyDirichlet) THEN

       OldValues = SystemMatrix % Values
       OldRHS = ForceVector
       
       ! manipulation of the matrix
       !---------------------------
       DO i=1,Model % Mesh % NumberOfNodes
         k = FreeSurfPerm(i)
         IF ((ActiveNode(i,1) .AND. ActiveNode(i,2))) &
              CALL FATAL(SolverName,"Upper as well as lower limiter active - this is a deadlock")
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
     IF( ApplyDirichlet ) THEN
       SystemMatrix % Values = OldValues
       ForceVector = OldRHS
       
       IF ( ParEnv % PEs > 1 ) THEN !!!!!!!!!!!!!!!!!!!!!! we have a parallel run
         CALL ParallelInitSolve( SystemMatrix, FreeSurf, ForceVector, ResidualVector )
         CALL ParallelMatrixVector( SystemMatrix, FreeSurf, StiffVector, .TRUE. )
         ResidualVector =  StiffVector - ForceVector
         CALL ParallelSumVector( SystemMatrix, ResidualVector )
       ELSE !!!!!!!!!!!!!!!!!!!!!! serial run 
         CALL CRS_MatrixVectorMultiply( SystemMatrix, FreeSurf, StiffVector)
         ResidualVector =  StiffVector - ForceVector
       END IF
       !-----------------------------
       ! determine "active" nodes set
       !-----------------------------
       numberofsurfacenodes = 0
       DO i=1,Model % NumberOfNodes
         l= FreeSurfPerm(i)  
         IF (l<1) CYCLE
         numberofsurfacenodes = numberofsurfacenodes + 1
         !---------------------------------------------------------
         ! if upper limit is exceeded, manipulate matrix in any case
         !----------------------------------------------------------
         IF ((LimitedSolution(i,1)).AND.(FreeSurf(l)-LowerLimit(i)<0.0_dp )) THEN
           ActiveNode(i,1) = .TRUE.
         END IF
         IF ((LimitedSolution(i,2)).AND.(FreeSurf(l)-UpperLimit(i)>0.0_dp )) THEN
           ActiveNode(i,2) = .TRUE.
         END IF
         
         IF ( LimitedSolution(i,1) .AND. ResidualVector(l) < -LinearTol & 
             .AND. iter>1 ) ActiveNode(i,1) = .FALSE.
         IF ( LimitedSolution(i,2) .AND. ResidualVector(l) >  LinearTol & 
             .AND. iter>1 ) ActiveNode(i,2) = .FALSE.
         
         IF( .NOT.ActiveNode(i,1) .AND. .NOT.ActiveNode(i,2) ) THEN
           PointerToResidualVector(VarSurfResidual % Perm(i)) = 0.0_dp
         ELSE
           PointerToResidualVector(VarSurfResidual % Perm(i)) = ResidualVector(l)
         END IF
       END DO
     END IF
     !------------------------------------------
     ! special treatment for periodic boundaries
     !------------------------------------------

     !------------------------------------------------------------------------------
     ! Relaxation
     !------------------------------------------------------------------------------
     IF(NeedOldValues) THEN
       IF(LimitDisp) THEN 
         maxdh = -HUGE(maxdh)         
         DO i=1, Model % NumberOfNodes
           j = FreeSurfPerm(i)
           IF(j > 0) THEN
             maxdh = MAX(maxdh, ABS(FreeSurf(j)-PreFreeSurf(j,1)))
           END IF
         END DO
         maxdh = ParallelReduction(maxdh,2)
         IF(maxdh > MaxDisp) THEN
           Relax = Relax * MaxDisp/maxdh
         END IF
         WRITE(Message,'(a,E9.2)') 'Maximum displacement ',maxdh
         CALL Info( SolverName, Message, Level=4 )
       END IF
       WRITE(Message,'(a,F9.2)') 'pp Relaxation factor',Relax
       CALL Info( SolverName, Message, Level=4 )
       DO i=1, Model % NumberOfNodes
         j = FreeSurfPerm(i)
         IF(j > 0) THEN
           FreeSurf(j) = Relax * FreeSurf(j) + (1-Relax) * PreFreeSurf(j,1)
         END IF
       END DO
     END IF
     
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
         'Max/min values surface:', MAXVAL(FreeSurf(:)),'/',MINVAL( FreeSurf(:))
     CALL Info(SolverName,Message,Level=4)

     IF (ApplyDirichlet) THEN
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
     END IF
   END DO ! End loop non-linear iterations
     !------------------------------------------------------------------------------
   CONTAINS

     !------------------------------------------------------------------------------
     !==============================================================================
     SUBROUTINE LocalMatrix( STIFF, MASS, FORCE,&
          SourceFunc, OldFreeSurf, Velo, MeshVelo, &
          Element, nCoord, Nodes, NodeIndexes, TransientSimulation,&
          Flux, NormalFlux, SubstantialSurface, ALEFormulation)
       !------------------------------------------------------------------------------
       !    INPUT:  SourceFunc(:)   nodal values of the accumulation/ablation function
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
            STIFF(:,:), MASS(:,:), FORCE(:), SourceFunc(:), &
            Velo(:,:), MeshVelo(:,:), OldFreeSurf(:), Flux(:,:)

       INTEGER :: nCoord, NodeIndexes(:)
       TYPE(Nodes_t) :: Nodes
       TYPE(Element_t), POINTER :: Element
       LOGICAL :: TransientSimulation,NormalFlux,SubstantialSurface,ALEFormulation
       !------------------------------------------------------------------------------
       !      internal variables:
       !      ------------------------------------------------------------------------
       REAL(KIND=dp) ::&
            Basis(2*nCoord),dBasisdx(2*nCoord,3), &
            Vgauss(3), VMeshGauss(3), Source, gradFreeSurf(3), normGradFreeSurf,&
            FluxGauss(3),X,Y,Z,U,V,W,S,SqrtElementMetric, SU(2*nCoord),SW(2*nCoord),Tau,hK,UNorm

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

       IF( ApplyDirichlet ) THEN
         UseLinear = UseLinear .OR. ANY(ActiveNode(NodeIndexes,:))
       END  IF

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

          gradFreeSurf=0.0_dp
          Vgauss=0.0_dp
          VMeshGauss=0.0_dp

          DO i=1,DIM-1
             gradFreeSurf(i) = SUM(dBasisdx(1:nCoord,i)*OldFreeSurf(1:nCoord))
          END DO

          gradFreeSurf(DIM) = 1.0_dp
          
          IF (.NOT.ALEFormulation) THEN
             DO i=1,DIM
                Vgauss(i) = SUM( Basis(1:nCoord)*(Velo(i,1:nCoord)))
             END DO
          ELSE
             MeshVelo(DIM,1:nCoord) = 0.0_dp
             DO i=1,DIM
                Vgauss(i) = SUM( Basis(1:nCoord)*(Velo(i,1:nCoord) - MeshVelo(i,1:nCoord)))
             END DO
          END IF

          IF (DIM==3) THEN
             normGradFreeSurf = SQRT(1.0_dp + gradFreeSurf(1)**2 + &
                  gradFreeSurf(2)**2)
          ELSE
             normGradFreeSurf = SQRT(1.0_dp + gradFreeSurf(1)**2)
          END IF

          UNorm = SQRT( SUM( Vgauss(1:dim-1)**2 ) )
          IF (UNorm .NE. 0.0_dp) THEN
             Tau = hK / ( 2*Unorm )
          ELSE
             Tau = 0.0_dp
          END IF

          IF ( .NOT. Bubbles ) THEN
             DO p=1,n
                SU(p) = 0.0_dp
                DO i=1,dim-1
                   SU(p) = SU(p) + Vgauss(i) * dBasisdx(p,i)
                END DO

                SW(p) = 0.0_dp
                DO i=1,dim-1
                   SW(p) = SW(p) + Vgauss(i) * dBasisdx(p,i)
                END DO
             END DO
          END IF

          !        Stiffness matrix:
          !        -----------------
          DO p=1,n
             DO q=1,n
                DO i=1,DIM-1
                   STIFF(p,q) = STIFF(p,q) + &
                        s * Vgauss(i) * dBasisdx(q,i) * Basis(p)
                END DO
                STIFF(p,q) =  STIFF(p,q) + s * Tau * SU(q) * SW(p)
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

          !        Get accumulation/ablation function if flux input is given
          !        (i.e., calculate vector product between flux and normal)
          !        --------------------------------------------------------- 
          IF (.NOT.(SubstantialSurface)) THEN
             IF (NormalFlux) THEN 
                Source = normGradFreeSurf * SUM( SourceFunc(1:nCoord) &
                     * Basis(1:nCoord) )
             ELSE
                DO i=1,dim
                   FluxGauss(i) = SUM(Basis(1:nCoord)*Flux(i,1:nCoord))
                END DO
                Source = SUM(FluxGauss(1:DIM)*gradFreeSurf(1:DIM))
             END IF
          ELSE
             Source = 0.0_dp
          END IF

          !        Assemble force vector:
          !        ---------------------
          FORCE(1:n) = FORCE(1:n) &
               + (Vgauss(dim)+Source) * (Basis(1:n) + Tau*SW(1:n)) * s
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
   END SUBROUTINE FreeSurfaceSolver
!------------------------------------------------------------------------------
