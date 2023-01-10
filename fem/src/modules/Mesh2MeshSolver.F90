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
! *  Module for performing nonconforming interpolation from mesh to mesh.
! *
! *  Author:  Peter Råback
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: June 2018
! *
! *****************************************************************************/
 
!> \ingroup Solvers
!> \{
!------------------------------------------------------------------------------
SUBROUTINE Mesh2MeshSolver_init0( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Found
    
  Params => GetSolverParams()

  CALL ListAddNewLogical( Params,'No Matrix',.TRUE.)
  CALL ListAddNewLogical( Params,'Mesh Enforce Local Copy',.TRUE.)
  
END SUBROUTINE Mesh2MeshSolver_Init0


!------------------------------------------------------------------------------
!> Interpolates fields from one mesh toanother. 
!------------------------------------------------------------------------------
SUBROUTINE Mesh2MeshSolver( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE AdditiveSol_t
    LOGICAL :: Additive
    REAL(KIND=dp), ALLOCATABLE :: Values(:)
  END TYPE AdditiveSol_t 
  TYPE(AdditiveSol_t), ALLOCATABLE :: BaseSol(:)
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh, ThisMesh, TargetMesh
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Variable_t), POINTER :: Var, pVar, FromVars
  CHARACTER(LEN=MAX_NAME_LEN) :: Name, VarName, MaskName
  INTEGER :: i, n, NoVar
  LOGICAL :: Found, GotMaskName, AdditiveMode, DoAdd, DoIt
  REAL(KIND=dp), ALLOCATABLE :: TmpSol(:)
  
  INTERFACE
    SUBROUTINE InterpolateMeshToMesh( OldMesh, NewMesh, OldVariables, &
        NewVariables, UseQuadrantTree, Projector, MaskName, UnfoundNodes )
      USE Lists
      USE SParIterComm
      USE Interpolation
      USE CoordinateSystems
      USE MeshUtils, ONLY: ReleaseMesh
      TYPE(Mesh_t), TARGET  :: OldMesh, NewMesh
      TYPE(Variable_t), POINTER, OPTIONAL :: OldVariables, NewVariables
      LOGICAL, OPTIONAL :: UseQuadrantTree
      TYPE(Projector_t), POINTER, OPTIONAL :: Projector
      CHARACTER(LEN=*),OPTIONAL :: MaskName
      LOGICAL, POINTER, OPTIONAL :: UnfoundNodes(:)
    END SUBROUTINE InterpolateMeshToMesh
  END INTERFACE

  
  CALL Info('Mesh2MeshSolver','Mapping result between meshes')
    
  ThisMesh => Getmesh()
  CALL Info('Mesh2MeshSolver','This mesh name is: '//TRIM(ThisMesh % Name),Level=7)   

  Params => GetSolverParams()

  TargetMesh => NULL()
  
  i = ListGetInteger( Params,'Target Mesh Solver Index',Found ) 
  IF( Found ) THEN
    ! Target mesh solver is explicitly given
    TargetMesh => CurrentModel % Solvers(i) % Mesh
    IF( ASSOCIATED( TargetMesh ) ) THEN
      CALL Info('Mesh2MeshSolver','Using target mesh as the mesh of Solver '//TRIM(I2S(i)),Level=8)
    ELSE
      CALL Fatal('Mesh2MeshSolver','Target Mesh for Solver not associated: '//TRIM(I2S(i)))
    END IF
  ELSE  
    ! Otherwise use the 1st mesh that is not this old data mesh
    Mesh => CurrentModel % Meshes      
    DO WHILE( ASSOCIATED(Mesh) ) 
      IF( .NOT. ASSOCIATED( Mesh, ThisMesh ) ) THEN
        TargetMesh => Mesh
        EXIT
      END IF
    END DO
    IF( ASSOCIATED( TargetMesh ) ) THEN
      CALL Info('Mesh2MeshSolver','Using target mesh as the first mesh different from this mesh',Level=8)
    ELSE
      CALL Fatal('Mesh2MeshSolver','Could not find target mesh that is not this mesh!')
    END IF
  END IF


  IF( ListGetLogical( Params,'Boundary Smoother', Found ) ) THEN
    CALL BoundarySmoother( Model,Solver,dt,Transient)
  END IF

  
  CALL Info('Mesh2MeshSolver','Target mesh name is: '//TRIM(TargetMesh % Name),Level=7)   

  AdditiveMode = .FALSE.
  NoVar = 0
  DO i = 1,100    
    WRITE (Name,'(A,I0)') 'Variable ',i
    IF( .NOT. ListCheckPresent( Params, Name ) ) EXIT
    NoVar = i
  END DO
   
  IF( NoVar > 0 ) THEN      
    CALL Info('Mesh2MeshSolver','Mapping '//TRIM(I2S(NoVar))//' fields as requested',Level=7)
        
    ALLOCATE( FromVars )
    pVar => FromVars

    AdditiveMode = ListCheckPrefix( Params,'Variable Additive')
    IF( AdditiveMode ) THEN
      ALLOCATE( BaseSol( NoVar ) )
      BaseSol(1:NoVar) % Additive = .FALSE.
      CALL Info('Mesh2MeshSolver','Interpolating in additive mode',Level=10)
    END IF
    
    DO i = 1,NoVar    
      WRITE (Name,'(A,I0)') 'Variable ',i
      VarName = GetString( Params, Name, Found )
      IF(.NOT. Found ) EXIT    

      IF( i > 1 ) THEN
        ALLOCATE( pVar % Next )
        pVar => pVar % Next
      END IF
            
      Var => VariableGet( ThisMesh % Variables, VarName, ThisOnly = .TRUE. )
      IF( ASSOCIATED( Var ) ) THEN
        pVar = Var
      ELSE
        pVar % Name = TRIM(VarName)
        pVar % dofs =  1
        ALLOCATE( PVar % Perm(0) )
      END IF
      pVar % Next => NULL()
      
      ! If the target variable has different name
      ! rename the primary variable also so that default subroutines can be used.
      !-------------------------------------------------------------------------
      WRITE (Name,'(A,I0)') 'Target Variable ',i
      VarName = GetString( Params, Name, Found )
      IF( Found ) THEN
        pVar % Name = VarName
        pVar % NameLen = LENTRIM( VarName ) 
      END IF

      ! If the variable is additive allocate and save for the base values
      !-------------------------------------------------------------------       
      WRITE (Name,'(A,I0)') 'Variable Additive ',i
      IF( ListGetLogical( Params,Name, Found ) ) THEN
        Var => VariableGet( TargetMesh % Variables, VarName, ThisOnly = .TRUE. )
        IF(.NOT. ASSOCIATED( Var ) ) THEN
          CALL Fatal('Mesh2MeshSolver','Only existing fields may be additive!')
        END IF
        BaseSol(i) % Additive = .TRUE.
        ALLOCATE( BaseSol(i) % Values( SIZE( Var % Values ) ) )
        BaseSol(i) % Values = Var % Values
      END IF
    END DO
  ELSE
    FromVars => ThisMesh % Variables 
  END IF
 
  ! For some reason this is not always active.
  ! If this is not set parallel interpolation could fail.
  !----------------------------------------------------------------------
  IF( ParEnv % PEs > 1 ) ParEnv % Active = .TRUE.

  ! Here is all the real interpolation work (serial & parallel)
  !---------------------------------------------------------------------
  CALL Info('Mesh2MeshSolver','Mapping all fields in primary mesh to target mesh!',Level=7)
  CALL InterpolateMeshToMesh( ThisMesh, TargetMesh, &
      FromVars, TargetMesh % Variables )
  CALL Info('Mesh2MeshSolver','Done mapping all fields in primary mesh to target mesh!',Level=15)
  !---------------------------------------------------------------------  

  ! If we has some additive variables then add the base values.
  !----------------------------------------------------------------
  IF( AdditiveMode ) THEN
    DO i = 1,NoVar
      IF( BaseSol(i) % Additive ) THEN
        WRITE (Name,'(A,I0)') 'Variable ',i
        VarName = GetString( Params, Name, Found )
        Var => VariableGet( TargetMesh % Variables, VarName, ThisOnly = .TRUE. )
        CALL Info('Mesh2MeshSolver','Adding variable to base field: '//TRIM(VarName),Level=7)
        Var % Values = Var % Values + BaseSol(i) % Values
        DEALLOCATE( BaseSol(i) % Values ) 
      END IF
    END DO
    DEALLOCATE( BaseSol ) 
  END IF
  
  ! Validate the mapped variables because otherwise we might accidentally be performing
  ! mapping again when calling GetVariable.
  !------------------------------------------------------------------------------------

  IF( NoVar > 0 ) THEN
    DO i = 1,NoVar    
      WRITE (Name,'(A,I0)') 'Variable ',i
      VarName = GetString( Params, Name, Found )
            
      Var => VariableGet( ThisMesh % Variables, VarName, ThisOnly = .TRUE. )
      pVar => VariableGet( TargetMesh % Variables, VarName, ThisOnly = .TRUE. )

      IF( ASSOCIATED( Var ) ) THEN
        IF( InfoActive(20) ) &
          PRINT *,'FromRange'//TRIM(I2S(ParEnv % Mype))//': ', TRIM(Var % Name),MAXVAL(Var % Values )      
      END IF
      IF( ASSOCIATED( pVar ) ) THEN
        pVar % Valid = .TRUE.
        pVar % ValuesChanged = .TRUE.
        IF( InfoActive(20) ) &
            PRINT *,'ToRange'//TRIM(I2S(ParEnv % Mype))//': ', TRIM(pVar % Name),MAXVAL(pVar % Values )      
      END IF
   END DO
  ELSE  
    Var => FromVars
    DO WHILE( ASSOCIATED( Var ) )
      ! Check that it exists as target variable
      pVar => VariableGet( TargetMesh % Variables, Var % Name(1:Var % NameLen), ThisOnly=.TRUE.)
      DoIt = ASSOCIATED( pVar )
      
      IF( DoIt ) THEN
        DoIt = ( SIZE( pVar % Values ) > Var % DOFs )  ! not global variable
        IF( DoIt ) DoIt = ( pVar % Name(1:10) /= 'coordinate' )  ! not coordinate
      END IF
      
      IF( DoIt ) THEN
        CALL Info('Mesh2MeshSolver','Validatig variable in target mesh: '//TRIM(Var % Name ),Level=7)
        pVar % Valid = .TRUE.
        pVar % ValuesChanged = .TRUE.
        
        IF( InfoActive(20) ) THEN
          PRINT *,'FromRange'//TRIM(I2S(ParEnv % Mype))//': ', TRIM(Var % Name),MAXVAL(Var % Values )      
          PRINT *,'ToRange'//TRIM(I2S(ParEnv % Mype))//': ', TRIM(pVar % Name),MAXVAL(pVar % Values )      
        END IF
      END IF
      
      Var => Var % Next
    END DO
  END IF
    
  CALL Info('Mesh2MeshSolver','Interpolated variables between meshes',Level=7)

CONTAINS


  !------------------------------------------------------------------------------
  !> A smoother that only uses topologilcal information of the mesh to detect the
  !> nodes that are not well joined.
  !------------------------------------------------------------------------------
  SUBROUTINE BoundarySmoother( Model,Solver,dt,Transient)
    !------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    TYPE(Model_t) :: Model
    REAL(KIND=dp) :: dt
    LOGICAL :: Transient
    !------------------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------------------
    TYPE(Element_t),POINTER :: Element
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Variable_t), POINTER :: Var
    INTEGER :: n,m,i,j,k,ii,jj,t,NoVar,BoundaryNodes,MaxDist,Loop,OperInd
    LOGICAL :: Visited=.FALSE., Found, Ready, Parallel
    INTEGER, POINTER :: BoundaryPerm(:), TopoDist(:), MinFriend(:), MaxFriend(:)
    INTEGER, POINTER :: TopoA(:), TopoB(:) 
    INTEGER, POINTER :: Indexes(:)
    REAL(KIND=dp) :: dx
    CHARACTER(LEN=MAX_NAME_LEN) :: VarName, OperName, Name
    TYPE (ParallelInfo_t), POINTER :: ParallelInfo
    CHARACTER(*), PARAMETER :: Caller = 'BoundarySmoother'
    !------------------------------------------------------------------------------

    CALL Info(Caller,'-------------------------------------------',Level=5)
    CALL Info(Caller,'Smoothing boundary',Level=5)
    CALL Info(Caller,'-------------------------------------------',Level=5)

    Mesh => GetMesh()
    Params => GetSolverParams()

    !For the moment the parallel communication is not activated
#define PARCOMP 0

    ParallelInfo => Mesh % ParallelInfo
    Parallel = ( ParEnv % PEs > 1) .AND. (.NOT. Mesh % SingleMesh)

    Visited = .TRUE.

    m = Solver % Mesh % NumberOfNodes
    ALLOCATE( BoundaryPerm(m) )
    BoundaryPerm = 0
    BoundaryNodes = 0
    CALL MakePermUsingMask( Model,Solver,Mesh,'smooth boundary', &
        .FALSE., BoundaryPerm, BoundaryNodes, ParallelComm = Parallel ) 

    ALLOCATE( TopoDist(m), MaxFriend(m), TopoA(m), TopoB(m) )

    MaxFriend = 0

    TopoDist = m
    WHERE(BoundaryPerm > 0)
      TopoDist = 0
    END WHERE

    j = ( COUNT(TopoDist == 0) )

    Ready = .FALSE.
    Loop = 0
    DO Loop = 1, m
      Ready = .TRUE.

      DO t = 1, Mesh % NumberOfBulkElements
        Element => Mesh % Elements(t)
        n = Element % Type % NumberOfNodes
        Indexes => Element % NodeIndexes
        DO i=1,n
          k = Indexes(i) 
          DO j=1,n
            IF(i==j) CYCLE
            IF(TopoDist(k) > TopoDist(Indexes(j))+1) THEN
              TopoDist(k) = TopoDist(Indexes(j))+1
              Ready = .FALSE.
            END IF
          END DO
        END DO
      END DO

#if PARCOMP
      IF( Parallel ) THEN
        CALL CommunicateParallelSystemTag(ParallelInfo,Itag = TopoDist,ParOper=1)
        i = 0
        IF( Ready ) i = 1
        i = ParallelReduction(i,1)
        IF(i==0) Ready = .FALSE.
      END IF
#endif
      
      IF( Ready ) EXIT
    END DO
    CALL Info(Caller,'Mesh coloring loops: '//TRIM(I2S(Loop)),Level=6)


    MaxDist = MAXVAL( TopoDist ) 
    CALL Info(Caller,'Maximum topology distance: '//TRIM(I2S(MaxDist)))

    DO i=0,MIN(10,MaxDist)
      j = COUNT(TopoDist==i)
      IF(j==0) EXIT
      CALL Info(Caller,'Node count with topology distance '//TRIM(I2S(i))//': '//TRIM(I2S(j)))
    END DO

    MaxFriend = 0
    DO t = 1, Mesh % NumberOfBulkElements
      Element => Mesh % Elements(t)        
      Indexes => Element % NodeIndexes
      n = Element % TYPE % NumberOfNodes
      DO i=1,n
        ii = Indexes(i)
        DO j=1,n
          IF(i==j) CYCLE
          jj = Indexes(j)
          MaxFriend(ii) = MAX(MaxFriend(ii),TopoDist(jj))
        END DO
      END DO
    END DO

#if PARCOMP
    IF( Parallel ) THEN
      CALL CommunicateParallelSystemTag(ParallelInfo,Itag = MaxFriend,ParOper=2)
    END IF
#endif
    
    DO i=0,5
      j = COUNT(MaxFriend==i)
      CALL Info(Caller,'Node count with max friend '//TRIM(I2S(i))//': '//TRIM(I2S(j)))
    END DO

    DO i=0,4
      DO j=0,4
        k = COUNT(TopoDist==i .AND. MaxFriend == j )
        IF(k>0) CALL Info(Caller,'Topo dist '//TRIM(I2S(i))//' with max friend '//TRIM(I2S(j))//': '//TRIM(I2S(k)))
      END DO
    END DO

    TopoA = 0
    TopoB = 0
    DO t = 1, Mesh % NumberOfBulkElements
      Element => Mesh % Elements(t)        
      Indexes => Element % NodeIndexes
      n = Element % TYPE % NumberOfNodes
      DO i=1,n
        ii = Indexes(i)
        DO j=1,n
          IF(i==j) CYCLE
          jj = Indexes(j)
          TopoA(ii) = TopoA(ii) + TopoDist(jj)
          TopoB(ii) = TopoB(ii) + TopoDist(jj)*MaxFriend(jj)**3
        END DO
      END DO
    END DO

    DO i=0,10
      j = COUNT(TopoA==i)
      IF(j==0) CYCLE
      CALL Info(Caller,'Topo measure A '//TRIM(I2S(i))//': '//TRIM(I2S(j)))
    END DO

    DO i=0,10
      j = COUNT(TopoB==i)
      IF(j==0) CYCLE
      CALL Info(Caller,'Topo measure B '//TRIM(I2S(i))//': '//TRIM(I2S(j)))
    END DO

    ! If we have "exported variables" for the fields then we can visualize them.
    Var => VariableGet( Mesh % Variables,'Topo Dist')
    IF( ASSOCIATED(Var) ) Var % Values = 1.0_dp * TopoDist

    Var => VariableGet( Mesh % Variables,'Max Topo Dist')
    IF( ASSOCIATED(Var) ) Var % Values = 1.0_dp * MaxFriend 

    Var => VariableGet( Mesh % Variables,'Topo A')
    IF( ASSOCIATED(Var) ) Var % Values = 1.0_dp * TopoA

    j = MINVAL(TopoB)
    Var => VariableGet( Mesh % Variables,'Topo B')
    IF( ASSOCIATED(Var) ) Var % Values = 1.0_dp * TopoB 

    NoVar = 0
    DO WHILE(.TRUE.)    
      NoVar = NoVar + 1

      WRITE (Name,'(A,I0)') 'Smooth Variable ',NoVar  
      VarName = ListGetString( Params, TRIM(Name), Found )
      IF(.NOT. Found ) EXIT

      CALL Info(Caller,'Smoothing variable '//TRIM(VarName))

      Var => VariableGet( Mesh % Variables, VarName )
      IF(.NOT. ASSOCIATED(Var) ) THEN
        CALL Fatal(Caller,'Could not access variable!')
      END IF
      
      ! Currently only one simple smoother     
      !WRITE (Name,'(A,I0)') 'Operator ',NoVar
      !OperName = ListGetString(Params,TRIM(Name),UnfoundFatal=.TRUE.)
      !CALL Info(Caller,'Smoothing variable "'//TRIM(VarName)//'" with operator "'//TRIM(OperName)//'"',Level=5)
      !OperInd = ListGetInteger(Params,TRIM(Name),UnfoundFatal=.TRUE.)
      !CALL Info(Caller,'Smoothing variable "'//TRIM(VarName)//'" with operator: '//I2S(OperInd),Level=5)

      DO t=1,Mesh % NumberOfBulkElements
        Element => Mesh % Elements(t)
        Indexes => Element % NodeIndexes
        n = Element % TYPE % NumberOfNodes
        DO i = 1,n        
          ii = Indexes(i)
          ! If TopoB==2 then we have a node that is only attached to one inland node that
          ! has no neighbours deeper inland. 
          IF(TopoB(ii) > 2) CYCLE 

          DO j=1,n
            jj = Indexes(j)
            IF(TopoDist(jj) == 0) CYCLE
            ! We are only taking the value, so no need to average or book keep etc. 
            IF( ASSOCIATED( Var % Perm ) ) THEN
              dx = Var % Values(Var % Perm(jj)) - Var % Values(Var % Perm(ii))
              IF( ABS(dx) > EPSILON(dx) ) THEN
                !PRINT *,'copy:',ii,jj,Var % Values(Var % Perm(ii)), dx
                Var % Values(Var % Perm(ii)) = Var % Values(Var % Perm(jj))
              END IF
            ELSE
              Var % Values(ii) = Var % Values(jj)
            END IF
          END DO
        END DO
      END DO
    END DO

    CALL Info(Caller,'Smoothing done!')
    
    DEALLOCATE(BoundaryPerm, TopoDist, MaxFriend, TopoA, TopoB ) 
    
    !------------------------------------------------------------------------------
  END SUBROUTINE BoundarySmoother
  !------------------------------------------------------------------------------
  
!------------------------------------------------------------------------------
END SUBROUTINE Mesh2MeshSolver
!------------------------------------------------------------------------------

!> \}
