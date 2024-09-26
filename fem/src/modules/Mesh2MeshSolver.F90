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
SUBROUTINE Mesh2MeshSolver_init0( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
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
SUBROUTINE Mesh2MeshSolver( Model,Solver,dt,TransientSimulation )
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
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh, ThisMesh, TargetMesh
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Variable_t), POINTER :: Var, pVar, FromVars
  CHARACTER(LEN=MAX_NAME_LEN) :: Name, VarName, MaskName
  INTEGER :: i, j, n, NoVar
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
      CALL Info('Mesh2MeshSolver','Using target mesh as the mesh of Solver '//I2S(i),Level=8)
    ELSE
      CALL Fatal('Mesh2MeshSolver','Target Mesh for Solver not associated: '//I2S(i))
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

  CALL Info('Mesh2MeshSolver','Target mesh name is: '//TRIM(TargetMesh % Name),Level=7)   

  AdditiveMode = .FALSE.
  NoVar = 0
  DO i = 1,100    
    WRITE (Name,'(A,I0)') 'Variable ',i
    IF( .NOT. ListCheckPresent( Params, Name ) ) EXIT
    NoVar = i
  END DO
   
  IF( NoVar > 0 ) THEN      
    CALL Info('Mesh2MeshSolver','Mapping '//I2S(NoVar)//' hand-picked fields as requested',Level=7)
        
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

      ! Add a new variable to the temporal list.
      IF( i > 1 ) THEN
        ALLOCATE( pVar % Next )
        pVar => pVar % Next
      END IF
            
      Var => VariableGet( ThisMesh % Variables, VarName, ThisOnly = .TRUE. )
      IF( ASSOCIATED( Var ) ) THEN
        pVar = Var
#if 0 
        IF(.NOT. ASSOCIATED(pVar % Perm) ) THEN
          ALLOCATE(pVar % Perm(ThisMesh % NumberOfNodes))
          DO j=1,ThisMesh % NumberOfNodes
            pVar % Perm(j) = j
          END DO

          PRINT *,'copied:',SIZE(Var % Values), TRIM(Var % Name), ASSOCIATED(Var % Perm), Var % Dofs
        END IF
#endif
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
        CALL Info('Mesh2MeshSolver','Renaming variable "'//TRIM(pVar % Name)//'" to "'//TRIM(VarName)//'"')
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
          PRINT *,'FromRange'//I2S(ParEnv % Mype)//': ', TRIM(Var % Name),MAXVAL(Var % Values )      
      END IF
      IF( ASSOCIATED( pVar ) ) THEN
        pVar % Valid = .TRUE.
        pVar % ValuesChanged = .TRUE.
        IF( InfoActive(20) ) &
            PRINT *,'ToRange'//I2S(ParEnv % Mype)//': ', TRIM(pVar % Name),MAXVAL(pVar % Values )      
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
        IF ( Doit .AND. LEN(pVar % Name)>=10  ) THEN
          DoIt = ( pVar % Name(1:10) /= 'coordinate' )  ! not coordinate
        ELSE
          DoIt = .FALSE.
        END IF
      END IF
      
      IF( DoIt ) THEN
        CALL Info('Mesh2MeshSolver','Validatig variable in target mesh: '//TRIM(Var % Name ),Level=7)
        pVar % Valid = .TRUE.
        pVar % ValuesChanged = .TRUE.
        
        IF( InfoActive(20) ) THEN
          PRINT *,'FromRange'//I2S(ParEnv % Mype)//': ', TRIM(Var % Name),MAXVAL(Var % Values )      
          PRINT *,'ToRange'//I2S(ParEnv % Mype)//': ', TRIM(pVar % Name),MAXVAL(pVar % Values )      
        END IF
      END IF
      
      Var => Var % Next
    END DO
  END IF
    
  CALL Info('Mesh2MeshSolver','Interpolated variables between meshes',Level=7)
  
!------------------------------------------------------------------------------
END SUBROUTINE Mesh2MeshSolver
!------------------------------------------------------------------------------

!> \}
