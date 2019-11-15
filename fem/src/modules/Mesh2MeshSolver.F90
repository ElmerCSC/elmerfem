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
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh, ThisMesh, TargetMesh
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Variable_t), POINTER :: Var
  CHARACTER(LEN=MAX_NAME_LEN) :: Name, VarName, MaskName
  INTEGER :: i, n
  LOGICAL :: Found, GotMaskName, AdditiveMode, DoAdd  
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

  CALL Info('Mesh2MeshSolver','Target mesh name is: '//TRIM(TargetMesh % Name),Level=7)   

  
  IF( ListCheckPresent( Params,'Variable 1') ) THEN
    CALL Info('Mesh2MeshSolver','Mapping field one at a time as requested',Level=7)
    
    CALL SetCurrentMesh( CurrentModel, TargetMesh )

    AdditiveMode = ListGetLogical( Params,'Interpolation Additive',Found ) 
    IF( AdditiveMode ) THEN
      CALL Info('Mesh2MeshSolver','Interpolating in additive mode',Level=15)
    END IF

    DO i = 1,100    
      WRITE (Name,'(A,I0)') 'Variable ',i
      VarName = GetString( Params, Name, Found )
      IF(.NOT. Found ) EXIT    

      WRITE (Name,'(A,I0)') 'Mask ',i
      MaskName = GetString( Params, Name, GotMaskName )

      ! Use namespace such that in principle each variable could have different set of
      ! interpolation rules attached to them. 
      CALL ListPushNameSpace('var'//TRIM(I2S(i))//':')

      ! Here we might invalidate the variable in the primary mesh so that it really needs to be interpolated
      Var => VariableGet( TargetMesh % Variables, VarName, ThisOnly = .TRUE. )
      IF( ASSOCIATED( Var ) ) THEN
        Var % Valid = .FALSE.
      END IF

      ! This is provided if we want the interpolation to be cumulative
      DoAdd = .FALSE.
      IF( ASSOCIATED( Var ) ) THEN      
        IF( AdditiveMode ) THEN
          n = SIZE( Var % Values ) 
          IF(ALLOCATED( TmpSol ) ) THEN
            IF( SIZE( TmpSol ) < n ) DEALLOCATE( TmpSol ) 
          END IF
          IF( .NOT. ALLOCATED( TmpSol ) ) THEN
            ALLOCATE( TmpSol( n ) )
          END IF
          DoAdd = .TRUE.
          TmpSol(1:n) = Var % Values(1:n)
        END IF
      END IF

      ! Try to find the variable in target mesh, this includes MeshToMesh interpolation by default
      IF( GotMaskName ) THEN      
        Var => VariableGet( TargetMesh % Variables, VarName, MaskName = MaskName )
      ELSE 
        Var => VariableGet( TargetMesh % Variables, VarName )
      END IF

      IF( DoAdd ) THEN
        Var % Values(1:n) = Var % Values(1:n) + TmpSol(1:n)
      END IF

      IF(.NOT. ASSOCIATED( Var ) ) THEN
        CALL Warn('Mesh2MeshSolver','Could not find variable '//TRIM(VarName)//' in part '//TRIM(I2S(ParEnv % MyPe)))
      END IF

      PRINT *,'Variable range:',MINVAL( Var % Values), MAXVAL( Var % Values ), &
          SUM( Var % Values ) / SIZE( Var % Values ) 

      CALL ListPopNamespace()
    END DO
    CALL SetCurrentMesh( CurrentModel, ThisMesh )
    CALL Info('Mesh2MeshSolver','Succesfully interpolated '//TRIM(I2S(i-1))//' variables',Level=7)
  ELSE
    CALL Info('Mesh2MeshSolver','Mapping all fields in primary mesh to target mesh!',Level=7)
    CALL InterpolateMeshToMesh( ThisMesh, TargetMesh, &
        ThisMesh % Variables, TargetMesh % Variables )
    CALL Info('Mesh2MeshSolver','Interpolated variables between meshes',Level=7)
  END IF
  
!------------------------------------------------------------------------------
END SUBROUTINE Mesh2MeshSolver
!------------------------------------------------------------------------------

!> \}
