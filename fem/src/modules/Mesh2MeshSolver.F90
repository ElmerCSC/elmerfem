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

!  CALL ListAddNewLogical( Params,'No Matrix',.TRUE.)

  
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
  CHARACTER(LEN=MAX_NAME_LEN) :: Name, VarName
  INTEGER :: i
  LOGICAL :: Found
  
  CALL Info('Mesh2MeshSolver','Mapping result between meshes')

  
  ThisMesh => Getmesh()
  Params => GetSolverParams()

  TargetMesh => NULL()

  i = ListGetInteger( Params,'Target Mesh Solver Index',Found ) 
  IF( Found ) THEN
    ! Target mesh solver is explicitely given
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
  
  CALL SetCurrentMesh( CurrentModel, TargetMesh )
  DO i = 1,100    
    WRITE (Name,'(A,I0)') 'Variable ',i
    VarName = GetString( Params, Name, Found )
    IF(.NOT. Found ) EXIT
    
    ! Use namespace such that in principle each variable could have different set of
    ! interpolation rules attached to them. 
    CALL ListPushNameSpace('var'//TRIM(I2S(i))//':')

    ! Here we might invalidate the variable in the primary mesh so that it really needs to be interpolated
    Var => VariableGet( TargetMesh % Variables, VarName, ThisOnly = .TRUE. )
    IF( ASSOCIATED( Var ) ) Var % Valid = .FALSE.
    
    ! Try to find the variable in target mesh, this includes MeshToMesh interpolation by default
    Var => VariableGet( TargetMesh % Variables, VarName )
    IF(.NOT. ASSOCIATED( Var ) ) THEN
      CALL Warn('Mesh2MeshSolver','Could not find variable '//TRIM(VarName)//' in part '//TRIM(I2S(ParEnv % MyPe)))
    ELSE
      PRINT *,'Variable range: ',i,parenv % mype, MINVAL(Var % values), MAXVAL(Var % values)
    END IF
    
    CALL ListPopNamespace()
  END DO
  CALL SetCurrentMesh( CurrentModel, ThisMesh )
  
  CALL Info('Mesh2MeshSolver','Succesfully interpolated '//TRIM(I2S(i-1))//' variables',Level=7)

  
!------------------------------------------------------------------------------
END SUBROUTINE Mesh2MeshSolver
!------------------------------------------------------------------------------

!> \}
