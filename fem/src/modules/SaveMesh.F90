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
! *  A simple solver wrapper around WriteMeshToDisk2 to write the mesh
! *  to disk.
! *
! ******************************************************************************
! *
! *  Authors: Joe Todd, Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 15.09.2014
! *
! *****************************************************************************/

SUBROUTINE SaveMesh( Model,Solver,dt,TransientSimulation )

  USE DefUtils
  USE MeshUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!---------------------------------------------------------------
  LOGICAL :: SaveMeshLogical, EveryTime, Found
  INTEGER :: i, parts, ierr
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Variable_t), POINTER :: TimestepVar
  CHARACTER(LEN=MAX_NAME_LEN):: MeshName, MeshDir, inty, tmp

  Mesh => Solver % Mesh
  SaveMeshLogical = ListGetLogical(Solver % Values, "Save Mesh", Found)
  IF(.NOT.Found) THEN
    Call Warn("Save Mesh","Can't find Save Mesh logical")
    SaveMeshLogical = .FALSE.
  END IF

  IF(.NOT. SaveMeshLogical) RETURN !nothing to do

  TimestepVar => VariableGet( Mesh % Variables, "Timestep", .TRUE. )

  MeshName = ListgetString( Solver % Values,"Mesh Name", Found)
  IF(.NOT. Found) CALL FATAL("SaveMesh","No name given for mesh")

  MeshDir = ListgetString( Solver % Values,&
       "Save Mesh Directory", Found)
  IF(.NOT. Found) CALL FATAL("SaveMesh",&
       "No directory given to save mesh")

  EveryTime = ListGetLogical( Solver % Values, "Save All Timesteps", Found)
  IF(.NOT. Found) EveryTime = .FALSE.

  IF(EveryTime) THEN
     WRITE(MeshDir, '(A,A,A,i4.4)') TRIM(MeshDir),TRIM(MeshName),"_",INT(TimestepVar % Values(1))
  ELSE
     WRITE(MeshDir, '(A,A)')TRIM(MeshDir)//TRIM(MeshName)
  END IF

  IF( ParEnv % PEs<=1 ) THEN !serial
     CALL SYSTEM("mkdir -p "//MeshDir)
     CALL WriteMeshToDisk2(Model, Mesh, MeshDir)
  ELSE !parallel
     
     parts = ParEnv % PEs

     tmp = TRIM(MeshDir)//"/partitioning."
     MeshDir = TRIM(tmp)

     WRITE (tmp, '(A,I0)') TRIM(MeshDir),parts
     MeshDir = TRIM(tmp)

     IF(ParEnv % MyPe==0) THEN
        PRINT *, 'Save Mesh, creating directory...' !TEST
        CALL SYSTEM("mkdir -p "//MeshDir)        
     END IF
     CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

     CALL WriteMeshToDisk2(Model, Mesh, MeshDir, ParEnv % MyPe)
  END IF

END SUBROUTINE SaveMesh
