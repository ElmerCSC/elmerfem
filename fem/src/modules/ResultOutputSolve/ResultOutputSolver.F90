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


!------------------------------------------------------------------------------
!>  Exports data to other suitable postprocessing software.
!>  Currently supported formats are ElmerPost, GiD, VTK legacy, VTK XML, and Open DX.
!> This is a dynamically loaded solver with a standard interface.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE ResultOutputSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------

  USE DefUtils
  USE SaveUtils
  USE AscBinOutputUtils
  
  IMPLICIT NONE
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation

  LOGICAL :: SaveGid, SaveVTK, SaveOpenDx, SaveGmsh, &
      SaveVTU, SaveEP, SaveAny, ListSet = .FALSE., ActiveMesh, &
      SomeMeshSaved, SaveAllMeshes
  INTEGER :: i,nInterval=1, nstep=0, OutputCount(6) = 0, MeshDim,&
      MinMeshDim,MaxMeshDim,MeshLevel,nlen,NoMeshes, m
  INTEGER, POINTER :: OutputIntervals(:), TimeSteps(:)

  TYPE(Mesh_t), POINTER :: Mesh, iMesh, MyMesh
  CHARACTER(10) :: OutputFormat
  CHARACTER(LEN=MAX_NAME_LEN) :: FilePrefix, MeshName, iMeshName, ListMeshName
  LOGICAL :: SubroutineVisited=.FALSE.,Found, SaveThisMesh, NowSave
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Variable_t), POINTER :: ModelVariables
  CHARACTER(*), PARAMETER :: Caller = 'ResultOutputSolver'
  INTEGER :: SaveSolverMeshIndex
  LOGICAL :: CalcNrm
  REAL(KIND=dp) :: Nrm
  REAL(KIND=dp), POINTER :: RefResults(:,:)
  
  
  SAVE SubroutineVisited, OutputCount, ListSet, MeshDim, ListMeshName

  INTERFACE
    RECURSIVE SUBROUTINE ElmerPostOutputSolver( Model, Solver,dt,TransientSimulation,ONOEfound )
      USE DefUtils
      IMPLICIT NONE
!------------------------------------------------------------------------------
      TYPE(Solver_t) :: Solver
      TYPE(Model_t) :: Model
      REAL(KIND=dp) :: dt
      LOGICAL, OPTIONAL :: ONOEfound
      LOGICAL :: TransientSimulation
    END SUBROUTINE ElmerPostOutputSolver
  END INTERFACE

  CALL Info( Caller, '-------------------------------------')

  Params => GetSolverParams()
  SaveGid = GetLogical(Params,'Gid Format',Found)
  SaveGmsh = GetLogical(Params,'Gmsh Format',Found)
  SaveVTK = GetLogical(Params,'VTK Format',Found)
  SaveVTU = GetLogical(Params,'VTU Format',Found)
  SaveOpenDx = GetLogical(Params,'Dx Format',Found)
  SaveEP = GetLogical(Params,'Elmerpost Format',Found)

  OutputFormat = GetString( Params, 'Output Format', Found )
  IF(Found) THEN
    IF( OutputFormat == "gid" )THEN
      SaveGid = .TRUE.
    ELSE IF( OutputFormat == "vtk" )THEN
      SaveVTK = .TRUE.
    ELSE IF( OutputFormat == "vtu" )THEN
      SaveVTU = .TRUE.
    ELSE IF( OutputFormat == "dx" )THEN
      SaveOpenDx = .TRUE.
    ELSE IF( OutputFormat == "gmsh" )THEN
      SaveGmsh = .TRUE.
    ELSE IF( OutputFormat == "elmerpost" )THEN
      SaveEP = .TRUE.
    ELSE
      CALL Warn( Caller, &
                 'Unknown output format "' // TRIM(OutputFormat) // '"' )
      CALL Warn( Caller, &
                 'Available formats are "GiD", "VTK", "VTU", "DX", "gmsh" and "elmerpost"' )
      RETURN
    END IF
  END IF

  MyMesh => GetMesh()

  IF( .NOT. SubroutineVisited ) THEN
    IF ( GetLogical(Params,'Show Variables',Found) ) THEN
      CALL CreateListForSaving( Model, Params,.TRUE. )    
    END IF
  END IF

  SaveAny = SaveGid .OR. SaveVTK .OR. SaveVTU .OR. SaveOpenDX .OR. SaveGmsh .OR. SaveEp
  IF(.NOT. SaveAny ) THEN
    CALL Warn(Caller,'No output format given, assuming VTU')
    SaveVTU = .TRUE.
  END IF

  FilePrefix = GetString( Params,'Output File Name',Found)
  IF(.NOT. Found) THEN
    FilePrefix = 'Case'
    CALL ListAddString( Params,'Output File Name',FilePrefix)
  END IF
  
  IF( .NOT. SubroutineVisited ) THEN 
    CALL Info(Caller,'Saving with prefix: '//TRIM(FilePrefix))
  END IF	

  ! The idea of this is that the independent subroutines may be called 
  ! with different data sets and still maintaining the standard output calling convention
  IF(SaveVtu)  OutputCount(1) = OutputCount(1) + 1
  IF(SaveGmsh) OutputCount(2) = OutputCount(2) + 1
  IF(SaveVTK)  OutputCount(3) = OutputCount(3) + 1
  IF(SaveEp)   OutputCount(4) = OutputCount(4) + 1
  IF(SaveGid)  OutputCount(5) = OutputCount(5) + 1
  IF(SaveOpenDx) OutputCount(6) = OutputCount(6) + 1
  
  ! Finally go for it and write desired data
  ! Some formats requite that the list of variables is explicitly given
  !-----------------------------------------

  MeshLevel = GetInteger( Params,'Output Mesh Level',Found)
  SomeMeshSaved = .FALSE.

  SaveAllMeshes = GetLogical( Params,'Save All Meshes',Found ) 
  SaveThisMesh = GetLogical( Params,'Save This Mesh Only',Found ) 
  SaveSolverMeshIndex = GetInteger( Params,'Save Solver Mesh Index',Found ) 
  
  MinMeshDim = ListGetInteger( Params,'Minimum Mesh Dimension',Found )
  MaxMeshDim = ListGetInteger( Params,'Maximum Mesh Dimension',Found )

  RefResults => ListGetConstRealArray( Params,'Reference Values',CalcNrm )
  CALL AscBinInitNorm(CalcNrm) 

  ! Loop over the meshes and save them using the selected format(s).
  ! First iteration just count the meshes. 
  !----------------------------------------------------------------------------------  
  NowSave = .FALSE.
1 NoMeshes = 0
  m = 1
  iMesh => Model % Meshes
  DO WHILE( ASSOCIATED(iMesh) )

    IF ( .NOT. SaveAllMeshes .AND. .NOT. iMesh % OutputActive ) THEN
      IF(NowSave) CALL Info(Caller,'Skipping inactive mesh: '//TRIM(iMesh % Name), Level=10 )
      iMesh => iMesh % next; m=m+1
      CYCLE 
    END IF    

    IF( SaveThisMesh ) THEN
      IF( .NOT. ASSOCIATED( iMesh, MyMesh ) ) THEN
        IF(NowSave) CALL Info(Caller,'Skipping not my mesh: '//TRIM(iMesh % Name), Level=10 )
        iMesh => iMesh % next; m=m+1
        CYCLE
      END IF
    END IF

    IF( SaveSolverMeshIndex > 0 ) THEN
      IF( .NOT. ASSOCIATED( iMesh, Model % Solvers(SaveSolverMeshIndex) % Mesh ) ) THEN
        iMesh => iMesh % next; m=m+1
        CYCLE
      END IF
    END IF
        
    IF(NowSave) CALL Info(Caller,'Dimension of mesh is: '//I2S(iMesh % MeshDim),Level=10)

    IF( MinMeshDim /= 0 .AND. iMesh % MeshDim < MinMeshDim ) THEN
      IF(NowSave) CALL Info(Caller,'Skipping lower dimensional mesh: '//TRIM(iMesh % Name), Level=10 )
      iMesh => iMesh % next; m=m+1
      CYCLE
    END IF

    IF( MaxMeshDim /= 0 .AND. iMesh % MeshDim > MaxMeshDim ) THEN
      IF(NowSave) CALL Info(Caller,'Skipping higher dimensional mesh: '//TRIM(iMesh % Name), Level=10 )
      iMesh => iMesh % next; m=m+1
      CYCLE
    END IF

    ! Optionally skip the writing of given meshes 
    !---------------------------------------------------------------    
    nlen = StringToLowerCase( iMeshName, iMesh % Name ) 
    MeshName = GetString( Params,'Mesh Name',Found )
    IF(Found) THEN      
      i = StringToLowerCase( MeshName, MeshName )      
      Found = ( i <= nlen ) 
      IF( Found ) Found = ( MeshName(1:i) == iMeshName(1:i) )
      
      IF( .NOT. Found ) THEN
        IF(NowSave) CALL Info(Caller,'Skipping mesh with mismatching name: '//TRIM(iMesh % Name), Level=10 )
        iMesh => iMesh % next; m=m+1
        CYCLE 
      END IF
    END IF
    
    IF(.NOT. NowSave ) THEN
      ! Discont mesh will get a separate prefix anyways so don't count that as
      ! a mesh competing from the same directory.     
      IF ( .NOT. iMesh % DiscontMesh ) THEN    
        NoMeshes = NoMeshes + 1
      END IF
      iMesh => iMesh % Next; m=m+1
      CYCLE
    END IF

    IF(NowSave) THEN
      CALL Info(Caller,'Working on mesh '//I2S(m)//': '//TRIM(iMesh % Name)//&
          ' with '//I2S(iMesh % NumberOfNodes)//' nodes', Level=10)
    END IF
      
    CALL SetCurrentMesh( Model, iMesh )
    ModelVariables => Model % Variables
    Model % Variables => iMesh % variables 

    IF( .NOT. ListSet ) THEN
      CALL Info(Caller,'Creating list for saving - if not present')
      CALL CreateListForSaving( Model, Params,.TRUE. )    
      ListSet = .TRUE.
    ELSE IF( MeshDim /= Model % Mesh % MeshDim .OR. (iMeshName(1:nlen) /= ListMeshName(1:nlen))) THEN
      CALL Info(Caller,'Recreating list for saving')
      CALL CreateListForSaving( Model, Params,.TRUE.,.TRUE.)
    END IF

    MeshDim = Model % Mesh % MeshDim
    nlen = StringToLowerCase( ListMeshName, iMesh % Name)

    Mesh => iMesh

    ! In case there are multiple mesh levels one may also save coarser ones
    !----------------------------------------------------------------------
    DO i=1,MeshLevel
      IF (.NOT.ASSOCIATED(Mesh % Parent)) EXIT
      Mesh => Mesh % Parent
    END DO

    IF ( ASSOCIATED(Mesh)) THEN
      CALL SetCurrentMesh( Model, Mesh )
      Model % Variables => Mesh % variables 
      SomeMeshSaved = .TRUE.

      IF( SaveVTU ) THEN
        CALL Info( Caller,'Saving in unstructured VTK XML (.vtu) format' )               
        CALL ListAddInteger( Params,'Output Count',OutputCount(1))
        CALL VtuOutputSolver( Model,Solver,dt,TransientSimulation )
      END IF
      IF( SaveGmsh ) THEN
        CALL Info( Caller,'Saving in gmsh 2.0 (.msh) format' )
        CALL ListAddInteger( Params,'Output Count',OutputCount(2))      
        ! For other call uses this recides in SaveUtils.
        CALL SaveGmshOutput( Model,Solver,dt,TransientSimulation )
      END IF
      IF( SaveVTK ) THEN
        CALL Info( Caller,'Saving in legacy VTK (.vtk) format' )
        CALL ListAddInteger( Params,'Output Count',OutputCount(3))
        CALL VtkOutputSolver( Model,Solver,dt,TransientSimulation )
      END IF
      IF( SaveEP ) THEN
        CALL Info( Caller,'Saving in ElmerPost (.ep) format' )
        CALL ListAddInteger( Params,'Output Count',OutputCount(4))
        CALL ElmerPostOutputSolver( Model,Solver,dt,TransientSimulation )
      END IF
      IF( SaveGid ) THEN
        CALL Info( Caller,'Saving in GiD format' )    
        CALL ListAddInteger( Params,'Output Count',OutputCount(5))
        CALL GiDOutputSolver( Model,Solver,dt,TransientSimulation )
      END IF
      IF( SaveOpenDx ) THEN
        CALL Info( Caller,'Saving in OpenDX (.dx) format' )
        CALL ListAddInteger( Params,'Output Count',OutputCount(6))
        CALL DXOutputSolver( Model,Solver,dt,TransientSimulation )
      END IF

      CALL Info( Caller, '-------------------------------------')
    END IF
    
    iMesh => iMesh % Next; m=m+1
  END DO

  IF( .NOT. NowSave ) THEN
    CALL ListAddInteger( Params,'Number of Output Meshes',NoMeshes)
    CALL Info(Caller,'Number of output meshes: '//I2S(NoMeshes),Level=12)
    NowSave = .TRUE.
    GOTO 1
  END IF
    
  IF( .NOT. SomeMeshSaved ) THEN
    OutputCount = OutputCount - 1
  END IF
  Model % Variables => ModelVariables

  SubroutineVisited = .TRUE.

  IF( CalcNrm ) THEN
    IF( SaveVtu ) THEN
      Nrm = AscBinCompareNorm(RefResults(:,1))
      Solver % Variable % Norm = Nrm
      WRITE( Message,'(A,ES15.6)' ) 'Calculate Pseudonorm:',Nrm
      CALL Info(Caller, Message)
    ELSE
      CALL Warn(Caller,'Reference norm computation implemented only for VTU format!')
    END IF
  END IF

  
END SUBROUTINE ResultOutputSolver
