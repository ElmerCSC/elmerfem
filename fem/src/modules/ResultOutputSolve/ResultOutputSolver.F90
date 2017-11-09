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

  IMPLICIT NONE
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation

  LOGICAL :: SaveGid, SaveVTK, SaveOpenDx, SaveGmsh, &
      SaveVTU, SaveEP, SaveAny, ListSet = .FALSE., ActiveMesh, &
      SomeMeshSaved, SaveAllMeshes
  INTEGER :: i,nInterval=1, nstep=0, OutputCount = 0, MeshDim,MeshLevel,nlen,NoMeshes
  INTEGER, POINTER :: OutputIntervals(:), TimeSteps(:)

  TYPE(Mesh_t), POINTER :: Mesh, iMesh, MyMesh
  CHARACTER(10) :: OutputFormat
  CHARACTER(LEN=MAX_NAME_LEN) :: FilePrefix, MeshName, iMeshName, ListMeshName
  LOGICAL :: SubroutineVisited=.FALSE.,Found, SaveThisMesh
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Variable_t), POINTER :: ModelVariables
    
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

  CALL Info( 'ResultOutputSolver', '-------------------------------------')

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
      CALL Warn( 'ResultOutputSolver', &
                 'Unknown output format "' // TRIM(OutputFormat) // '"' )
      CALL Warn( 'ResultOutputSolver', &
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
    CALL Warn('ResultOutputSolver','No output format given, assuming VTU')
    SaveVTU = .TRUE.
  END IF

  FilePrefix = GetString( Params,'Output File Name',Found)
  IF(.NOT. Found) THEN
    FilePrefix = 'Case'
    CALL ListAddString( Params,'Output File Name',FilePrefix)
  END IF
  
  IF( .NOT. SubroutineVisited ) THEN 
    CALL Info('ResultOutputSolve','Saving with prefix: '//TRIM(FilePrefix))
  END IF	


  ! The idea of this is that the independent subroutines may be called 
  ! with different data sets and still maintaining the standard output calling convention
  OutputCount = OutputCount + 1
  CALL ListAddInteger( Params,'Output Count',OutputCount)

  ! Finally go for it and write desired data
  ! Some formats requite that the list of variables is explicitely given
  !-----------------------------------------

  MeshLevel = GetInteger( Params,'Output Mesh Level',Found)
  SomeMeshSaved = .FALSE.

  SaveAllMeshes = GetLogical( Params,'Save All Meshes',Found ) 
  SaveThisMesh = GetLogical( Params,'Save This Mesh Only',Found ) 


  ! This is similar cycle as below but only to count the number of meshes to save
  NoMeshes = 0 
  iMesh => Model % Meshes
  DO WHILE( ASSOCIATED(iMesh) )
    IF ( .NOT. SaveAllMeshes .AND. .NOT. iMesh % OutputActive ) THEN
      iMesh => iMesh % next
      CYCLE 
    END IF    

    IF( SaveThisMesh ) THEN
      IF( .NOT. ASSOCIATED( iMesh, MyMesh ) ) THEN
        iMesh => iMesh % next
        CYCLE
      END IF
    END IF

    IF( iMesh % MeshDim < ListGetInteger( Params,'Minimum Mesh Dimension',Found )  ) THEN
      iMesh => iMesh % next
      CYCLE
    END IF

    nlen = StringToLowerCase( iMeshName, iMesh % Name ) 
    MeshName = GetString( Params,'Mesh Name',Found )
    IF(Found) THEN
      nlen = LEN_TRIM(MeshName)
      IF( MeshName(1:nlen) /= iMeshName(1:nlen) ) THEN
        iMesh => iMesh % next
        CYCLE
      END IF
    END IF

    ! Discont mesh will get a separate prefix anyways so don't count that as
    ! a mesh competing from the same directory. 
    IF ( .NOT. iMesh % DiscontMesh ) THEN    
      NoMeshes = NoMeshes + 1
    END IF
    
    iMesh => iMesh % next
  END DO
  CALL ListAddInteger( Params,'Number of Output Meshes',NoMeshes)
  CALL Info('ResultOutputSolve','Number of output meshes: '//TRIM(I2S(NoMeshes)),Level=12)
  

  ! Now this is the real thing. Now we cycle over the meshes and save them
  ! using the selected format(s).
  !----------------------------------------------------------------------------------  
  iMesh => Model % Meshes
  DO WHILE( ASSOCIATED(iMesh) )
    
    CALL Info('ResultOutputSolver','Working on mesh: '//TRIM(iMesh % Name), Level=7 )

    IF ( .NOT. SaveAllMeshes .AND. .NOT. iMesh % OutputActive ) THEN
      CALL Info('ResultOutputSolver','Skipping mesh: '//TRIM(iMesh % Name), Level=7 )
      iMesh => iMesh % next
      CYCLE 
    END IF    

    IF( SaveThisMesh ) THEN
      IF( .NOT. ASSOCIATED( iMesh, MyMesh ) ) THEN
        CALL Info('ResultOutputSolver','Skipping mesh: '//TRIM(iMesh % Name), Level=7 )
        iMesh => iMesh % next
        CYCLE
      END IF
    END IF

    CALL Info('ResultOutputSolver','Dimension of mesh is: '//TRIM(I2S(iMesh % MeshDim)),Level=7)
    IF( iMesh % MeshDim < ListGetInteger( Params,'Minimum Mesh Dimension',Found )  ) THEN
      CALL Info('ResultOutputSolver','Skipping meshes with too low dimension')
      iMesh => iMesh % next
      CYCLE
    END IF


    ! Optionally skip the writing of given meshes 
    !---------------------------------------------------------------    
    nlen = StringToLowerCase( iMeshName, iMesh % Name ) 
    MeshName = GetString( Params,'Mesh Name',Found )
    IF(Found) THEN
      nlen = LEN_TRIM(MeshName)
      IF( MeshName(1:nlen) /= iMeshName(1:nlen) ) THEN
        iMesh => iMesh % next
        CYCLE
      END IF
    END IF

        
    CALL SetCurrentMesh( Model, iMesh )
    ModelVariables => Model % Variables
    Model % Variables => iMesh % variables 

    IF( .NOT. ListSet ) THEN
      CALL Info('ResultOutputSolver','Creating list for saving - if not present')
      CALL CreateListForSaving( Model, Params,.TRUE. )    
      ListSet = .TRUE.
    ELSE IF( MeshDim /= Model % Mesh % MeshDim .OR. (iMeshName(1:nlen) /= ListMeshName(1:nlen))) THEN
      CALL Info('ResultOutputSolver','Recreating list for saving')
      CALL CreateListForSaving( Model, Params,.TRUE.,.TRUE.)
    END IF

    MeshDim = Model % Mesh % MeshDim
    nlen = StringToLowerCase( ListMeshName, iMesh % Name)

    ! In case there are multiple mesh levels one may also save coarser ones
    !----------------------------------------------------------------------
    Mesh => iMesh
    DO i=1,MeshLevel
      Mesh => Mesh % Parent
      IF (.NOT.ASSOCIATED(Mesh)) EXIT
    END DO
    IF ( ASSOCIATED(Mesh)) THEN

    CALL SetCurrentMesh( Model, Mesh )
    Model % Variables => Mesh % variables 
    SomeMeshSaved = .TRUE.

    
    IF( SaveGid ) THEN
      CALL Info( 'ResultOutputSolver','Saving in GiD format' )    
      CALL GiDOutputSolver( Model,Solver,dt,TransientSimulation )
    END IF
    IF( SaveGmsh ) THEN
      CALL Info( 'ResultOutputSolver','Saving in new gmsh format' )        
      CALL GmshOutputSolver( Model,Solver,dt,TransientSimulation )
    END IF
    IF( SaveVTK ) THEN
      CALL Info( 'ResultOutputSolver','Saving in legacy VTK format' )            
      CALL VtkOutputSolver( Model,Solver,dt,TransientSimulation )
    END IF
    IF( SaveVTU ) THEN
      CALL Info( 'ResultOutputSolver','Saving in unstructured VTK XML (.vtu) format' )               
      CALL VtuOutputSolver( Model,Solver,dt,TransientSimulation )
    END IF
    IF( SaveOpenDx ) THEN
      CALL Info( 'ResultOutputSolver','Saving in OpenDX format' )                   
      CALL DXOutputSolver( Model,Solver,dt,TransientSimulation )
    END IF
    IF( SaveEP ) THEN
      CALL Info( 'ResultOutputSolver','Saving in ElmerPost format' )                   
      CALL ElmerPostOutputSolver( Model,Solver,dt,TransientSimulation )
    END IF

    CALL Info( 'ResultOutputSolver', '-------------------------------------')
    END IF

    iMesh => iMesh % Next
  END DO

  IF( .NOT. SomeMeshSaved ) THEN
    OutputCount = OutputCount - 1
  END IF
  Model % Variables => ModelVariables

  SubroutineVisited = .TRUE.

END SUBROUTINE ResultOutputSolver
