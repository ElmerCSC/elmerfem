!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
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
! ******************************************************************************
! *
! *  Authors: F. Gillet-Chaulet (IGE-France)
! *  Web:     http://elmerice.elmerfem.org
! *  Original Date: 05/2022
! * 
! *****************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE UGRIDDataReader_init( Model,Solver,dt,TransientSimulation )
      USE DefUtils
      IMPLICIT NONE
      !------------------------------------------------------------------------------
      TYPE(Solver_t), TARGET :: Solver
      TYPE(Model_t) :: Model
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation
      !--------------------------------------------------------------------------
      CHARACTER(LEN=MAX_NAME_LEN) :: Name
      TYPE(ValueList_t), POINTER :: SolverParams
      LOGICAL :: GotIt

      SolverParams => Solver % Values

      ! A solver variable is required 
      ! if we need to create a variable here using Exported Variable; 
      Name = ListGetString( SolverParams, 'Equation',GotIt)
      IF( .NOT. ListCheckPresent( SolverParams,'Variable') ) THEN
        IF( ListCheckPresent( SolverParams,'Exported Variable 1') ) THEN
          CALL ListAddString( SolverParams,'Variable',&
           '-nooutput '//TRIM(Name)//'_var')
        END IF
      ENDIF

      IF(.NOT. ListCheckPresent(SolverParams,'Optimize Bandwidth')) &
        CALL ListAddLogical(SolverParams,'Optimize Bandwidth',.FALSE.)

      END SUBROUTINE UGRIDDataReader_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      SUBROUTINE UGRIDDataReader( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
      USE DefUtils
      USE SolverUtils
      USE NETCDF
      IMPLICIT NONE
!------------------------------------------------------------------------------
      TYPE(Solver_t), TARGET :: Solver
      TYPE(Model_t) :: Model
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
      CHARACTER(*), PARAMETER :: SolverName="UGRIDDataReader"
      TYPE(ValueList_t), POINTER :: SolverParams
      TYPE(Mesh_t), POINTER :: ThisMesh, TargetMesh,Mesh
      TYPE(Projector_t), POINTER :: Projector
      TYPE(Variable_t),POINTER :: Var,pVar
      TYPE(Element_t), POINTER :: Element
      INTEGER :: i,k
      INTEGER :: NN,nf
      CHARACTER (len=MAX_STRING_LEN) :: FName
      CHARACTER (len=MAX_NAME_LEN) :: VarName,TVarName,T2VarName
      CHARACTER (len=MAX_NAME_LEN) :: Txt
      CHARACTER (len=MAX_NAME_LEN) :: MeshName

      INTEGER :: VarType
      INTEGER :: varid,nvals,ntime,ncid,ndims,TVarID
      INTEGER :: NetCDFstatus
      INTEGER :: dimids(2) 
      REAL(KIND=dp), ALLOCATABLE :: Values(:)
      REAL(KIND=dp) :: Time
      INTEGER :: TimeIndex,TimePoint,TimeOffset
      INTEGER :: EIndex,NIndex,VarIndex
      LOGICAL :: Parallel,Found,VarExist
      INTEGER, SAVE :: VisitedTimes=0
      LOGICAL, POINTER :: UnFoundNodes(:) => NULL()
      LOGICAL :: UnFoundNodesFatal
      LOGICAL :: DoInterp
      REAL(KIND=dp) :: MinDist
      INTEGER :: MinNode
      REAL(KIND=dp), ALLOCATABLE,SAVE :: coord(:)
      INTEGER :: dim
      !------------------------------------------------------------------------------
      INTERFACE
        SUBROUTINE InterpolateMeshToMesh( OldMesh, NewMesh, OldVariables, &
            NewVariables, UseQuadrantTree, Projector, MaskName, UnfoundNodes )
          USE Types
          TYPE(Variable_t), POINTER, OPTIONAL :: OldVariables, NewVariables
          TYPE(Mesh_t), TARGET  :: OldMesh, NewMesh
          LOGICAL, OPTIONAL :: UseQuadrantTree
          LOGICAL, POINTER, OPTIONAL :: UnfoundNodes(:)
          CHARACTER(LEN=*),OPTIONAL :: MaskName
          TYPE(Projector_t), POINTER, OPTIONAL :: Projector
        END SUBROUTINE InterpolateMeshToMesh
      END INTERFACE
!------------------------------------------------------------------------------
! get parameters
      SolverParams => GetSolverParams()

! get mesh
      ThisMesh => GetMesh(Solver)

! check if this is a parallel run
      Parallel=(ParEnv % PEs > 1) .AND. ( .NOT. ThisMesh % SingleMesh ) 

! check if a mesh ha been defined for this solver; in this case will
! will interpolate nodal variables to the target mesh
      MeshName = ListGetString(SolverParams,'mesh',DoInterp)
      IF (DoInterp) THEN
        dim = CoordinateSystemDimension()
        IF (.NOT.ALLOCATED(coord)) ALLOCATE(coord(dim))
        ! target mesh if interpolation required
        TargetMesh => NULL()
        i = ListGetInteger( SolverParams,'Target Mesh Solver Index',Found )
        IF( Found ) THEN
         ! Target mesh solver is explicitly given
         TargetMesh => CurrentModel % Solvers(i) % Mesh
         IF( ASSOCIATED( TargetMesh ) ) THEN
           CALL Info(SolverName,'Using target mesh as the mesh of Solver '//I2S(i),Level=8)
         ELSE
          CALL Fatal(SolverName,'Target Mesh for Solver not associated: '//I2S(i))
         END IF
        ELSE
          ! Otherwise use the 1st mesh that is not this old data mesh
          Mesh => CurrentModel % Meshes
          DO WHILE( ASSOCIATED(Mesh) )
            IF( .NOT. ASSOCIATED( Mesh, ThisMesh ) ) THEN
              TargetMesh => Mesh
              EXIT
            END IF
            Mesh => Mesh % Next
          END DO
          IF( ASSOCIATED( TargetMesh ) ) THEN
           CALL Info(SolverName,'Using target mesh as the first mesh different from this mesh',Level=8)
          END IF
        END IF
      END IF
      DoInterp=DoInterp.AND.ASSOCIATED( TargetMesh )

      ! should we kill the run if nodes are not found at interpolation?
      UnFoundNodesFatal = ListGetLogical(SolverParams,'UnFoundNodes Fatal',Found )
      IF (.NOT.Found) UnFoundNodesFatal = .TRUE.


      ! get time index
      VisitedTimes = VisitedTimes + 1
      IF( ListGetLogical( SolverParams, "Is Time Counter", Found ) ) THEN
        TimeOffset=ListGetInteger( SolverParams, "Time Counter start", Found )
        IF (Found) THEN
          TimePoint = VisitedTimes + TimeOffset - 1
        ELSE
          TimePoint = VisitedTimes
        ENDIF
      ELSE
        TimePoint = ListGetInteger( SolverParams, "Time Index", Found )
        IF (.NOT.Found) THEN
          Time = ListGetCReal( SolverParams, "Time Point", Found )
          IF (.NOT.Found) Time=GetTime()
          dt = GetTimeStepSize()
          TimePoint = floor(time-dt/2) + 1
        END IF
      END IF


      VarIndex=1
      VarName = ListGetString(SolverParams,'Variable Name 1',UnFoundFatal=.TRUE.)
      VarExist=.TRUE.
      DO WHILE (VarExist)

         WRITE(Txt,'(A,I0)') 'File Name ',VarIndex
         FName = ListGetString(SolverParams,TRIM(Txt),Found)
         IF (.NOT.Found) &
            FName = ListGetString(SolverParams,'File Name',UnFoundFatal=.TRUE.)

         write(Message,'(a,a)') 'File name: ',Trim(FName)
        CALL INFO(SolverName,Message,Level=4)

!! open NETCDF File
      NetCDFstatus = NF90_OPEN(trim(FName),NF90_NOWRITE,ncid)
      IF ( NetCDFstatus /= NF90_NOERR ) &
         CALL FATAL(SolverName,"file open failed")

        ! get variable ID
        NetCDFstatus = nf90_inq_varid(ncid,trim(VarName),TVarId)
        IF ( NetCDFstatus /= NF90_NOERR ) &
          CALL FATAL(SolverName,"unable to get varid "//TRIM(VarName))

        ! variable dimensions
        NetCDFstatus = nf90_inquire_variable(ncid, TVarId, ndims=ndims)
        IF ( NetCDFstatus /= NF90_NOERR ) &
          CALL FATAL(SolverName,"unable to get variable dimensions "//TRIM(VarName))

        ! Get the dim ids
        NetCDFstatus = nf90_inquire_variable(ncid, Tvarid, dimids = dimids(1:ndims))
        IF ( NetCDFstatus /= NF90_NOERR ) &
           CALL FATAL(SolverName,"unable to get dimids "//TRIM(VarName))
       
        ! here we assume that the first dim is the size of the variable
        ! and that the second should be time
        NetCDFstatus = nf90_inquire_dimension(ncid,dimids(1),len=nvals)
        IF ( NetCDFstatus /= NF90_NOERR ) &
           CALL FATAL(SolverName,"unable to get nvals"//TRIM(VarName))


        WRITE(Message,'(A,I0)') &
                    TRIM(VarName)//' ndims: ',ndims
        CALL INFO(SolverName,Trim(Message),level=4)
        WRITE(Message,'(A,I0)') &
                    TRIM(VarName)//' dimension: ',nvals
        CALL INFO(SolverName,Trim(Message),level=4)

        allocate(Values(nvals))

        SELECT CASE(ndims)
         CASE(1)
           NetCDFstatus=nf90_get_var(ncid,TVarId,Values,start=(/1/),count=(/nvals/))
           IF ( NetCDFstatus /= NF90_NOERR ) &
              CALL FATAL(SolverName,"unable to get variable "//TRIM(VarName))

            WRITE(Message,'(A)') &
                    'reading constant var '//TRIM(VarName)
            CALL INFO(SolverName,Trim(Message),level=4)
         CASE(2)
            NetCDFstatus = nf90_inquire_dimension(ncid,dimids(2),len=ntime)
            IF ( NetCDFstatus /= NF90_NOERR ) &
              CALL FATAL(SolverName,"unable to get ntime "//TRIM(VarName))

            IF (TimePoint.LT.0) THEN
              TimeIndex = ntime
            ELSE
              TimeIndex = max(1,min(TimePoint,ntime))
            ENDIF

            WRITE(Message,'(A,I0)') &
                    TRIM(VarName)//', reading time step: ',TimeIndex
            CALL INFO(SolverName,Trim(Message),level=4)

            NetCDFstatus=nf90_get_var(ncid,TVarId,Values,start=(/1,TimeIndex/),count=(/nvals/))
            IF ( NetCDFstatus /= NF90_NOERR ) &
               CALL FATAL(SolverName,"unable to get variable "//TRIM(VarName))
         CASE DEFAULT
           CALL FATAL(SolverName,"wrong ndims "//TRIM(VarName))
        END SELECT
        
     ! get variable
        WRITE(Txt,'(A,I0)') 'Target Variable ',VarIndex
        TVarName = ListGetString(SolverParams,TRIM(Txt),Found)
        IF (.NOT.Found) TVarName=TRIM(VarName)

        Var => VariableGet( ThisMesh % Variables,TRIM(TVarName),UnFoundFatal=.TRUE.)
        VarType=Var % TYPE

        ! special cases.... time do not seems to be a global variable by
        ! default
        IF ( Var % Name  == 'time') VarType=Variable_global

        SELECT CASE(VarType)
         CASE(Variable_global)
           ! for a global variable nvals should be the time dimension...
           TimeIndex = max(1,min(TimePoint,nvals))
           Var % Values(1)=Values(TimeIndex)
           IF (ASSOCIATED(TargetMesh)) THEN
             pVar => VariableGet( TargetMesh % Variables,TRIM(TVarName),ThisOnly=.TRUE.,UnFoundFatal=.TRUE.)
             pVar % Values(1)=Values(TimeIndex)
           END IF
         CASE(Variable_on_elements)

           NN = GetNOFActive()
           DO i=1,NN
             Element => GetActiveElement(i)
             IF (Parallel) THEN 
               EIndex=Element % GElementIndex
             ELSE
               EIndex=Element % ElementIndex
             ENDIF
             IF (EIndex.GT.nvals) &
                     CALL FATAL(SolverName,"EIndex larger than nvals "//TRIM(VarName))
             IF (ASSOCIATED(Var % Perm)) THEN
               k=Var % Perm(Element % ElementIndex)
             ELSE
               k=Element % ElementIndex
             ENDIF
             IF (k==0) CYCLE
             Var % Values(k)=Values(EIndex)
           END DO

         CASE(Variable_on_nodes)
           NN = Solver % Mesh % NumberOfNodes
           DO i=1,NN
             IF (Parallel) THEN 
               NIndex=Solver % Mesh % ParallelInfo % GlobalDOFs(i)
             ELSE
              NIndex=i
             ENDIF
             IF (ASSOCIATED(Var % Perm)) THEN
               k=Var % Perm(i)
             ELSE
               k=i
             ENDIF
             IF (k==0) CYCLE
             !IF NIndex>nvals assume the mesh is structured
             ! and nodenumbering is  by layers
             IF (NIndex.GT.nvals) THEN
                     NIndex=MOD(NIndex,nvals)
                     IF (NIndex.EQ.0) NIndex=nvals
             ENDIF
             IF ((NIndex.GT.nvals).OR.(NIndex.LT.1)) &
                CALL FATAL(SolverName,"Wrong NIndex for "//TRIM(VarName)//" "//I2S(NIndex))

             Var%Values(k)=Values(NIndex)
           END DO

         CASE DEFAULT
            CALL FATAL(SolverName,"wrong var type "//TRIM(VarName))

        END SELECT

        DEALLOCATE(Values)

        ! close file
        NetCDFstatus=nf90_close(ncid)

        IF (DoInterp.AND.(VarType.EQ.Variable_on_nodes)) THEN

           ! Rename variable in this mesh with the name in the target mesh
	   ! to do the interpolation
	   WRITE(Txt,'(A,I0)') 'Target Mesh Variable ',VarIndex
           T2VarName = ListGetString(SolverParams,TRIM(Txt),Found)
           IF (.NOT.Found) T2VarName=TRIM(VarName)
           Var % NameLen = StringToLowerCase( Var % Name,T2VarName)

          CALL InterpolateMeshToMesh( ThisMesh, &
                  TargetMesh, Var, TargetMesh % Variables,&
                  UnfoundNodes=UnfoundNodes)

          ! Validate the mapped variables
          pVar => VariableGet( TargetMesh % Variables, TRIM(Var%Name), ThisOnly = .TRUE.,UnFoundFatal=.TRUE. )
          pVar % Valid = .TRUE.
          pVar % ValuesChanged = .TRUE.

          ! better to do something if there is unfoundnodes
          nf = COUNT(UnfoundNodes)
          IF (nf.GT.0) THEN
            IF (UnFoundNodesFatal) THEN
              CALL FATAL(SolverName,"There is unfound nodes : "//I2S(nf))
            ELSE
              CALL WARN(SolverName,TRIM(TVarName)//"; there is "//I2S(nf)//" unfound nodes; get closest node in input mesh")
              IF (Parallel) &
                CALL FATAL(SolverName,"dealing with unfound nodes only for serial meshes; add -single ")

              DO i=1,TargetMesh % NumberOfNodes
                IF (UnfoundNodes(i)) THEN
                   coord(1)=TargetMesh % Nodes % x (i)
                   IF (dim.EQ.2) coord(2)=TargetMesh % Nodes % y (i)
                   IF (dim.EQ.3) coord(3)=TargetMesh % Nodes % z (i)
                   CALL FindClosestNode(ThisMesh,coord,MinDist,MinNode,Parallel)
                   pVar % Values (pVar % Perm(i)) = Var%Values(Var % Perm(MinNode))
                ENDIF
              END DO

            END IF
          END IF



        END IF

        VarIndex=VarIndex+1
        WRITE(Txt,'(A,I0)') 'Variable Name ',VarIndex
        VarName = ListGetString(SolverParams,TRIM(Txt),VarExist)
      END DO


      END SUBROUTINE UGRIDDataReader
