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
      CHARACTER (len=MAX_NAME_LEN) :: FName
      CHARACTER (len=MAX_NAME_LEN) :: VarName,TVarName,T2VarName
      CHARACTER (len=MAX_NAME_LEN) :: Txt

      INTEGER :: VarType
      INTEGER :: varid,nvals,ntime,ncid,ndims,TVarID
      INTEGER :: NetCDFstatus
      INTEGER :: dimids(2) 
      REAL(KIND=dp), ALLOCATABLE :: Values(:)
      REAL(KIND=dp) :: Time
      INTEGER :: TimeIndex,TimePoint
      INTEGER :: EIndex,NIndex,VarIndex
      LOGICAL :: Parallel,Found,VarExist
      INTEGER, SAVE :: VisitedTimes=0
      LOGICAL, POINTER :: UnFoundNodes(:) => NULL()
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

! check if this is a paralell run
      Parallel=(ParEnv % PEs > 1)

! get mesh
      ThisMesh => GetMesh(Solver)

! target mesh if interpolation required; copy form Mesh2MeshSolver
      TargetMesh => NULL()
      i = ListGetInteger( SolverParams,'Target Mesh Solver Index',Found )
      IF( Found ) THEN
       ! Target mesh solver is explicitly given
       TargetMesh => CurrentModel % Solvers(i) % Mesh
       IF( ASSOCIATED( TargetMesh ) ) THEN
         CALL Info(SolverName,'Using target mesh as the mesh of Solver '//TRIM(I2S(i)),Level=8)
       ELSE
        CALL Fatal(SolverName,'Target Mesh for Solver not associated: '//TRIM(I2S(i)))
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


      ! get time index
      VisitedTimes = VisitedTimes + 1
      IF( ListGetLogical( SolverParams, "Is Time Counter", Found ) ) THEN
        TimePoint = VisitedTimes
      ELSE
        Time = ListGetConstReal( SolverParams, "Time Point", Found )
        IF (.NOT.Found) Time=GetTime()
        dt = GetTimeStepSize()
        TimePoint = floor(time-dt/2) + 1
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

            TimeIndex = max(1,min(TimePoint,ntime))

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

        Var => VariableGet( ThisMesh % Variables,TRIM(TVarName),ThisOnly=.TRUE.,UnFoundFatal=.TRUE.)
        VarType=Var % TYPE

        ! special cases.... time do not seems to be a gloabl variable by
        ! default
        IF ( Var % Name  == 'time') VarType=Variable_global

        SELECT CASE(VarType)
         CASE(Variable_global)
           ! for a global variable nvals should be the time dimension...
           TimeIndex = max(1,min(TimePoint,nvals))
           Var % Values(1)=Values(TimeIndex)
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
             IF (i.GT.nvals) &
                CALL FATAL(SolverName,"Too many nodes "//TRIM(VarName))
             Var%Values(k)=Values(NIndex)
           END DO

         CASE DEFAULT
            CALL FATAL(SolverName,"wrong var type "//TRIM(VarName))

        END SELECT

        DEALLOCATE(Values)

        ! close file
        NetCDFstatus=nf90_close(ncid)

        WRITE(Txt,'(A,I0)') 'Interpolated target Variable ',VarIndex
        T2VarName = ListGetString(SolverParams,TRIM(Txt),Found)
        IF (Found) THEN
          IF (VarType.NE.Variable_on_nodes) &
            CALL FATAL(SolverName,"Interpolation possible only with nodal variables")
          IF (.NOT.ASSOCIATED(TargetMesh)) &
            CALL FATAL(SolverName,"Target mesh not found")

          Var % Name = T2VarName
          Var % NameLen = LENTRIM( T2VarName )
          CALL InterpolateMeshToMesh( ThisMesh, &
            TargetMesh, Var, TargetMesh % Variables,&
            Projector=Projector, UnfoundNodes=UnfoundNodes)

          ! better to do something if there is unfoundnodes
          nf = COUNT(UnfoundNodes)
          IF (nf.GT.0) &
            CALL FATAL(SolverName,"There is unfound nodes"//TRIM(I2S(nf)))

          ! Validate the mapped variables
          pVar => VariableGet( TargetMesh % Variables, TRIM(T2VarName), ThisOnly = .TRUE.,UnFoundFatal=.TRUE. )
          pVar % Valid = .TRUE.
          pVar % ValuesChanged = .TRUE.

          ! back to initial name
          Var % Name = TVarName
          Var % NameLen = LENTRIM( TVarName )
          Var => VariableGet( ThisMesh % Variables,TRIM(TVarName),ThisOnly=.TRUE.,UnFoundFatal=.TRUE.)
        END IF

        VarIndex=VarIndex+1
        WRITE(Txt,'(A,I0)') 'Variable Name ',VarIndex
        VarName = ListGetString(SolverParams,TRIM(Txt),VarExist)
      END DO


      END SUBROUTINE UGRIDDataReader
