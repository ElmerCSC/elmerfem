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
! *  Authors: f. Gillet-Chaulet (IGE, Grenoble,France)
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! * 
! *****************************************************************************
SUBROUTINE Adjoint_CostDiscSolver_init0(Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  CHARACTER(LEN=MAX_NAME_LEN) :: Name
  
  Name = ListGetString( Solver % Values, 'Equation',UnFoundFatal=.TRUE.)
  CALL ListAddNewString( Solver % Values,'Variable',&
          '-nooutput '//TRIM(Name)//'_var')
  CALL ListAddLogical(Solver % Values, 'Optimize Bandwidth',.FALSE.)

  CALL ListAddInteger(Solver % Values, 'Nonlinear System Norm Degree',0)

END SUBROUTINE Adjoint_CostDiscSolver_init0
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE Adjoint_CostDiscSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Discrete cost function evaluated as a sum of the differences at observations points
! See documentation under ELMER_SRC/elmerice/Solvers/Documentation/Adjoint_CostDiscSolver.md
!******************************************************************************
  USE ParticleUtils
  USE GeneralUtils
  USE DefUtils
  USE Interpolation
#ifdef HAVE_NETCDF
  USE Netcdf
#endif
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!  
  TYPE(Mesh_t),POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: SolverParams

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName="Adjoint_CostDiscSolver"
  CHARACTER(LEN=MAX_NAME_LEN),SAVE :: CostSolName,VariableName
  CHARACTER(LEN=MAX_NAME_LEN) :: FVarName
  CHARACTER(LEN=MAX_NAME_LEN):: ObsFileName
  CHARACTER(LEN=MAX_NAME_LEN):: SideFile,SideParFile
  CHARACTER(LEN=MAX_NAME_LEN),SAVE :: CostFile
  CHARACTER(:), ALLOCATABLE :: OutPutDirectory
  CHARACTER(*), PARAMETER :: DefaultCostFile = 'CostOfT.dat'

  TYPE(Variable_t), POINTER :: TimeVar ! Time
  TYPE(Variable_t), POINTER :: CostVar ! Cost Value
  TYPE(Variable_t), POINTER :: Variable ! Active variable
  TYPE(Variable_t), POINTER :: VbSol   ! The sensitivity
  REAL(KIND=dp), POINTER :: Vb(:),Values(:)
  INTEGER, POINTER :: VbPerm(:),Perm(:)
  INTEGER,SAVE :: VDOFs

  ! Variable related to elements
  TYPE(Element_t),POINTER ::  Element
  TYPE(Nodes_t),SAVE :: ElementNodes
  INTEGER :: ElementIndex
  INTEGER, POINTER :: NodeIndexes(:)
  REAL(KIND=dp),ALLOCATABLE,SAVE :: Basis(:), dBasisdx(:,:)
  REAL(KIND=dp) :: SqrtElementMetric
  INTEGER :: n
  LOGICAL :: stat

  LOGICAL,SAVE :: Firsttime=.true.,Parallel
  LOGICAL,SAVE :: FirstRound=.True.
  LOGICAL,SAVE :: SAVE_USED_DATA=.False.
  LOGICAL,SAVE :: HaveSTD=.FALSE.
  LOGICAL,SAVE :: CompNorm=.FALSE.
  LOGICAL :: Found
  LOGICAL :: ParallelFile,ProcessedFile,ActiveNumbering

  INTEGER :: CoordDIM ! Coordinate system dimension for the observations points
  INTEGER,SAVE :: VarDIM   ! Dimension of the observed Variable
  REAL(KIND=dp),SAVE :: Lambda

  integer :: i,j,s
  INTEGER :: k
  real(kind=dp) :: StdMin
  real(kind=dp) :: Cost,Cost_S
  real(kind=dp) :: rms,rms_s
  real(kind=dp) :: misfit
  REAL(Kind=dp) :: VmNorm,VoNorm,StdNorm,VmNormb,um,umb
  real(kind=dp) :: Coord(3),UVW(3)
  REAL(KIND=dp) :: V(3)
  INTEGER,SAVE :: NTOT=-1,NTOT_S
  INTEGER :: AI
  CHARACTER*10 :: date,temps
  INTEGER,PARAMETER :: IO=12
  INTEGER :: ok
  LOGICAL :: WarnActive
  LOGICAL :: BoundarySolver

  ! Variables with obs.
  REAL(KIND=dp),ALLOCATABLE,SAVE :: xobs(:,:),Vobs(:,:),StdObs(:,:)
  INTEGER,ALLOCATABLE,SAVE :: InElement(:)
  INTEGER,SAVE :: nobs

  ! Variables related to netcdf
  CHARACTER(LEN=MAX_NAME_LEN) :: Xdim,Ydim,VarName
  REAL(KIND=dp),ALLOCATABLE :: Fillv(:)
  REAL(KIND=dp),ALLOCATABLE :: xx(:),yy(:),ObsArray(:,:,:),StdObsArray(:,:,:)
  INTEGER :: XMinIndex,XMaxIndex,YMinIndex,YMaxIndex,nx,ny
  INTEGER :: NetcdfStatus,varid,ncid,ierr
  LOGICAL :: NETCDFFormat
  REAL(KIND=dp) :: xmin,ymin,xmax,ymax


  SolverParams => GetSolverParams()

  !! check if we are on a boundary or in the bulk
  BoundarySolver = ( Solver % ActiveElements(1) > Model % Mesh % NumberOfBulkElements )
  IF (BoundarySolver) THEN
    CoordDIM = CoordinateSystemDimension() - 1
  ELSE
    CoordDIM = CoordinateSystemDimension()
  ENDIF

! Get min/max coordinates of current mesh
  xmin=MINVAL(Solver%Mesh%Nodes%x(:))
  xmax=MAXVAL(Solver%Mesh%Nodes%x(:))
  IF (CoordDIM.GE.2) THEN
    ymin=MINVAL(Solver%Mesh%Nodes%y(:))
    ymax=MAXVAL(Solver%Mesh%Nodes%y(:))
  ENDIF


  If (Firsttime) then
    N = model % MaxElementNodes
    allocate(ElementNodes % x(N), ElementNodes % y(N), ElementNodes % z(N))
    allocate(Basis(N), dBasisdx(N,3))

!!!!!!! Check for parallel run 
    Parallel = .FALSE.
    IF ( ASSOCIATED( Solver % Matrix % ParMatrix ) ) THEN
            IF ( Solver %  Matrix % ParMatrix % ParEnv % PEs > 1 )  THEN
                    Parallel = .TRUE.
            END IF
    END IF

   Lambda =  GetConstReal( SolverParams,'Lambda', Found)
   IF(.NOT.Found) THEN
      CALL WARN(SolverName,'Keyword >Lambda< not found  in section >Solver<')
      CALL WARN(SolverName,'Taking default value Lambda=1.0')
      Lambda = 1.0_dp
      CALL ListAddConstReal( SolverParams,'Lambda',Lambda)
   End if

   ! error is the difference of the norm?
   CompNorm = ListGetLogical(SolverParams,'Compute norm difference', Found)
   IF(.NOT.Found) CompNorm=.FALSE.

   !! Give information about standard error
   HaveSTD=ListGetLogical(SolverParams,'Use standard error',Found)
   IF (.NOT.Found) HaveSTD=.False.
   IF (HaveSTD) THEN
      StdMin=ListGetConstReal(SolverParams,'standard error Min Value',Found)
      IF (.NOT.Found) StdMin=100*AEPS
   END IF

!!!!!!!!!!! get Solver Variables

  CostFile = ListGetString(Solver % Values,'Cost Filename',Found )
  IF (.NOT. Found) CostFile = DefaultCostFile
  CALL DATE_AND_TIME(date,temps)
  If (Parallel) then
     if (ParEnv % MyPe.EQ.0) then
        OPEN (IO, FILE=CostFile)
             write(IO,1000) date(5:6),date(7:8),date(1:4),temps(1:2),temps(3:4),temps(5:6)
             write(IO,1001) Lambda
             IF (HaveSTD) THEN
               write(IO,'(A)') '# iter, Lambda*J0, sqrt(2*J0/Nobs),rms'
             ELSE
               write(IO,'(A)') '# iter, Lambda*J0, sqrt(2*J0/Nobs)'
             END IF
        CLOSE(IO)
     End if
  Else
       OPEN (IO, FILE=CostFile)
             write(IO,1000) date(5:6),date(7:8),date(1:4),temps(1:2),temps(3:4),temps(5:6)
             write(IO,1001) Lambda
             IF (HaveSTD) THEN
               write(IO,'(A)') '# iter, Lambda*J0, sqrt(2*J0/Nobs),rms'
             ELSE
               write(IO,'(A)') '# iter, Lambda*J0, sqrt(2*J0/Nobs)'
             END IF
       CLOSE(IO)
  End if

   CostSolName =  GetString( SolverParams,'Cost Variable Name', Found)
   IF(.NOT.Found) THEN
         CALL WARN(SolverName,'Keyword >Cost Variable Name< not found  in section >Solver<')
         CALL WARN(SolverName,'Taking default value >CostValue<')
         WRITE(CostSolName,'(A)') 'CostValue'
   END IF

   VariableName =  ListGetString( SolverParams,'Observed Variable Name', UnFoundFatal=.TRUE.)
   Variable => VariableGet( Solver % Mesh % Variables, Trim(VariableName) ,UnFoundFatal=.TRUE. )
   VDOFs=Variable%DOFs

   VarDIM=GetInteger(SolverParams ,'Observed Variable Dimension',Found)
   IF (.NOT.Found) THEN
    VarDIM=MIN(VDOFs,CoordDIM)
    WRITE(message,*) 'Keyword >Observed Variable Dimension< not found, assume VarDIM = ',VarDIM
    CALL WARN(SolverName,message)
   ENDIF

!! Get the obs
   ObsFileName =  ListGetString( SolverParams,'Observation File Name', UnFoundFatal=.TRUE.)

!! I the file netcf?
   k = INDEX(ObsFileName ,'.nc' )
   NETCDFFormat = ( k /= 0 )

   IF (NETCDFFormat) then

#ifdef HAVE_NETCDF

      IF (CoordDIM.NE.2) CALL FATAL(Trim(SolverName),'Netcdf only supported for 2D')

      CALL INFO(Trim(SolverName),'Data File is in netcdf format', Level=5)

      NetCDFstatus = NF90_OPEN(trim(ObsFileName),NF90_NOWRITE,ncid)
      IF ( NetCDFstatus /= NF90_NOERR ) &
          CALL FATAL(Trim(SolverName), 'Unable to open NETCDF File')

      Xdim=ListGetString( SolverParams,"X Dim Name", Found )
      if (.NOT.Found) Xdim='x'

      NetCDFstatus = nf90_inq_dimid(ncid, trim(Xdim) , varid)
      IF ( NetCDFstatus /= NF90_NOERR ) &
          CALL FATAL(Trim(SolverName),'Unable to  get netcdf x-dim Id')

      NetCDFstatus = nf90_inquire_dimension(ncid,varid,len=nx)
      IF ( NetCDFstatus /= NF90_NOERR ) &
           CALL FATAL(Trim(SolverName), 'Unable to  get netcdf nx')

      Ydim=ListGetString( SolverParams,"Y Dim Name", Found )
      if (.NOT.Found) Ydim='y'

      NetCDFstatus = nf90_inq_dimid(ncid, trim(Ydim) , varid)
      IF ( NetCDFstatus /= NF90_NOERR ) &
          CALL FATAL(Trim(SolverName),'Unable to  get netcdf y-dim Id')

      NetCDFstatus = nf90_inquire_dimension(ncid,varid,len=ny)
      IF ( NetCDFstatus /= NF90_NOERR ) &
           CALL FATAL(Trim(SolverName), 'Unable to  get netcdf ny')

      ALLOCATE(xx(nx),yy(ny))

      !! Get X and Y variable
      Xdim=ListGetString( SolverParams, "X Var Name", Found )
      if (.NOT.Found) Xdim='x'
      NetCDFstatus = nf90_inq_varid(ncid,trim(Xdim),varid)
      IF ( NetCDFstatus /= NF90_NOERR ) &
           CALL FATAL(Trim(SolverName),'Unable to get netcdf x-variable id')

      NetCDFstatus = nf90_get_var(ncid, varid,xx)
      IF ( NetCDFstatus /= NF90_NOERR ) &
           CALL FATAL(Trim(SolverName), 'Unable to get netcdf x-variable ')

      Ydim=ListGetString( SolverParams, "Y Var Name", Found )
      if (.NOT.Found) Ydim='y'
      NetCDFstatus = nf90_inq_varid(ncid,trim(Ydim),varid)
      IF ( NetCDFstatus /= NF90_NOERR ) &
           CALL FATAL(Trim(SolverName),'Unable to get netcdf x-variable id')

      NetCDFstatus = nf90_get_var(ncid, varid,yy)
      IF ( NetCDFstatus /= NF90_NOERR ) &
           CALL FATAL(Trim(SolverName), 'Unable to get netcdf x-variable ')

      !! Check that there is data within the domain
      IF ((MAXVAL(xx).LT.xmin).OR.(MINVAL(xx).GT.xmax)&
            .OR.(MAXVAL(yy).LT.ymin).OR.(MINVAL(yy).GT.ymax)) &
             CALL Fatal(Trim(SolverName), &
                        'No data within model domain')
      !!! get the index of values that are within the mesh bounding box
      CALL MinMaxIndex(xx,nx,xmin,xmax,XMinIndex,XMaxIndex)
      CALL MinMaxIndex(yy,ny,ymin,ymax,YMinIndex,YMaxIndex)


     !! get ux and uy
      nx=XMaxIndex-XMinIndex+1
      ny=YMaxIndex-YMinIndex+1
      ALLOCATE(ObsArray(nx,ny,VarDIM),fillv(VarDIM))
      IF (HaveSTD) ALLOCATE(StdObsArray(nx,ny,VarDIM))

      DO k=1,VarDIM
        IF (VarDIM==1) THEN
          VarName=ListGetString( SolverParams, "Netcdf Var Name", Found )
          IF (.NOT.Found) VarName=TRIM(VariableName)
        ELSE
          VarName=ListGetString( SolverParams, TRIM(ComponentName("Netcdf Var",k)) // " Name", Found )
          IF (.NOT.Found) THEN
              CALL WARN(SolverName,'keyword <' // &
                      TRIM(ComponentName("Netcdf Var",k)) // " Name" // ' not found, taking default')
              select case (k)
                case (1)
                  VarName='vx'
                case (2)
                  VarName='vy'
              end select
          ENDIF
        ENDIF
        NetCDFstatus = nf90_inq_varid(ncid,trim(VarName),varid)
        IF ( NetCDFstatus /= NF90_NOERR ) &
           CALL Fatal(Trim(SolverName),'Unable to get netcdf variable id for ' // TRIM(VarName))
        NetCDFstatus = nf90_get_var(ncid, Varid, ObsArray(:, :,k), &
                            start = (/ XMinIndex, YMinIndex /),     &
                            count = (/ nx,ny/))
        IF ( NetCDFstatus /= NF90_NOERR ) &
            CALL Fatal(Trim(SolverName),'Unable to get netcdf variable' // TRIM(VarName))
        NetCDFstatus = nf90_get_att(ncid, varid,"_FillValue",fillv(k))
        IF ( NetCDFstatus /= NF90_NOERR ) &
           CALL FATAL(Trim(SolverName),TRIM(VarName) // 'has no attribute _FillValue')
      END DO

      IF (HaveSTD) THEN
        DO k=1,VarDIM
          IF (VarDIM==1) THEN
            VarName=ListGetString( SolverParams, "Netcdf Std Var Name", UnFoundFatal=.True. )
          ELSE
            VarName=ListGetString( SolverParams, TRIM(ComponentName("Netcdf Std Var",k)) // " Name", Found )
            IF (.NOT.Found) THEN
              CALL WARN(SolverName,'keyword <' // &
                      TRIM(ComponentName("Netcdf Std Var",k)) // " Name" // ' not found, taking default')
              select case (k)
                case (1)
                  VarName='ex'
                case (2)
                  VarName='ey'
              end select
            ENDIF
          ENDIF
          NetCDFstatus = nf90_inq_varid(ncid,trim(VarName),varid)
          IF ( NetCDFstatus /= NF90_NOERR ) &
           CALL Fatal(Trim(SolverName),'Unable to get netcdf variable id for ' // TRIM(VarName))
          NetCDFstatus = nf90_get_var(ncid, Varid, StdObsArray(:, :,k), &
                            start = (/ XMinIndex, YMinIndex /),     &
                            count = (/ nx,ny/))
          IF ( NetCDFstatus /= NF90_NOERR ) &
            CALL Fatal(Trim(SolverName),'Unable to get netcdf variable' // TRIM(VarName))
        END DO
      END IF

      !! close netcdf
      NetCDFstatus = nf90_close(ncid)

      !!
      ALLOCATE(xobs(nx*ny,3),Vobs(nx*ny,VarDIM))
      IF (HaveSTD) ALLOCATE(StdObs(nx*ny,VarDIM))

      ALLOCATE(InElement(nx*ny))
      InElement(:)=-1
      xobs=0._dp
      Vobs=0._dp

      nobs=0
      Do i=1,nx
        Do j=1,ny
         IF (ANY(abs(ObsArray(i,j,:)-fillv(:)).LT.AEPS)) CYCLE
         nobs=nobs+1
         xobs(nobs,1)=xx(XMinIndex+i-1)
         xobs(nobs,2)=yy(YMinIndex+j-1)
         Do k=1,VarDIM
          Vobs(nobs,k)=ObsArray(i,j,k)
          IF (HaveSTD) & 
            StdObs(nobs,k)=max(StdObsArray(i,j,k),StdMin)
         END DO
        End do
      End do

      DEALLOCATE(xx,yy,ObsArray,fillv)
      IF (HaveSTD) DEALLOCATE(StdObsArray)

#else

    CALL FATAL(Trim(SolverName),'Elmer/Ice has not been installed with netcdf -convert your file to ASCII table--')

#endif

   Else !.not.NETCDFFormat

     ParallelFile = .False.
     ParallelFile = GetLogical(SolverParams,'Parallel Observation Files', Found)
     if (Parallel.AND.ParallelFile) THEN
       SideParFile = AddFilenameParSuffix(ObsFileName,'dat',Parallel,ParEnv % MyPe)
     ELSE
       SideParFile = trim(ObsFileName)
     ENDIF

     ProcessedFile = ListGetLogical(SolverParams,'Pre-Processed File')
     IF (ProcessedFile) THEN
        ActiveNumbering=ListGetLogical(SolverParams,"Element Number is Active Number")
     ENDIF
               
     open(IO,file=trim(SideParFile),status = 'old',iostat = ok)
     if(ok /= 0) then
       write(message,'(A,A)') 'Unable to open file ',TRIM(ObsFileName)
       CALL Fatal(Trim(SolverName),Trim(message))
     end if
     nobs=0
     do while(ok == 0)
       read(io,*,iostat = ok)
       if (ok == 0) nobs = nobs + 1
     end do
     close(IO)

     ALLOCATE(xobs(nobs,3),Vobs(nobs,VarDIM))
     IF (HaveSTD) ALLOCATE(StdObs(nobs,VarDIM))
     ALLOCATE(InElement(nobs))
     InElement(:)=-1

     Vobs=0.0_dp
     xobs=0.0_dp
     open(IO,file=trim(SideParFile))
     do i=1,nobs
       IF (ProcessedFile) THEN
         IF (HaveSTD) THEN
           read(IO,*) (xobs(i,j),j=1,CoordDIM),(Vobs(i,j),j=1,VarDIM),(StdObs(i,j),j=1,VarDIM),AI
         ELSE
           read(IO,*) (xobs(i,j),j=1,CoordDIM),(Vobs(i,j),j=1,VarDIM),AI
         END IF

         IF (ActiveNumbering) THEN
           Element => GetActiveElement(AI)
           IF (.NOT.ASSOCIATED(Element)) &
               CALL FATAL(SolverName,'No Active Element')
           InElement(i) = Element % ElementIndex
         ELSE
           InElement(i) = AI
         ENDIF
       ELSE
         IF (HaveSTD) THEN
           read(IO,*) (xobs(i,j),j=1,CoordDIM),(Vobs(i,j),j=1,VarDIM),(StdObs(i,j),j=1,VarDIM)
         ELSE
           read(IO,*) (xobs(i,j),j=1,CoordDIM),(Vobs(i,j),j=1,VarDIM)
         ENDIF
       ENDIF
     end do
     close(IO)
   ENDIF

   SAVE_USED_DATA = GetLogical(SolverParams,'Save used data', Found)
   IF (SAVE_USED_DATA) THEN
     SideFile =ListGetString(SolverParams,"output data file name",UnfoundFatal=.TRUE.)
     CALL SolverOutputDirectory( Solver, SideFile, OutputDirectory )
     SideFile = TRIM(OutputDirectory)// '/' //TRIM(SideFile)
     SideParFile = AddFilenameParSuffix(SideFile,'dat',Parallel,ParEnv % MyPe)
   END IF

  !!! End of First visit
    Firsttime=.false.
  Endif

! Get the variable solution
  Variable => VariableGet( Solver % Mesh % Variables, Trim(VariableName), UnFoundFatal=.TRUE.  )
  Values => Variable % Values
  Perm => Variable % Perm

! Get the sensitivity variable
  IF ( SEQL(TRIM(VariableName),'flow solution').OR.SEQL(TRIM(VariableName),'ssavelocity')) THEN
    FVarName = 'Velocityb'
  ELSE
    FVarName = TRIM(VariableName) // 'b'
  ENDIF
  VbSol => VariableGet( Solver % Mesh % Variables,TRIM(FVarName), UnFoundFatal=.TRUE.  )
  Vb => VbSol % Values
  VbPerm => VbSol % Perm
  Vb=0.0_dp

  IF (VbSol%DOFs.NE.Variable%DOFs) & 
        CALL FATAL(Trim(SolverName),'DOFs do not correspond')
  IF (VarDIM.GT.Variable%DOFs) & 
          CALL FATAL(Trim(SolverName),'Observed Variable DIM too large')
!
 Mesh => GetMesh()
  
! Temporary trick to disable Warnings as LocateParticleInMeshOctree
! will send a lot of warning for data points not found in the emsh
  IF (FirstRound) THEN
    WarnActive=OutputLevelMask(2)
    OutputLevelMask(2)=.FALSE.
  ENDIF

! Compute the cost
    Cost=0._dp
    rms=0._dp

    Element => GetActiveElement(1)
    ElementIndex = Element % ElementIndex

    CALL StartAdvanceOutput(SolverName,'Compute cost')
    Do s=1,nobs
     
     CALL AdvanceOutput(s,nobs)
     
     IF (FirstRound.AND.(InElement(s)<0)) then
       !Need to find in which Element the data point resides
       ! First check that it is within the computation domain
       Coord=0._dp
       Coord(1:CoordDIM)=xobs(s,1:CoordDIM)
       IF ((Coord(1).GT.xmax).OR.(Coord(1).LT.xmin)) CYCLE
       IF (CoordDIM.GE.2) THEN
        IF ((Coord(2).GT.ymax).OR.(Coord(2).LT.ymin)) CYCLE
       END IF

       IF (CoordDIM.LT.CoordinateSystemDimension()) THEN
       ! we are in a lower dimension than the mesh
       ! => use brute force search

       ! first check in the current elmement
        Element => Mesh % Elements(ElementIndex)
        n = GetElementNOFNodes(Element)
        NodeIndexes => Element % NodeIndexes
        ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
        IF (CoordDIM == 1) THEN
              ElementNodes % y(1:n) = 0.0_dp
              ElementNodes % z(1:n) = 0.0_dp
        ELSE IF (CoordDIM == 2) THEN
              ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
              ElementNodes % z(1:n) = 0.0_dp
        ENDIF
        IF (PointInElement(Element,ElementNodes,Coord(:),UVW)) THEN
          IF (.NOT.CheckPassiveElement(Element)) &
               InElement(s)=ElementIndex

        ELSE
         ! go through all elements
          Do i=1,Solver%NumberOfActiveElements
           Element => GetActiveElement(i)
           n = GetElementNOFNodes(Element)
           NodeIndexes => Element % NodeIndexes
           ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
           IF (CoordDIM == 1) THEN
              ElementNodes % y(1:n) = 0.0_dp
              ElementNodes % z(1:n) = 0.0_dp
           ELSE IF (CoordDIM == 2) THEN
              ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
              ElementNodes % z(1:n) = 0.0_dp
           ENDIF
           IF (PointInElement(Element,ElementNodes,Coord(:),UVW)) THEN
               IF (CheckPassiveElement(Element)) Exit
               ElementIndex= Element % ElementIndex
               InElement(s)=ElementIndex
               Exit
           END IF
          End Do
        END IF

      ELSE
      !> use particles meshoctree
       ElementIndex=-1  ! don't know why but if i don't reset ElmentIndex it fails
       CALL LocateParticleInMeshOctree( ElementIndex,Coord)
       If (ElementIndex.NE.0) THEN
         InElement(s)=ElementIndex  
         Element => Mesh % Elements(ElementIndex)
         IF (CheckPassiveElement(Element).OR. &
             (Element % PartIndex /= ParEnv % myPE)) THEN  
            InElement(s)=-1
         END IF
       END IF

      END IF
     ENDIF !End if FirstRound
      

    ! Data Point has been found in one element
      IF (InElement(s)>0) THEN
         Element => Mesh % Elements(InElement(s))
         n = GetElementNOFNodes(Element)
         NodeIndexes => Element % NodeIndexes
    ! set coords of highest occurring dimension to zero (to get correct path element)
          !-------------------------------------------------------------------------------
         ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
         IF (CoordDIM == 1) THEN 
            ElementNodes % y(1:n) = 0.0_dp
            ElementNodes % z(1:n) = 0.0_dp
         ELSE IF (CoordDIM == 2) THEN 
            ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
            ElementNodes % z(1:n) = 0.0_dp
         ELSE IF (CoordDIM == 3) THEN
            ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
            ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)
         END IF

         IF (.NOT.PointInElement(Element,ElementNodes,xobs(s,1:3),UVW))  THEN
              PRINT *,ParEnv%MyPE,'eindex',Element % ElementIndex
              PRINT *,ParEnv%MyPE,'obs. coords',xobs(s,1:2)
              PRINT *,ParEnv%MyPE,'node indexes',NodeIndexes(1:n)
              PRINT *,ParEnv%MyPE,'nodex',ElementNodes % x(1:n)
              PRINT *,ParEnv%MyPE,'nodey',ElementNodes % y(1:n)

              CALL FATAL(SolverName,'Point was supposed to be found in this element')
         ELSE
            stat = ElementInfo( Element,ElementNodes,UVW(1),UVW(2),UVW(3),SqrtElementMetric, &
                              Basis,dBasisdx )

            IF (CompNorm) THEN
              VmNorm=0._dp
              VoNorm=0._dp
              StdNorm=0._dp
              Do i=1,VarDIM
                um=SUM(Values(VDOFs*(Perm(NodeIndexes(1:n))-1)+i)*basis(1:n))
                VmNorm=VmNorm+um**2
                VoNorm=VoNorm+Vobs(s,i)**2
                IF (HaveSTD) &
                  StdNorm=StdNorm+StdObs(s,i)**2
              END DO
              misfit=sqrt(VmNorm)-sqrt(VoNorm)
              IF (HaveSTD) then
                rms=rms+misfit**2
                misfit=misfit/sqrt(StdNorm)
              ENDIF
              Cost = Cost + 0.5*Lambda*misfit**2
              VmNormb=Lambda*misfit
              IF (HaveSTD) &
                 VmNormb=VmNormb/sqrt(StdNorm)
              IF (VmNorm.GT.10*AEPS) THEN
                VmNormb=0.5*VmNormb*VmNorm**(-0.5)
              ELSE
                VmNormb=0._dp
              END IF

              Do i=1,VarDIM
                um=SUM(Values(VDOFs*(Perm(NodeIndexes(1:n))-1)+i)*basis(1:n))
                umb=2*um*VmNormb
                Do j=1,n
                  k=VDOFS*(VbPerm(NodeIndexes(j))-1)+i
                  Vb(k)=Vb(k)+umb*basis(j)
                End do
              END DO


            ELSE
              ! Update cost
              Do i=1,VarDIM
                V(i)=SUM(Values(VDOFs*(Perm(NodeIndexes(1:n))-1)+i)*basis(1:n))
                misfit=(V(i)-Vobs(s,i))
                IF (HaveSTD) then
                 rms=rms+misfit**2
                 misfit=misfit/StdObs(s,i)
                END IF
                Cost=Cost+0.5*Lambda*misfit**2
              End do

              !Derivative of Cost at nodal Point
              ! derivative of (V(i)-Vobs(s,i))**2/std**2=
              !  2*(V(i)-Vobs(s,i))/std**2 * (dV(i)/Vnodal)
              Do j=1,n
                Do i=1,VarDIM
                 k=VDOFS*(VbPerm(NodeIndexes(j))-1)+i
                 misfit=(V(i)-Vobs(s,i))
                 IF (HaveSTD) misfit=misfit/StdObs(s,i)**2
                 Vb(k)=Vb(k)+Lambda*misfit*basis(j)
                End do
              End Do

            ENDIF

          END IF

         ELSE

            WRITE(Message,'(a,I0,ES20.11E3,ES20.11E3,ES20.11E3,a)')&
                'Data Point ',s,xobs(s,1:3),' found in no element'
            CALL Info( SolverName, Message,level=15)
         END IF

    END DO !Do on s


    if (NTOT < 0) NTOT=COUNT(InElement(:)>0)
!! Save Cost as a function of time

    TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )

    IF (Parallel) THEN
           CALL MPI_ALLREDUCE(NTOT,NTOT_S,1,&
                  MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)
           CALL MPI_ALLREDUCE(Cost,Cost_S,1,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
          IF (HaveSTD) &
             CALL MPI_ALLREDUCE(rms,rms_S,1,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
          CostVar => VariableGet( Solver % Mesh % Variables, CostSolName )
          IF (ASSOCIATED(CostVar)) THEN
                 CostVar % Values(1)=Cost_S
          END IF
         IF (Solver % Matrix % ParMatrix % ParEnv % MyPE == 0) then
                 OPEN (IO, FILE=CostFile,POSITION='APPEND')
                 IF (HaveSTD) THEN
                   write(IO,'(4(ES20.11E3))') TimeVar % Values(1),Cost_S,& 
                           sqrt(2.0*Cost_S/(NTOT_S*Lambda)),sqrt(rms_S/NTOT_S)
                 ELSE
                   write(IO,'(3(ES20.11E3))') TimeVar % Values(1),Cost_S,sqrt(2.0*Cost_S/(NTOT_S*Lambda))
                 END IF
                 CLOSE(IO)
                 write(Message,'(A,A,I0)') trim(SolverName),'total number of data points:',NTOT_S
                 call INFO(SolverName,Message,level=3)
         End if
   ELSE
            CostVar => VariableGet( Solver % Mesh % Variables, CostSolName )
            IF (ASSOCIATED(CostVar)) THEN
                    CostVar % Values(1)=Cost
            END IF
            OPEN (IO, FILE=CostFile,POSITION='APPEND')
            IF (HaveSTD) THEN
              write(IO,'(4(ES20.11E3))') TimeVar % Values(1),Cost,&
                      sqrt(2.0*Cost/(NTOT*Lambda)),sqrt(rms/NTOT)
            ELSE 
              write(IO,'(3(ES20.11E3))') TimeVar % Values(1),Cost,sqrt(2.0*Cost/(NTOT*Lambda))
            END IF
            close(IO)
            write(Message,'(A,A,I0)') trim(SolverName),'total number of data points:',NTOT
            call INFO(SolverName,Message,level=3)
            Cost_S=Cost
   END IF

   Solver % Variable % Values(1)=Cost_S
   
   If (FirstRound.AND.SAVE_USED_DATA) THEN


     open(IO,File=SideParFile)
     Do s=1,nobs
        If (InElement(s)>0) then
           DO i=1,CoordDIM
             write(IO,'(ES20.11E3)',ADVANCE='NO') xobs(s,i)
           END DO
           DO i=1,VarDIM
             write(IO,'(ES20.11E3)',ADVANCE='NO') Vobs(s,i)
           END DO
           IF (HaveSTD) THEN
             DO i=1,VarDIM
               write(IO,'(ES20.11E3)',ADVANCE='NO') StdObs(s,i)
             END DO
           END IF
           write(IO,'(x,I0)') InElement(s)
        End if
      End do
     close(IO)

    END IF

! Reset Warning level to previous value
  IF (FirstRound) THEN
    OutputLevelMask(2)=WarnActive
  ENDIF

   FirstRound=.False.
   Return

 1000  format('#date,time,',a2,'/',a2,'/',a4,',',a2,':',a2,':',a2)
 1001  format('#lambda,',e15.8)

 CONTAINS
 ! Find the min and max indexes of values within bBox
 SUBROUTINE MinMaxIndex(x,n,minx,maxx,MinIndex,MaxIndex)
 IMPLICIT NONE
 REAL(KIND=dp),INTENT(IN) :: x(:),minx,maxx
 INTEGER,INTENT(IN) :: n
 INTEGER,INTENT(OUT) :: MinIndex,MaxIndex

 IF (n.EQ.1) THEN
         MinIndex=1
         MaxIndex=1
         Return
 ENDIF

 ! coordinates should be  monotonically increasing or
 ! decreasing 
 IF ((x(2)>x(1)).AND.(x(n)>x(n-1))) THEN
   MinIndex = MAXLOC(x, DIM=1,MASK=(x < minx))
   MinIndex=MAX(MinIndex,1)

   MaxIndex = MINLOC(x, DIM=1,MASK=(x > maxx))
   IF (MaxIndex.LE.1) MaxIndex=n
 ELSE IF ((x(2)<x(1)).AND.(x(n)<x(n-1))) THEN
   MaxIndex = MAXLOC(x, DIM=1,MASK=(x < minx))
   IF (MaxIndex.LE.1) MaxIndex=n

   MinIndex = MINLOC(x, DIM=1,MASK=(x > maxx))
   MinIndex=MAX(MinIndex,1)
 ELSE
   CALL FATAL(SolverName,'coordinate is neither monotonically increasing or decreasing')
 ENDIF

 IF (MaxIndex.LT.MinIndex) THEN
    WRITE(message,*) 'Error Min Max Indexes ',MinIndex,MaxIndex
    CALL WARN(SolverName,message)
    WRITE(message,*) 'Min values ',MINVAL(x),minx
    CALL WARN(SolverName,message)
    WRITE(message,*) 'Max values ',MAXVAL(x),maxx
    CALL WARN(SolverName,message)
    CALL FATAL(SolverName,'This is a fatal error')
 ENDIF

 END SUBROUTINE MinMaxIndex

!------------------------------------------------------------------------------
END SUBROUTINE Adjoint_CostDiscSolver
!------------------------------------------------------------------------------
! *****************************************************************************
