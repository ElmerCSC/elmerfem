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
! *  Authors: f. Gillet-Chaulet (LGGE, Grenoble,France)
! *  Email:   gillet-chaulet@lgge.obs.ujf-grenoble.fr
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! * 
! *****************************************************************************
SUBROUTINE AdjointSSA_CostDiscSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!Compute the Cost function of the SSA Adjoint inverse Problem  as: SUM_1^Nobs (u-u^obs)2
!            where u is evaluated at observation location
!   Serial/Parallel    2D/3D
!
!     OUTPUT are : J and DJDu (==Velocityb variable used as forcing of the SSA adjoint problem)
!
!    !! Be careful this solver will reset the cost and DJDu to 0; so it has to
!    be used as the first cost solver if regularistaion of flux divergence cost
!    solvers are used in the simulation!!
!
!
!     INPUT PARAMETERS are:
!
!      In solver section:
!               Problem Dimension = Integer (default = Coordinate system dimension)
!               Cost Filename = File (default = 'CostOfT.dat')
!               Cost Variable Name = String (default= 'CostValue')
!               Lambda = Real (default= '1.0')
!               Observed Variable Name = String 
!               Observation File Name = File
!               Parallel observation Files = Logical
!               Save use data = logical
!
!      Variables
!                Velocityb (forcing for the adjoint pb;DOFs== Pb dimension)
!
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
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: DefaultCostFile = 'CostOfT.dat'
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName,CostFile,UsedDataFile
  CHARACTER(LEN=MAX_NAME_LEN) :: CostSolName,VariableName,ObsFileName
  TYPE(Element_t),POINTER ::  Element
  TYPE(Variable_t), POINTER :: TimeVar,CostVar
  TYPE(Variable_t), POINTER :: VelocitybSol,Variable
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Nodes_t) :: ElementNodes
  REAL(KIND=dp), POINTER :: Vb(:),Values(:)
  INTEGER, POINTER :: NodeIndexes(:)
  INTEGER, POINTER :: VbPerm(:),Perm(:)
  Logical :: Firsttime=.true.,Found,Parallel,ParallelFile,stat,Gotit
  integer :: i,j,k,l,s,t,n,NMAX,DIM,ierr,c,ok
  real(kind=dp) :: Cost,Cost_S,Lambda
  real(kind=dp) :: Coord(3),UVW(3),coeff,SqrtElementMetric,x
  REAL(KIND=dp) :: Basis(Model % MaxElementNodes), dBasisdx(Model % MaxElementNodes,3)
  REAL(KIND=dp) :: xmin,ymin,xmax,ymax
  INTEGER,SAVE :: VDOFs
  REAL(KIND=dp),ALLOCATABLE,SAVE :: V(:)
  REAL(KIND=dp),ALLOCATABLE,SAVE :: xobs(:,:),Vobs(:,:)
  REAL(KIND=dp),ALLOCATABLE :: xx(:),yy(:),ux(:,:),uy(:,:)
  INTEGER,ALLOCATABLE,SAVE :: InElement(:)
  INTEGER :: ElementIndex
  INTEGER :: XMinIndex,XMaxIndex,YMinIndex,YMaxIndex,nx,ny
  INTEGER,SAVE :: NTOT=-1,NTOT_S
  integer,SAVE :: nobs
  LOGICAL,SAVE :: FirstRound=.True.,SAVE_USED_DATA=.False.
  CHARACTER*10 :: date,temps
  INTEGER,PARAMETER :: IO=12
  LOGICAL :: WarnActive

  LOGICAL :: NETCDFFormat
  INTEGER :: NetcdfStatus,varid,ncid
  CHARACTER(LEN=MAX_NAME_LEN) :: Xdim,Ydim,VarName
  REAL(KIND=dp):: UxFillv,UyFillv
  REAL(KIND=dp):: NodalU,NodalV

  save Firsttime,Parallel 
  save SolverName,CostSolName,CostFile,VariableName
  save ElementNodes


  SolverParams => GetSolverParams()
  DIM=GetInteger(SolverParams ,'Problem Dimension',Found)
  If (.NOT.Found) then
     CALL WARN(SolverName,'Keyword >Problem Dimension< not found, assume DIM = CoordinateSystemDimension()')
     DIM = CoordinateSystemDimension()
  Endif

! Get min/max coordinates of current mesh
  xmin=MINVAL(Solver%Mesh%Nodes%x(:))
  xmax=MAXVAL(Solver%Mesh%Nodes%x(:))
  IF (DIM.EQ.2) THEN
    ymin=MINVAL(Solver%Mesh%Nodes%y(:))
    ymax=MAXVAL(Solver%Mesh%Nodes%y(:))
  ENDIF

  Lambda =  GetConstReal( SolverParams,'Lambda', Found)
  IF(.NOT.Found) THEN
      CALL WARN(SolverName,'Keyword >Lambda< not found  in section >Solver<')
      CALL WARN(SolverName,'Taking default value Lambda=1.0')
      Lambda = 1.0_dp
      CALL ListAddConstReal( SolverParams,'Lambda',Lambda)
  End if

  If (Firsttime) then
    N = model % MaxElementNodes
    allocate(ElementNodes % x(N), ElementNodes % y(N), ElementNodes % z(N))

!!!!!!! Check for parallel run 
    Parallel = .FALSE.
    IF ( ASSOCIATED( Solver % Matrix % ParMatrix ) ) THEN
            IF ( Solver %  Matrix % ParMatrix % ParEnv % PEs > 1 )  THEN
                    Parallel = .TRUE.
            END IF
    END IF

     WRITE(SolverName, '(A)') 'CostSolver_Adjoint'

!!!!!!!!!!! get Solver Variables

  CostFile = ListGetString(Solver % Values,'Cost Filename',Found )
  IF (.NOT. Found) CostFile = DefaultCostFile
  CALL DATE_AND_TIME(date,temps)
  If (Parallel) then
     if (ParEnv % MyPe.EQ.0) then
        OPEN (IO, FILE=CostFile)
             write(IO,1000) date(5:6),date(7:8),date(1:4),temps(1:2),temps(3:4),temps(5:6)
             write(IO,1001) Lambda
             write(IO,'(A)') '# iter, Lambda*J0, rms'
        CLOSE(IO)
     End if
  Else
       OPEN (IO, FILE=CostFile)
             write(IO,1000) date(5:6),date(7:8),date(1:4),temps(1:2),temps(3:4),temps(5:6)
             write(IO,1001) Lambda
             write(IO,*)'# iter, Lambda*J0, rms'
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
   ALLOCATE(V(VDOFs))

!! Get the obs
   ObsFileName =  ListGetString( SolverParams,'Observation File Name', UnFoundFatal=.TRUE.)

!! I the file netcf?
   k = INDEX(ObsFileName ,'.nc' )
   NETCDFFormat = ( k /= 0 )

#ifdef HAVE_NETCDF
   IF (NETCDFFormat) then
      IF (VDOFs.NE.2) CALL FATAL(Trim(SolverName),'Netcdf only supported for 2D')

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


      !! get the index of values that are within the mesh bounding box
      XMinIndex = MINLOC(xx, DIM=1,MASK=(xx >= xmin))  
      XMaxIndex = MAXLOC(xx, DIM=1,MASK=(xx <= xmax))
      IF (XMinIndex.GT.XMaxIndex) &
         CALL FATAL(Trim(SolverName),'x values in decreasing order -- not supported --')
      YMinIndex = MINLOC(yy, DIM=1,MASK=(yy >= ymin))  
      YMaxIndex = MAXLOC(yy, DIM=1,MASK=(yy <= ymax))
      IF (YMinIndex.GT.YMaxIndex) &
         CALL FATAL(Trim(SolverName),'y values in decreasing order -- not supported --')
 

     !! get ux and uy
      nx=XMaxIndex-XMinIndex+1
      ny=YMaxIndex-YMinIndex+1
      ALLOCATE(ux(nx,ny),uy(nx,ny))

      VarName=ListGetString( SolverParams, "Ux Var Name", Found )
      if (.NOT.Found) VarName='vx'
      NetCDFstatus = nf90_inq_varid(ncid,trim(VarName),varid)
      IF ( NetCDFstatus /= NF90_NOERR ) &
           CALL Fatal(Trim(SolverName),'Unable to get netcdf ux variable id')
      NetCDFstatus = nf90_get_var(ncid, Varid, ux(:, :), &
                            start = (/ XMinIndex, YMinIndex /),     &
                            count = (/ nx,ny/))
      IF ( NetCDFstatus /= NF90_NOERR ) &
           CALL Fatal(Trim(SolverName),'Unable to get netcdf ux variable')
      NetCDFstatus = nf90_get_att(ncid, varid,"_FillValue",Uxfillv)
      IF ( NetCDFstatus /= NF90_NOERR ) &
           CALL FATAL(Trim(SolverName),'Ux has no attribute _FillValue')

      VarName=ListGetString( SolverParams, "Uy Var Name", Found )
      if (.NOT.Found) VarName='vy'
      NetCDFstatus = nf90_inq_varid(ncid,trim(VarName),varid)
      IF ( NetCDFstatus /= NF90_NOERR ) &
           CALL Fatal(Trim(SolverName),'Unable to get netcdf uy variable id')
      NetCDFstatus = nf90_get_var(ncid, Varid, uy(:, :), &
                            start = (/ XMinIndex, YMinIndex /),     &
                            count = (/ nx,ny /))
      IF ( NetCDFstatus /= NF90_NOERR ) &
           CALL Fatal(Trim(SolverName),'Unable to get netcdf uy variable')
      NetCDFstatus = nf90_get_att(ncid, varid,"_FillValue",Uyfillv)
      IF ( NetCDFstatus /= NF90_NOERR ) &
           CALL FATAL(Trim(SolverName),'Uy has no attribute _FillValue')

      !! close netcdf
      NetCDFstatus = nf90_close(ncid)

      !!
      ALLOCATE(xobs(nx*ny,3),Vobs(nx*ny,VDOFs))
      xobs=0._dp
      Vobs=0._dp

      nobs=0
      Do i=1,nx
        Do j=1,ny
          NodalU=ux(i,j)
          NodalV=uy(i,j)
          IF ((NodalU.EQ.Uxfillv).OR.(NodalV.EQ.Uyfillv)) CYCLE
          nobs=nobs+1
          xobs(nobs,1)=xx(XMinIndex+i-1)
          xobs(nobs,2)=yy(YMinIndex+j-1)
          Vobs(nobs,1)=NodalU
          Vobs(nobs,2)=NodalV
        End do
      End do

      DEALLOCATE(xx,yy,ux,uy)
   Endif
#else
   IF (NETCDFFormat) &
      CALL FATAL(Trim(SolverName),'Elmer/Ice has not been installed with netcdf -convert your file to ASCII table--')

   ParallelFile = .False.
   ParallelFile = GetLogical(SolverParams,'Parallel Observation Files', Found)
   if (Parallel.AND.ParallelFile) &
    write(ObsFileName,'(A,A,I0)') trim(ObsFileName),'.',ParEnv % MyPE  
               
   open(IO,file=trim(ObsFileName),status = 'old',iostat = ok)
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


   ALLOCATE(xobs(nobs,3),Vobs(nobs,VDOFs))

   Vobs=0.0_dp
   xobs=0.0_dp
   open(IO,file=trim(ObsFileName))
   do i=1,nobs
     read(IO,*) (xobs(i,j),j=1,DIM),(Vobs(i,j),j=1,VDOFs)
   end do
   close(IO)
#endif

  ALLOCATE(InElement(nobs))
  InElement(:)=-1

   SAVE_USED_DATA = GetLogical(SolverParams,'Save used data', Found)
  !!! End of First visit
    Firsttime=.false.
  Endif

! Get the variable solution
  Variable => VariableGet( Solver % Mesh % Variables, Trim(VariableName), UnFoundFatal=.TRUE.  )
  Values => Variable % Values
  Perm => Variable % Perm

! Get the adjoint
  VelocitybSol => VariableGet( Solver % Mesh % Variables, 'Velocityb', UnFoundFatal=.TRUE.  )
  Vb => VelocitybSol % Values
  VbPerm => VelocitybSol % Perm
  Vb=0.0_dp

  IF (VelocitybSol%DOFs.NE.Variable%DOFs) & 
        CALL FATAL(Trim(SolverName),'DOFs do not correspond')

! Temporary trick to disable Warnings as LocateParticleInMeshOctree
! will send a lot of warning for data points not found in the emsh
  IF (FirstRound) THEN
    WarnActive=OutputLevelMask(2)
    OutputLevelMask(2)=.FALSE.
  ENDIF

! Compute the cost
    Cost=0._dp

    CALL StartAdvanceOutput(SolverName,'Compute cost')
    Do s=1,nobs

     CALL AdvanceOutput(s,nobs)
     
     IF (FirstRound) then
      !Need to find in which Element the data point resides
      ElementIndex=-1  ! don't know why but if i don't reset ElmentIndex it fails
      Coord=0._dp
      Coord(1:DIM)=xobs(s,1:DIM)
      IF ((Coord(1).GT.xmax).OR.(Coord(1).LT.xmin)) CYCLE
      IF (DIM.EQ.2) THEN
        IF ((Coord(2).GT.ymax).OR.(Coord(2).LT.ymin)) CYCLE
      END IF
      CALL LocateParticleInMeshOctree( ElementIndex,Coord)
      If (ElementIndex.NE.0) THEN
         InElement(s)=ElementIndex
         Element => GetActiveElement(InElement(s))
         IF (CheckPassiveElement(Element)) THEN
            InElement(s)=-1
         END IF
      END IF
     ENDIF !End if FirstRound
      

    ! Data Point has been found in one element
      IF (InElement(s)>0) THEN
         Element => GetActiveElement(InElement(s))
         n = GetElementNOFNodes()
         NodeIndexes => Element % NodeIndexes
    ! set coords of highest occurring dimension to zero (to get correct path element)
          !-------------------------------------------------------------------------------
         ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
         IF (DIM == 1) THEN !1D SSA
            ElementNodes % y(1:n) = 0.0_dp
            ElementNodes % z(1:n) = 0.0_dp
         ELSE IF (DIM == 2) THEN !2D SSA
            ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
            ElementNodes % z(1:n) = 0.0_dp
         ELSE
            WRITE(Message,'(a,i1,a)')&
                'It is not possible to compute SSA problems with DOFs=',&
                DIM, ' . Aborting'
            CALL Fatal( SolverName, Message)
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

            ! Variable at obs point
            Do i=1,VDOFs
              V(i)=SUM(Values(VDOFs*(Perm(NodeIndexes(1:n))-1)+i)*basis(1:n))
            End do

            ! Update cost
            Do i=1,VDOFs
              Cost=Cost+0.5*Lambda*(V(i)-Vobs(s,i))*(V(i)-Vobs(s,i))
            End do
            !PRINT *,V,Vobs

            !Derivative of Cost at nodal Point
            Do j=1,n
              Do i=1,VDOFs
               k=VDOFS*(VbPerm(NodeIndexes(j))-1)+i
               Vb(k)=Vb(k)+Lambda*(V(i)-Vobs(s,i))*basis(j)
              End do
            End Do

          END IF

         ELSE

            WRITE(Message,'(a,I0,a)')&
                'Data Point',s,'found in no element'
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
          CostVar => VariableGet( Solver % Mesh % Variables, CostSolName )
          IF (ASSOCIATED(CostVar)) THEN
                 CostVar % Values(1)=Cost_S
          END IF
         IF (Solver % Matrix % ParMatrix % ParEnv % MyPE == 0) then
                 OPEN (IO, FILE=CostFile,POSITION='APPEND')
                 write(IO,'(e13.5,2x,e15.8,2x,e15.8)') TimeVar % Values(1),Cost_S,sqrt(2.0*Cost_S/(NTOT_S*Lambda))
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
            write(IO,'(e13.5,2x,e15.8,2x,e15.8)') TimeVar % Values(1),Cost,sqrt(2.0*Cost/(NTOT*Lambda))
            close(IO)
            write(Message,'(A,A,I0)') trim(SolverName),'total number of data points:',NTOT
            call INFO(SolverName,Message,level=3)
   END IF
   
   If (FirstRound) THEN
     IF (SAVE_USED_DATA) THEN
      If (Parallel) then
       write(UsedDataFile,'(A,A,I0)') trim(SolverName),'.useddata.',ParEnv % MyPE
      Else
       write(UsedDataFile,'(A,A)') trim(SolverName),'.useddata'
      End if

       open(IO,File=UsedDataFile)
       Do s=1,nobs
         If (InElement(s)>0) then
            if ((DIM.eq.1).AND.(VDOFS.EQ.1)) Then
               write(IO,'(e13.5,2x,e15.8)') (xobs(s,i),i=1,DIM),(Vobs(s,i),i=1,VDOFs)
            Else if ((DIM.eq.2).AND.(VDOFS.EQ.2)) Then
               write(IO,'(e15.8,2x,e15.8,2x,e15.8,2x,e15.8)') (xobs(s,i),i=1,DIM),(Vobs(s,i),i=1,VDOFs)
            End if
         End if
       End do
       close(IO)
     END IF
   End if

! Reset Warning level to previous value
  IF (FirstRound) THEN
    OutputLevelMask(2)=WarnActive
  ENDIF

   FirstRound=.False.
   Return

 1000  format('#date,time,',a1,'/',a1,'/',a4,',',a2,':',a2,':',a2)
 1001  format('#lambda,',e15.8)

!------------------------------------------------------------------------------
END SUBROUTINE AdjointSSA_CostDiscSolver
!------------------------------------------------------------------------------
! *****************************************************************************
