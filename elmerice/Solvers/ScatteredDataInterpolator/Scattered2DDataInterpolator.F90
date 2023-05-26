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
! *  Authors: 
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! * 
! *****************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Interpolate scattered 2D data read from an ASCII file (x y Value)
!    in the mesh nodes using:
!     - nn-c library (http://code.google.com/p/nn-c/): linear and Natural
!               Neighbours interpolation
!     - csa-c library (http://code.google.com/p/csa-c): cubic spline approximation
!
! TODO: - add possibility to prescibe std for csa method
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE Scattered2DDataInterpolator( Model,Solver,dt,TransientSimulation )

        USE Netcdf
        USE DefUtils
        USE NearestNeighbour

        IMPLICIT NONE

        TYPE(Solver_t), TARGET :: Solver
        TYPE(Model_t) :: Model
        REAL(KIND=dp) :: dt
        LOGICAL :: TransientSimulation

        TYPE(ValueList_t), POINTER :: Params
        TYPE(Variable_t), POINTER :: Var
        TYPE(Element_t), POINTER :: Element
        REAL(KIND=dp), POINTER :: Values(:)
        INTEGER, POINTER :: Perm(:)

        REAL(KIND=DP) :: Win,NaNVal,fillv
        REAL(KIND=DP) :: Val,MinMaxVals(2)
        REAL(KIND=DP) :: x,y,z,MinD,D
        REAL(KIND=DP) :: xmin,xmax,ymin,ymax,BBoxdx
        REAL(KIND=DP),allocatable :: xx(:),yy(:),DEM(:,:)

        INTEGER,parameter :: io=20
        INTEGER :: ok,nNaN
        INTEGER :: i,j,k,t,kmin,NoVar
        INTEGER :: nppcin,csakin
        INTEGER :: NetcdfStatus,varid,ncid
        INTEGER :: compt,nx,ny,nval
        INTEGER :: nodeind,Varind
        INTEGER :: nlines
        INTEGER :: XMinIndex,XMaxIndex
        INTEGER :: YMinIndex,YMaxIndex

        CHARACTER(LEN=MAX_STRING_LEN) :: DataF
        CHARACTER(LEN=MAX_NAME_LEN) :: TargetVariableName,VariableName
        CHARACTER(LEN=MAX_NAME_LEN) :: Name,FName,WName,MName,CSVName,tmpName
        CHARACTER(LEN=MAX_NAME_LEN),parameter :: &
                         SolverName='Scattered2DDataInterpolator'
        CHARACTER(LEN=MAX_NAME_LEN) :: Xdim,dimName,FillName
        CHARACTER(2) :: method

        LOGICAL :: GotVar,GotTVar,Found,UnFoundFatal=.TRUE.
        LOGICAL :: HaveMin,HaveMax
        LOGICAL :: INBOX
        LOGICAL :: CheckBBox,CheckNaN,ReplaceNaN,NETCDFFormat
        LOGICAL :: GoodVal,HaveFillv
        LOGICAL,dimension(:,:), allocatable :: mask
        LOGICAL,dimension(:), allocatable :: BBox
        LOGICAL :: Debug=.False.

        ! Variables to pass to the nn C library
        type(POINT),dimension(:),allocatable :: pout,pin,ptmp
        INTEGER(C_INT) :: nout,nin,nppc,csak
        REAL(C_DOUBLE) :: W

       CALL INFO(SolverName, &
           '-----Initialise Variables using nn Library----------',Level=5)

       Params => GetSolverParams()
       
       CheckBBox=.False.
       BBoxdx = 0._dp
       BBoxdx = ListGetConstReal(Params, 'Bounding Box dx', Found)
       If (Found) CheckBBox=.True.
 
       xmin=MINVAL(Model%Mesh%Nodes%x)-BBoxdx
       xmax=MAXVAL(Model%Mesh%Nodes%x)+BBoxdx
       ymin=MINVAL(Model%Mesh%Nodes%y)-BBoxdx
       ymax=MAXVAL(Model%Mesh%Nodes%y)+BBoxdx

       CheckNaN = ListGetLogical(Params, 'Look For NaN',Found)
       If(.NOT.Found) CheckNaN=.True.
       NaNVal = ListGetConstReal( Params, 'Replace NaN by', ReplaceNaN )

       ! Mesh nodes coordinates for the NN iterpolation
       nout=Model % Mesh % NumberOfNodes
       allocate(pout(nout))
   
       IF (DEBUG) open(10,file='MeshNodes.dat') !tmp

        Do i=1,Model % Mesh % NumberOfNodes
           pout(i)%x = Model % Mesh % Nodes % x(i)
           pout(i)%y = Model % Mesh % Nodes % y(i)
           IF (DEBUG) write(10,*) pout(i)%x,pout(i)%y
        End Do

        IF (DEBUG) close(10) !tmp

       ! Read variable to initialize and Data
        NoVar=0
        GotVar=.True.

        DO WHILE(GotVar)
            NoVar = NoVar + 1
            WRITE(Name,'(A,I0)') 'Variable ',NoVar

            VariableName = ListGetString( Params, TRIM(Name), GotVar )
            IF (.NOT.GotVar) exit
            
            WRITE(Name,'(A,I0)') 'Target Variable ',NoVar
            TargetVariableName=ListGetString( Params, TRIM(Name), GotTVar)
            IF (.NOT.GotTVar) TargetVariableName=VariableName

            Var => VariableGet(Model %  Mesh % Variables, TargetVariableName )
            IF(.NOT.ASSOCIATED(Var)) Then
               CALL VariableAddVector(Model % Mesh % Variables,Model % Mesh,Solver,TargetVariableName,1)
               Var => VariableGet(Model %  Mesh % Variables, TargetVariableName )
            ENDIF
            Values => Var % Values
            Perm => Var % Perm

            WRITE (FName,'(A,I0,A)') 'Variable ',NoVar,' Data File'
            DataF = ListGetString( Params, TRIM(FName), Found, UnFoundFatal )

            k = INDEX( DataF,'.nc' )
            NETCDFFormat = ( k /= 0 )

            IF (NETCDFFormat) then
              CALL INFO(SolverName,'Data File is in netcdf format', Level=5)
            Else
              CALL INFO(SolverName,'Data File is in ascii', Level=5)
            Endif
            

            WRITE (WName,'(A,I0,A)') 'Variable ',NoVar,' W'
            Win = ListGetConstReal( Params, TRIM(WName), Found )
            if (.NOT.Found) then
                W=-HUGE(0.0d0)
            else
                W=Win
            endif

            WRITE (MName,'(A,I0,A)') 'Variable ',NoVar,' method'
            method = ListGetString( Params, TRIM(MName), Found )
            if (.NOT.Found) then
                method='n'
            endif

            WRITE (MName,'(A,I0,A)') 'Variable ',NoVar,&
                                     ' Valid Min Value'
            Val = ListGetCReal( Params, TRIM(MName), Found )
            IF (Found) THEN
              MinMaxVals(1)=Val
            ELSE
              MinMaxVals(1)=-HUGE(MinMaxVals(1))
            ENDIF
            HaveMin=Found
            WRITE (MName,'(A,I0,A)') 'Variable ',NoVar,&
                                     ' Valid Max Value'
            Val = ListGetCReal( Params, TRIM(MName), Found )
            IF (Found) THEN
              MinMaxVals(2)=Val
            ELSE
              MinMaxVals(2)=HUGE(MinMaxVals(1))
            ENDIF
            HaveMax=Found
            
           IF (.NOT.NETCDFFormat) Then
               open(unit = io, file = TRIM(DataF), status = 'old',iostat = ok)

               if(ok /= 0) then
                  write(message,'(A,A)') 'Unable to open file ',TRIM(DataF)
                  CALL Fatal(SolverName,Trim(message))
               end if
            
               nin=0
               nlines=0
               !count the line number in the file
               do while(ok == 0) 
                 read(io,*,iostat = ok) x,y,z
                 if (ok == 0) THEN
                   nlines = nlines + 1
                   IF (CheckBBox) THEN
                    INBOX=(x.GE.xmin).AND.(x.LE.xmax).AND.&
                          (y.GE.ymin).AND.(y.LE.ymax)
                    IF (.NOT.INBOX) CYCLE
                   ENDIF
                   nin = nin + 1
                  ENDIF
               end do 

               IF (nin.eq.0) then
                 WRITE(message,'(A,A)') 'No Data within BBox found in',TRIM(DataF)
                 CALL Fatal(SolverName,TRIM(message))
               ENDIF
                                        
          
               allocate(pin(nin))

               ! comes back to beginning of file
               rewind(unit=io,iostat=ok)

               ! read data
               compt=0
               do i = 1, nlines
                   read(io,*,iostat = ok) x,y,z
                   IF (CheckBBox) THEN
                    INBOX=(x.GE.xmin).AND.(x.LE.xmax).AND.&
                          (y.GE.ymin).AND.(y.LE.ymax)
                    IF (.NOT.INBOX) CYCLE
                   ENDIF
                   compt = compt + 1
                   pin(compt)%x=x
                   pin(compt)%y=y
                   pin(compt)%z=z
                   IF (HaveMin) pin(compt)%z=Max(MinMaxVals(1),pin(compt)%z)
                   IF (HAVEMax) pin(compt)%z=Min(MinMaxVals(2),pin(compt)%z)
               end do
               IF (compt /= nin) &
                   CALL Fatal(SolverName,'error in reading values')
               close(io)
           ELSE
               NetCDFstatus = NF90_OPEN(trim(DataF),NF90_NOWRITE,ncid)
               IF ( NetCDFstatus /= NF90_NOERR ) THEN
                 CALL Fatal(SolverName,'Unable to open NETCDF File')
               END IF
               WRITE (dimName,'(A,I0,A)') 'Variable ',&
                         NoVar,' x-dim Name'
               Xdim=ListGetString( Params, TRIM(dimName), Found )
               if (.NOT.Found) Xdim='x'
               NetCDFstatus = nf90_inq_dimid(ncid, trim(Xdim) , varid)
               IF ( NetCDFstatus /= NF90_NOERR ) THEN
                 CALL Fatal(SolverName,'Unable to  get netcdf x-dim Id')
               ENDIF
               NetCDFstatus = nf90_inquire_dimension(ncid,varid,len=nx)
               IF ( NetCDFstatus /= NF90_NOERR ) THEN
                 CALL Fatal(SolverName,'Unable to  get netcdf nx')
               ENDIF
               WRITE (dimName,'(A,I0,A)') 'Variable ',&
                         NoVar,' y-dim Name'
               Xdim=ListGetString( Params, TRIM(dimName), Found )
               if (.NOT.Found) Xdim='y'
               NetCDFstatus = nf90_inq_dimid(ncid, trim(Xdim) , varid)
               IF ( NetCDFstatus /= NF90_NOERR ) THEN
                 CALL Fatal(SolverName,'Unable to  get netcdf y-dim Id')
               ENDIF
               NetCDFstatus = nf90_inquire_dimension(ncid,varid,len=ny)
               IF ( NetCDFstatus /= NF90_NOERR ) THEN
                 CALL Fatal(SolverName,'Unable to  get netcdf ny')
               ENDIF

               !! allocate good size
               allocate(xx(nx),yy(ny))

               !! Get X variable
               WRITE (dimName,'(A,I0,A)') 'Variable ',&
                       NoVar,' x-Var Name'
               Xdim=ListGetString( Params, TRIM(dimName), Found )
               if (.NOT.Found) Xdim='x'
               NetCDFstatus = nf90_inq_varid(ncid,trim(Xdim),varid)
               IF ( NetCDFstatus /= NF90_NOERR ) THEN
                 CALL Fatal(SolverName, 'Unable to get netcdf x-variable id')
               ENDIF
               NetCDFstatus = nf90_get_var(ncid, varid,xx)
               IF ( NetCDFstatus /= NF90_NOERR ) THEN
                 CALL Fatal(SolverName, 'Unable to get netcdf x-variable ')
               ENDIF
               !! Get Y variable
               WRITE (dimName,'(A,I0,A)') 'Variable ',&
                       NoVar,' y-Var Name'
               Xdim=ListGetString( Params, TRIM(dimName), Found )
               if (.NOT.Found) Xdim='y'
               NetCDFstatus = nf90_inq_varid(ncid,trim(Xdim),varid)
               IF ( NetCDFstatus /= NF90_NOERR ) THEN
                 CALL Fatal(SolverName,'Unable to get netcdf y-variable id')
               ENDIF
               NetCDFstatus = nf90_get_var(ncid, varid,yy)
               IF ( NetCDFstatus /= NF90_NOERR ) THEN
                 CALL Fatal(SolverName,'Unable to get netcdf y-variable')
               ENDIF

               !! Check that there is data within the domain
               IF ((MAXVAL(xx).LT.xmin).OR.(MINVAL(xx).GT.xmax)&
                .OR.(MAXVAL(yy).LT.ymin).OR.(MINVAL(yy).GT.ymax)) &
                 CALL Fatal(SolverName,'No data within model domain')

               !! only get Vars within BBox
               IF (CheckBBox) THEN
                 CALL MinMaxIndex(xx,nx,xmin,xmax,XMinIndex,XMaxIndex)
                 CALL MinMaxIndex(yy,ny,ymin,ymax,YMinIndex,YMaxIndex)
               ELSE
                 XMinIndex = 1
                 XMaxIndex = nx
                 YMinIndex = 1
                 YMaxIndex = ny
               END IF
               nx=XMaxIndex-XMinIndex+1
               ny=YMaxIndex-YMinIndex+1

               write(message,'(A,I0,A,I0,A)') 'NETCDF: reading nx=',nx,&
                     ' and ny=',ny,' data points'
              CALL INFO(SolverName,Trim(message),Level=5)
               write(message,*) 'X Indexes: ',&
                XMinIndex,XMaxIndex,xx(XMinIndex),xx(XMaxIndex)
              CALL INFO(SolverName,Trim(message),Level=10)
               write(message,*) 'Y Indexes: ', &
                YMinIndex,YMaxIndex,yy(YMinIndex),yy(YMaxIndex)
              CALL INFO(SolverName,Trim(message),Level=10)

               allocate(DEM(nx,ny),mask(nx,ny))

               !! Get the variable
               NetCDFstatus = nf90_inq_varid(ncid,TRIM(VariableName),varid)
               IF ( NetCDFstatus /= NF90_NOERR ) THEN
                   CALL Fatal(SolverName,'Unable to get netcdf variable id')
               ENDIF
               NetCDFstatus = nf90_get_var(ncid, varid,DEM(:,:),&
                       start = (/ XMinIndex, YMinIndex /),     &
                       count = (/ nx,ny/))
               IF ( NetCDFstatus /= NF90_NOERR ) THEN
                   CALL Fatal(SolverName, 'Unable to get netcdf variable')
               ENDIF
               HaveFillV=.True.
               NetCDFstatus = nf90_get_att(ncid, varid,"_FillValue",fillv)
               IF ( NetCDFstatus /= NF90_NOERR ) THEN
                   HaveFillV=.False.
               ENDIF
               WRITE (FillName,'(A,I0,A)') 'Variable ',NoVar,' Fill Value'
               Val=ListGetConstReal(Params, TRIM(FillName) , Found)
               IF (Found) THEN
                HaveFillV=.TRUE.
                fillv=Val
               ENDIF


               !! Close NETCDF
               NetCDFstatus = nf90_close(ncid)

               !! Fill Data values
               IF (HaveFillV) then
                  mask(:,:)=(DEM(:,:).NE.fillv)
                  nin=COUNT(mask)
               Else
                  nin=nx*ny
               Endif

               write(message,'(A,I0,A,I0)') 'NETCDF: I Found ',nin,&
                     ' data points with non-missing  value over ',nx*ny
              CALL INFO(SolverName,Trim(message),Level=5)

               if (nin == 0) CALL Fatal(SolverName, &
                                  'no data with non-missing value??')

               allocate(pin(nin))
                
               IF (DEBUG) write(tmpName,'(A,A)') &
                          TRIM(VariableName),'Data.dat'
               IF (DEBUG) open(10,file=trim(tmpName)) !tmp

               compt=0
               Do i=1,nx
                  Do j=1,ny
                     GoodVal=.True.
                     IF (HaveFillV) GoodVal=(DEM(i,j).NE.fillv)
                     if (GoodVal) then
                       compt=compt+1
                       IF (compt.GT.nin) then
                          CALL Fatal(SolverName,&
                               'get more Non-Nan values than expected')
                       ENDIF
                       pin(compt)%x = xx(XMinIndex+i-1)
                       pin(compt)%y = yy(YMinIndex+j-1)
                       pin(compt)%z = DEM(i,j)
                       IF (HaveMin) &
                         pin(compt)%z=Max(MinMaxVals(1),pin(compt)%z)
                       IF (HAVEMax) &
                         pin(compt)%z=Min(MinMaxVals(2),pin(compt)%z)
                       IF (DEBUG) write(10,*) pin(compt)%x,pin(compt)%y,pin(compt)%z
                     endif         
                  End do
               End do
               IF (DEBUG) close(10) !tmp
               IF (compt /= nin) CALL Fatal(SolverName,&
                       'sorry I didn t found the good number of values')
               deallocate(xx,yy,DEM,mask)
           ENDIF

            ! call the nn C library
            SELECT CASE (method(1:1))
                  CASE ('c')
                   WRITE (CSVName,'(A,I0,A)') 'Variable ',NoVar,' nppc'
                   nppcin = ListGetInteger( Params, TRIM(CSVName), Found )
                   if (.NOT.Found) then  
                           nppc=-1
                   else
                           nppc=nppcin
                   endif
                   WRITE (CSVName,'(A,I0,A)') 'Variable ',NoVar,' k'
                   csakin = ListGetInteger( Params, TRIM(CSVName), Found )
                   if (.NOT.Found) then  
                           csak=-1
                   else
                           csak=csakin
                   endif
                       
                       write(message,'(A,A,A)') 'Initialise ',&
                                                  Trim(TargetVariablename),&
                                     ' using cubic spline interpolation'
                        CALL INFO(SolverName,Trim(message),Level=5)
                        call csa_interpolate_points(nin, pin, nout, pout, nppc ,csak)
                  CASE ('l')
                      write(message,'(A,A,A)') 'Initialise ',&
                                                  Trim(TargetVariablename),&
                                    ' using Linear interpolation'
                       CALL INFO(SolverName,Trim(message),Level=5)
                      call lpi_interpolate_points(nin, pin, nout, pout)
                  CASE DEFAULT
                    SELECT CASE (method(2:2))
                        CASE ('s')
                         nn_rule=1
                         write(message,'(A,A,A)') 'Initialise ',&
                                                  Trim(TargetVariablename),&
                 ' using Natural Neighbours Non-Sibsonian interpolation'
                        CASE DEFAULT
                         write(message,'(A,A,A)') 'Initialise ',&
                                                  Trim(TargetVariablename),&
                       ' using Natural Neighbours Sibson interpolation'
                     END SELECT
                   CALL INFO(SolverName,Trim(message), Level=5)
                   call nnpi_interpolate_points(nin,pin,w,nout,pout)
           END SELECT
            
             CALL INFO(SolverName,'-----Interpolation Done---', Level=5)
            !update variable value
            nNaN=0
            If (CheckNaN) Then
              Do i=1,nout
               z=pout(i)%z
                  If (isnan(z)) Then
                     nNaN=nNaN+1
                     If (ReplaceNaN) then
                         z=NaNVal
                     Else
                         x=pout(i)%x
                         y=pout(i)%y
                         kmin=1
                         MinD=sqrt((x-pin(1)%x)*(x-pin(1)%x)+ &
                                (y-pin(1)%y)*(y-pin(1)%y))
                         Do k=2,nin
                            D=sqrt((x-pin(k)%x)*(x-pin(k)%x)+ &
                               (y-pin(k)%y)*(y-pin(k)%y))
                          If (D.LT.MinD) then
                           kmin=k
                           MinD=D
                          End if
                         End Do
                         z=pin(kmin)%z
                      End IF
                   End IF
                 pout(i)%z=z
               !Values(Perm(i))=z
              End do
            End IF
            If (nNaN.GT.0) then
                If (ReplaceNaN) then
                   write(message,'(I0,A,A,e14.7)') nNaN, &
                                 ' values where NaN and have', &
                                 ' been replaced by ',NaNVal
                Else
                   write(message,'(I0,A,A)') nNaN, &
                                 ' values where NaN and have', &
                                 ' been replaced by nearest data value'
                EndIf
                CALL INFO(SolverName,Trim(message), Level=5)
             End If
             
             DO t=1,Model%Mesh%NumberOfBulkElements
                Element => Model%Mesh%Elements(t)
                DO i = 1, Element % TYPE % NumberOfNodes
                   nodeind = Element % NodeIndexes(i)
                   IF (Var % TYPE == Variable_on_nodes_on_elements) THEN
                      Varind=Element % DGIndexes(i)
                   ELSE
                      Varind=Element % NodeIndexes(i)
                   ENDIF
                   k=Perm(Varind)
                   IF (k.GT.0) THEN
                     z = pout(nodeind)%z
                     IF (HaveMin) z=Max(MinMaxVals(1),z)
                     IF (HaveMax) z=Min(MinMaxVals(2),z)
                     Values(Perm(Varind))= z 
                   END IF
                END DO
             END DO
            deallocate(pin)
        END DO

        deallocate(pout)

       CALL INFO(SolverName,'-----ALL DONE----------',Level=5)

       CONTAINS

       ! Find the min and max indexes of values within bBox
       SUBROUTINE MinMaxIndex(x,n,minx,maxx,MinIndex,MaxIndex)
       IMPLICIT NONE
       REAL(KIND=dp),INTENT(IN) :: x(:),minx,maxx
       INTEGER,INTENT(IN) :: n
       INTEGER,INTENT(OUT) :: MinIndex,MaxIndex

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
             
       END SUBROUTINE Scattered2DDataInterpolator
