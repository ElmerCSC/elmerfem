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
!  Interpolate scattered 2D data readed from an ASCII file (x y Value)
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

        CHARACTER(LEN=MAX_NAME_LEN) :: TargetVariableName,VariableName,DataF
        CHARACTER(LEN=MAX_NAME_LEN) :: Name,FName,WName,MName,CSVName,tmpName
        CHARACTER(LEN=MAX_NAME_LEN),parameter :: &
                         SolverName='Scattered2DDataInterpolator'
        CHARACTER(LEN=MAX_NAME_LEN) :: Xdim,dimName,FillName
        CHARACTER(2) :: method

        LOGICAL :: GotVar,GotTVar,Found,UnFoundFatal=.TRUE.
        LOGICAL :: CheckBBox,CheckNaN,ReplaceNaN,NETCDFFormat
        LOGICAL :: GoodVal,HaveFillv
        LOGICAL,dimension(:,:), allocatable :: mask
        LOGICAL,dimension(:), allocatable :: BBox
        LOGICAL :: Debug=.False.

        ! Variables to pass to the nn C library
        type(POINT),dimension(:),allocatable :: pout,pin,ptmp
        INTEGER(C_INT) :: nout,nin,nppc,csak
        REAL(C_DOUBLE) :: W

       CALL INFO(Trim(SolverName), &
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
              CALL INFO(Trim(SolverName),'Data File is in netcdf format', Level=5)
            Else
               CALL INFO(Trim(SolverName),'Data File is in ascii', Level=5)
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
            
           IF (.NOT.NETCDFFormat) Then
               open(unit = io, file = TRIM(DataF), status = 'old',iostat = ok)

               if(ok /= 0) then
                  write(message,'(A,A)') 'Unable to open file ',TRIM(DataF)
                  CALL Fatal(Trim(SolverName),Trim(message))
               end if
            
               nin=0
               !count the line number in the file
               do while(ok == 0) 
                 read(io,*,iostat = ok)
                 if (ok == 0) nin = nin + 1
               end do 

               IF (nin.eq.0) then
                   write(message,'(A,A)') 'No Data found in',TRIM(DataF)
                   CALL Fatal(Trim(SolverName),Trim(message))
               ENDIF
                                        
          
               allocate(pin(nin))

               ! comes back to beginning of file
               rewind(unit=io,iostat=ok)

               ! read datas
               do i = 1, nin
                   read(io,*,iostat = ok) pin(i)%x,pin(i)%y,pin(i)%z
               end do
               close(io)
           ELSE
               NetCDFstatus = NF90_OPEN(trim(DataF),NF90_NOWRITE,ncid)
               IF ( NetCDFstatus /= NF90_NOERR ) THEN
                   CALL Fatal(Trim(SolverName), &
                       'Unable to open NETCDF File')
               END IF
               WRITE (dimName,'(A,I0,A)') 'Variable ',&
                         NoVar,' x-dim Name'
               Xdim=ListGetString( Params, TRIM(dimName), Found )
               if (.NOT.Found) Xdim='x'
               NetCDFstatus = nf90_inq_dimid(ncid, trim(Xdim) , varid)
               IF ( NetCDFstatus /= NF90_NOERR ) THEN
                    CALL Fatal(Trim(SolverName), &
                        'Unable to  get netcdf x-dim Id')
               ENDIF
               NetCDFstatus = nf90_inquire_dimension(ncid,varid,len=nx)
               IF ( NetCDFstatus /= NF90_NOERR ) THEN
                    CALL Fatal(Trim(SolverName), &
                        'Unable to  get netcdf nx')
               ENDIF
               WRITE (dimName,'(A,I0,A)') 'Variable ',&
                         NoVar,' y-dim Name'
               Xdim=ListGetString( Params, TRIM(dimName), Found )
               if (.NOT.Found) Xdim='y'
               NetCDFstatus = nf90_inq_dimid(ncid, trim(Xdim) , varid)
               IF ( NetCDFstatus /= NF90_NOERR ) THEN
                    CALL Fatal(Trim(SolverName), &
                        'Unable to  get netcdf y-dim Id')
               ENDIF
               NetCDFstatus = nf90_inquire_dimension(ncid,varid,len=ny)
               IF ( NetCDFstatus /= NF90_NOERR ) THEN
                    CALL Fatal(Trim(SolverName), &
                        'Unable to  get netcdf ny')
               ENDIF
               !! allocate good size
               allocate(xx(nx),yy(ny),DEM(nx,ny),mask(nx,ny))
               !! Get the variable
               NetCDFstatus = nf90_inq_varid(ncid,TRIM(VariableName),varid)
               IF ( NetCDFstatus /= NF90_NOERR ) THEN
                   CALL Fatal(Trim(SolverName), &
                        'Unable to get netcdf variable id')
               ENDIF
               NetCDFstatus = nf90_get_var(ncid, varid,DEM)
               IF ( NetCDFstatus /= NF90_NOERR ) THEN
                   CALL Fatal(Trim(SolverName), &
                        'Unable to get netcdf variable')
               ENDIF
               HaveFillV=.True.
               NetCDFstatus = nf90_get_att(ncid, varid,"_FillValue",fillv)
               IF ( NetCDFstatus /= NF90_NOERR ) THEN
                   WRITE (FillName,'(A,I0,A)') 'Variable ',&
                           NoVar,' Fill Value'
                   fillv=ListGetConstReal(Params, TRIM(FillName) , Found)
                   if (.NOT.Found) HaveFillV=.False.
               ENDIF

               !! Get X variable
               WRITE (dimName,'(A,I0,A)') 'Variable ',&
                       NoVar,'x-Var Name'
               Xdim=ListGetString( Params, TRIM(dimName), Found )
               if (.NOT.Found) Xdim='x'
               NetCDFstatus = nf90_inq_varid(ncid,trim(Xdim),varid)
               IF ( NetCDFstatus /= NF90_NOERR ) THEN
                   CALL Fatal(Trim(SolverName), &
                        'Unable to get netcdf x-variable id')
               ENDIF
               NetCDFstatus = nf90_get_var(ncid, varid,xx)
               IF ( NetCDFstatus /= NF90_NOERR ) THEN
                   CALL Fatal(Trim(SolverName), &
                        'Unable to get netcdf x-variable ')
               ENDIF
               !! Get Y variable
               WRITE (dimName,'(A,I0,A)') 'Variable ',&
                       NoVar,'y-Var Name'
               Xdim=ListGetString( Params, TRIM(dimName), Found )
               if (.NOT.Found) Xdim='y'
               NetCDFstatus = nf90_inq_varid(ncid,trim(Xdim),varid)
               IF ( NetCDFstatus /= NF90_NOERR ) THEN
                   CALL Fatal(Trim(SolverName), &
                        'Unable to get netcdf y-variable id')
               ENDIF
               NetCDFstatus = nf90_get_var(ncid, varid,yy)
               IF ( NetCDFstatus /= NF90_NOERR ) THEN
                   CALL Fatal(Trim(SolverName), &
                        'Unable to get netcdf y-variable')
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
              CALL INFO(Trim(SolverName),Trim(message),Level=5)

               if (nin.eq.0) CALL Fatal(Trim(SolverName), &
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
                          CALL Fatal(Trim(SolverName),&
                               'get more Non-Nan values than expected')
                       ENDIF
                       pin(compt)%x = xx(i)
                       pin(compt)%y = yy(j)
                       pin(compt)%z = DEM(i,j)
                       IF (DEBUG) write(10,*) pin(compt)%x,pin(compt)%y,pin(compt)%z
                     endif         
                  End do
               End do
               IF (DEBUG) close(10) !tmp
               if (compt.ne.nin) CALL Fatal(Trim(SolverName),&
                       'sorry I didn t found the good number of values')
               deallocate(xx,yy,DEM,mask)
           ENDIF

           IF (CheckBBox) THEN
            allocate(Bbox(nin),ptmp(nin))
            ptmp=pin
            deallocate(pin)
            Bbox(:)=((ptmp(:)%x.ge.xmin).and.(ptmp(:)%x.le.xmax)& 
                       .and.(ptmp(:)%y.ge.ymin).and.(ptmp(:)%y.le.ymax))
            nval=COUNT(Bbox)
            if (nval.eq.0) CALL Fatal(Trim(SolverName), &
                                  'no data within bounding box??')
            write(message,'(A,I0,A,I0)') 'I Found ',nval,&
                     ' data points within bounding box over ',nin
              CALL INFO(Trim(SolverName),Trim(message),Level=5)
            allocate(pin(nval))
            compt=0
            Do i=1,nin
               if ((ptmp(i)%x.ge.xmin).and.(ptmp(i)%x.le.xmax)&
                   .and.(ptmp(i)%y.ge.ymin).and.(ptmp(i)%y.le.ymax)) then
                  compt=compt+1
                  pin(compt)=ptmp(i)
                endif
            End do
            nin=nval
            if (compt.ne.nin) CALL Fatal(Trim(SolverName),&
                 'BBox: sorry I didn t found the good number of values')
            deallocate(Bbox,ptmp)
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
                        CALL INFO(Trim(SolverName),Trim(message),Level=5)
                        call csa_interpolate_points(nin, pin, nout, pout, nppc ,csak)
                  CASE ('l')
                      write(message,'(A,A,A)') 'Initialise ',&
                                                  Trim(TargetVariablename),&
                                    ' using Linear interpolation'
                       CALL INFO(Trim(SolverName),Trim(message),Level=5)
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
                   CALL INFO(Trim(SolverName),Trim(message), Level=5)
                   call nnpi_interpolate_points(nin,pin,w,nout,pout)
           END SELECT
            
             CALL INFO(Trim(SolverName),'-----Interpolation Done---', Level=5)
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
                CALL INFO(Trim(SolverName),Trim(message), Level=5)
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
                   Values(Perm(Varind))=pout(nodeind)%z
                END DO
             END DO
            deallocate(pin)
        END DO

        deallocate(pout)

       CALL INFO(Trim(SolverName), &
           '-----ALL DONE----------',Level=5)

        END SUBROUTINE Scattered2DDataInterpolator
