        FUNCTION UIni( Model, nodenumber, dumy) RESULT(U)
        USE types

        implicit none
        TYPE(Model_t) :: Model
        Real(kind=dp) :: dumy,U
        INTEGER :: nodenumber

        Real(kind=dp),allocatable :: dem(:,:),xx(:),yy(:)
        Real(kind=dp) :: x,y
        Real(kind=dp) :: LinearInterp

        integer :: nx,ny
        integer :: i,j

        character(len=MAX_NAME_LEN) :: filin='PROG/UDEM.xy'

        logical :: Firsttime=.true.

        SAVE dem,xx,yy,nx,ny
        SAVE Firsttime

        if (Firsttime) then
                Firsttime=.False.

        ! open file
                open(10,file=trim(filin))
                Read(10,*) nx
                Read(10,*) ny
                allocate(xx(nx),yy(ny))
                Allocate(dem(nx,ny))
                Do i=1,nx
                   Do j=1,ny
                      read(10,*) xx(i),yy(j),dem(i,j)
                   End Do
                End do
                close(10)
        End if

        ! position current point
        x=Model % Nodes % x (nodenumber)
        y=Model % Nodes % y (nodenumber)

        U=LinearInterp(dem,xx,yy,nx,ny,x,y)

        Return 
        End

        FUNCTION VIni( Model, nodenumber, dumy) RESULT(U)
        USE types

        implicit none
        TYPE(Model_t) :: Model
        Real(kind=dp) :: dumy,U
        INTEGER :: nodenumber

        Real(kind=dp),allocatable :: dem(:,:),xx(:),yy(:)
        Real(kind=dp) :: x,y
        Real(kind=dp) :: LinearInterp

        integer :: nx,ny
        integer :: i,j

        character(len=MAX_NAME_LEN) :: filin='PROG/VDEM.xy'

        logical :: Firsttime=.true.

        SAVE dem,xx,yy,nx,ny
        SAVE Firsttime

        if (Firsttime) then
                Firsttime=.False.

        ! open file
                open(10,file=trim(filin))
                Read(10,*) nx
                Read(10,*) ny
                allocate(xx(nx),yy(ny))
                Allocate(dem(nx,ny))
                Do i=1,nx
                   Do j=1,ny
                      read(10,*) xx(i),yy(j),dem(i,j)
                   End Do
                End do
                close(10)
        End if

        ! position current point
        x=Model % Nodes % x (nodenumber)
        y=Model % Nodes % y (nodenumber)

        U=LinearInterp(dem,xx,yy,nx,ny,x,y)

        Return 
        End

        FUNCTION zsIni( Model, nodenumber, dumy) RESULT(U)
        USE types

        implicit none
        TYPE(Model_t) :: Model
        Real(kind=dp) :: dumy,U
        INTEGER :: nodenumber

        Real(kind=dp),allocatable :: dem(:,:),xx(:),yy(:)
        Real(kind=dp) :: x,y
        Real(kind=dp) :: LinearInterp

        integer :: nx,ny
        integer :: i,j

        character(len=MAX_NAME_LEN) :: filin='PROG/zsDEM.xy'

        logical :: Firsttime=.true.

        SAVE dem,xx,yy,nx,ny
        SAVE Firsttime

        if (Firsttime) then
                Firsttime=.False.

        ! open file
                open(10,file=trim(filin))
                Read(10,*) nx
                Read(10,*) ny
                allocate(xx(nx),yy(ny))
                Allocate(dem(nx,ny))
                Do i=1,nx
                   Do j=1,ny
                      read(10,*) xx(i),yy(j),dem(i,j)
                   End Do
                End do
                close(10)
        End if

        ! position current point
        x=Model % Nodes % x (nodenumber)
        y=Model % Nodes % y (nodenumber)

        U=LinearInterp(dem,xx,yy,nx,ny,x,y)

        Return 
        End

        FUNCTION zbIni( Model, nodenumber, dumy) RESULT(U)
        USE types

        implicit none
        TYPE(Model_t) :: Model
        Real(kind=dp) :: dumy,U
        INTEGER :: nodenumber

        Real(kind=dp),allocatable :: dem(:,:),xx(:),yy(:)
        Real(kind=dp) :: x,y
        Real(kind=dp) :: LinearInterp

        integer :: nx,ny
        integer :: i,j

        character(len=MAX_NAME_LEN) :: filin='PROG/zbDEM.xy'

        logical :: Firsttime=.true.

        SAVE dem,xx,yy,nx,ny
        SAVE Firsttime

        if (Firsttime) then
                Firsttime=.False.

        ! open file
                open(10,file=trim(filin))
                Read(10,*) nx
                Read(10,*) ny
                allocate(xx(nx),yy(ny))
                Allocate(dem(nx,ny))
                Do i=1,nx
                   Do j=1,ny
                      read(10,*) xx(i),yy(j),dem(i,j)
                   End Do
                End do
                close(10)
        End if

        ! position current point
        x=Model % Nodes % x (nodenumber)
        y=Model % Nodes % y (nodenumber)

        U=LinearInterp(dem,xx,yy,nx,ny,x,y)

        Return 
        End


        Function LinearInterp(dem,xx,yy,nx,ny,x,y) Result(InterP1)

        USE TYPES
        implicit none
        
        REAL(KIND=dp) :: dem(nx,ny),xx(nx),yy(ny)
        REAL(KIND=dp) :: Dx,Dy,DxDy
        Real(kind=dp) :: x,y,x_1,y_1,B(4)
        Real(kind=dp) :: InterP1
        integer :: nx,ny
        integer :: nx_1,ny_1

        Dx=(xx(nx)-xx(1))/(nx-1)
        Dy=(yy(ny)-yy(1))/(ny-1)
        DxDy=Dx*Dy

        ! lower left point in DEM
        nx_1=floor((x-xx(1))/Dx) + 1
        ny_1=floor((y-yy(1))/Dy) + 1
        nx_1=min(nx_1,nx-1)
        ny_1=min(ny_1,ny-1)

        x_1=xx(nx_1)
        y_1=yy(ny_1)


        ! DEM Value in surroundings points
        !       4 ----- 3
        !       |       |
        !       1 ----- 2
        B(1)=dem(nx_1,ny_1)
        B(2)=dem(nx_1+1,ny_1)
        B(3)=dem(nx_1+1,ny_1+1)
        B(4)=dem(nx_1,ny_1+1)

        
        ! Linear Interpolation at Point x,y
        InterP1=(x-x_1)*(y-y_1)*(B(3)+B(1)-B(2)-B(4))/DxDy
        InterP1=InterP1+(x-x_1)*(B(2)-B(1))/Dx+(y-y_1)*(B(4)-B(1))/Dy+B(1)

        Return
        End
