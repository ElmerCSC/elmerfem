program Compare

  !***************************************************
  !                    Comparison outputs
  ! Format outputs En.m
  ! n = #caracteres
  ! m = nber digits
  !
  ! parm1, parm2, parm3, pamr4 -> shell - inputs 
  ! 
  !
  ! Results in file difference.txt for each test-case:
  ! if x<1E-08, Absolute deviation
  !   if difference>1E-8 print "ERROR abs" end
  ! else relative gap
  !   if difference>1E-8 print "ERROR " end
  ! end
  !
  !***************************************************
  
  implicit none
  
  !Variables
  CHARACTER(LEN=150) :: parm1,parm2,parm3,parm4,file1,file2
  CHARACTER(LEN=2) :: str4
  CHARACTER(LEN=150) :: filename1, filename2, my_fmt
  INTEGER :: ios, ios2, ierr, ierr2
  INTEGER :: n_rows, n_rows2
  INTEGER :: i, j, n
  INTEGER :: nb_arg, n_col_cpu
  REAL :: read_value, read_value2, x_L_null
  REAL, DIMENSION(:,:), ALLOCATABLE :: X_L, X_m
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ecart

  call getarg(1,parm1) !name file computer
  read (parm1,*) file1
  call getarg(2,parm2) !name file
  read (parm2,*) file2
  call getarg(3,parm3) !nber of inputs to compare
  read (parm3,*) nb_arg
  call getarg(4,parm4) !column nber of time_cpu
  read (parm4,*) n_col_cpu

  !Concatenation file.txt
  filename1=trim(file1)//'.txt'
  filename2=trim(file2)//'.txt'
  
  !Format of digits
  write(str4,'(I2)') nb_arg
  my_fmt='(' // trim(str4) // 'E20.12)'

  !Number of lines in each file
  OPEN(unit=20,file=trim(filename1),status='old',action='read')
  ios = 0
  n_rows = 0
  do while (.true.)
     read(unit=20,fmt=*,iostat=ios) read_value
     if (ios .ne. 0) exit
     n_rows = n_rows + 1
  enddo
  CLOSE(20)

  OPEN(unit=30,file='./DATA/'//trim(filename2),status='old')
  ios2 = 0
  n_rows2 = 0
  do while (.true.)
     read(unit=30,fmt=*,iostat=ios2) read_value2
     if (ios2 .ne. 0) exit
     n_rows2 = n_rows2 + 1
  enddo
  CLOSE(30)
  if (n_rows /= n_rows2) then
     print*, 'ERROR: nb_line is not the same'
  end if

  !Allocation
  Allocate(x_L(1:n_rows,1:nb_arg))
  Allocate(x_m(1:n_rows,1:nb_arg))
  Allocate(ecart(1:n_rows,1:nb_arg))

  !Reading File - 
  OPEN(unit=40,file='./DATA/'//trim(filename2),iostat=ierr)
  if (ierr/=0) then
     print*, 'ERROR: can not open file '
     stop
  endif
  rewind(unit=40,iostat=ierr)
  Do i=1,n_rows
     read(40,*) (x_L(i,j),j=1,nb_arg)
  End Do
  CLOSE(40)
  
  !Reading File - your computer
  OPEN(unit=50,file=trim(filename1),iostat=ierr2)
  if (ierr2/=0) then
     print*, 'ERROR: le fichier '
     stop
  endif
  rewind(unit=50,iostat=ierr2)
  Do i=1,n_rows    
     read(50,*) (x_m(i,n), n=1,nb_arg)
  End do
  CLOSE(50)
  
  !Comparisons of outputs
  OPEN(60, file='difference.txt', action='write')
  Do i = 1,n_rows
     Do j=1,nb_arg
        if (j == n_col_cpu) then
           ecart(i,j)= abs((x_L(i,j)-x_m(i,j))/x_L(i,j))
           if (ecart(i,j)>0.9) then
              print*, 'DIFF-TIME'
           end if
        else
           if (x_L(i,j)<1e-8) then
              ecart(i,j)= abs(x_L(i,j)-x_m(i,j))
              if (ecart(i,j) > 1e-8) then
                 print*, 'ERROR abs ',ecart(i,j),'ligne',i,'colonne ',j
              end if
           else
              ecart(i,j)= abs((x_L(i,j)-x_m(i,j))/x_L(i,j))
           end if
           if (ecart(i,j) > 1e-8) then
              print*, 'ERROR ',ecart(i,j),'ligne',i,'colonne ',j
           end if
        end if
     End do
     WRITE(60,my_fmt)(ecart(i,j), j=1,nb_arg)
  End Do
  CLOSE(60)  
end program Compare
