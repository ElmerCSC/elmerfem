! MinGW/MPICH2:
!  gfortran-4.4.exe -fsecond-underscore -I/c/MPICH2/include -L/c/MPICH2/lib mpitest.f90 -o mpitest -lfmpich2g
! mpiexec -localonly 2 mpitest

program main
  include 'mpif.h'
  integer rc
  integer num_procs
  integer id

  call MPI_Init(rc)
  call MPI_Comm_size(MPI_COMM_WORLD, num_procs, rc)
  call MPI_Comm_rank(MPI_COMM_WORLD, id, rc)

  if(id == 0) print *, 'Num procs: ', num_procs

  print *, 'Proc Id:   ', id
  
  call MPI_Finalize(rc)
  
end program main
