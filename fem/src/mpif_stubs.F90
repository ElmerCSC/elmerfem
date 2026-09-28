MODULE mpi_stub
  IMPLICIT NONE
  LOGICAL :: mpi_stub_init = .FALSE.
END MODULE mpi_stub


SUBROUTINE mpi_init(ierr)
  USE mpi_stub
  IMPLICIT NONE
  INTEGER :: ierr
  ierr = 0
  mpi_stub_init = .TRUE.
END SUBROUTINE mpi_init

SUBROUTINE mpi_initialized(init, ierr)
  IMPLICIT NONE
  LOGICAL :: init
  INTEGER :: ierr
  ierr = 0
  init = .FALSE.
END SUBROUTINE mpi_initialized

SUBROUTINE mpi_init_thread(req, prov, ierr)
  IMPLICIT NONE
  INTEGER :: req, prov, ierr
  prov = req
  ierr = 0
END SUBROUTINE mpi_init_thread

SUBROUTINE mpi_finalize(ierr)
  USE mpi_stub
  IMPLICIT NONE
  INTEGER :: ierr
  ierr = 0
  mpi_stub_init = .FALSE.
END SUBROUTINE mpi_finalize

SUBROUTINE mpi_comm_size(comm, csize, ierr)
  IMPLICIT NONE
  INTEGER :: comm, csize, ierr
  ierr = 0
  csize = 1
END SUBROUTINE mpi_comm_size

SUBROUTINE mpi_comm_rank(comm, rank, ierr)
  IMPLICIT NONE
  INTEGER :: comm, rank, ierr
  rank = 0
  ierr = 0
END SUBROUTINE mpi_comm_rank

SUBROUTINE mpi_comm_split(comm, color, key, newcomm, ierr)
  IMPLICIT NONE
  INTEGER :: comm, color, key, newcomm, ierr
  newcomm = comm
  ierr = 0
END SUBROUTINE mpi_comm_split

SUBROUTINE mpi_scan
  IMPLICIT NONE
END SUBROUTINE mpi_scan

SUBROUTINE mpi_allreduce
  IMPLICIT NONE
  RETURN
END SUBROUTINE mpi_allreduce

SUBROUTINE mpi_buffer_detach
  IMPLICIT NONE
  RETURN
END SUBROUTINE mpi_buffer_detach

SUBROUTINE mpi_recv
  IMPLICIT NONE
  RETURN
END SUBROUTINE mpi_recv

SUBROUTINE mpi_buffer_attach
  IMPLICIT NONE
  RETURN
END SUBROUTINE mpi_buffer_attach

SUBROUTINE mpi_send
  IMPLICIT NONE
  RETURN
END SUBROUTINE mpi_send

SUBROUTINE mpi_barrier
  IMPLICIT NONE
  RETURN
END SUBROUTINE mpi_barrier

SUBROUTINE mpi_wait
  IMPLICIT NONE
  RETURN
END SUBROUTINE mpi_wait

SUBROUTINE mpi_waitany
  IMPLICIT NONE
  RETURN
END SUBROUTINE mpi_waitany

SUBROUTINE mpi_bsend
  IMPLICIT NONE
  RETURN
END SUBROUTINE mpi_bsend

SUBROUTINE mpi_comm_free
  IMPLICIT NONE
  RETURN
END SUBROUTINE mpi_comm_free

SUBROUTINE mpi_waitall
  IMPLICIT NONE
  RETURN
END SUBROUTINE mpi_waitall

SUBROUTINE mpi_comm_group
  IMPLICIT NONE
  RETURN
END SUBROUTINE mpi_comm_group

SUBROUTINE mpi_group_incl
  IMPLICIT NONE
  RETURN
END SUBROUTINE mpi_group_incl

SUBROUTINE mpi_comm_create
  IMPLICIT NONE
  RETURN
END SUBROUTINE mpi_comm_create

SUBROUTINE mpi_irecv
  IMPLICIT NONE
  RETURN
END SUBROUTINE mpi_irecv

SUBROUTINE mpi_isend
  IMPLICIT NONE
  RETURN
END SUBROUTINE mpi_isend

SUBROUTINE mpi_bcast
  IMPLICIT NONE
  RETURN
END SUBROUTINE mpi_bcast

SUBROUTINE mpi_allgather
  IMPLICIT NONE
  RETURN
END SUBROUTINE mpi_allgather

SUBROUTINE mpi_allgatherv
  IMPLICIT NONE
  RETURN
END SUBROUTINE mpi_allgatherv

SUBROUTINE mpi_alltoallv
  IMPLICIT NONE
  RETURN
END SUBROUTINE mpi_alltoallv

SUBROUTINE mpi_alltoall
  IMPLICIT NONE
  RETURN
END SUBROUTINE mpi_alltoall

SUBROUTINE mpi_gatherv
  IMPLICIT NONE
  RETURN
END SUBROUTINE mpi_gatherv

SUBROUTINE mpi_gather
  IMPLICIT NONE
  RETURN
END SUBROUTINE mpi_gather

SUBROUTINE mpi_reduce
  IMPLICIT NONE
  RETURN
END SUBROUTINE mpi_reduce

SUBROUTINE mpi_get_count
  IMPLICIT NONE
  RETURN
END SUBROUTINE mpi_get_count

! Parpack 
SUBROUTINE pdseupd
  IMPLICIT NONE
  RETURN
END SUBROUTINE pdseupd

SUBROUTINE pdneupd
  IMPLICIT NONE
  RETURN
END SUBROUTINE pdneupd

SUBROUTINE pzneupd
  IMPLICIT NONE
  RETURN
END SUBROUTINE pzneupd

SUBROUTINE pdsaupd
  IMPLICIT NONE
  RETURN
END SUBROUTINE pdsaupd

SUBROUTINE pdnaupd
  IMPLICIT NONE
  RETURN
END SUBROUTINE pdnaupd

SUBROUTINE pznaupd
  IMPLICIT NONE
  RETURN
END SUBROUTINE pznaupd


