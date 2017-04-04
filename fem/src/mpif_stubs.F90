SUBROUTINE mpi_init(ierr)
  INTEGER :: ierr
  ierr = 0
END SUBROUTINE mpi_init

SUBROUTINE mpi_init_thread(req, prov, ierr)
  INTEGER :: req, prov, ierr
  prov = req
  ierr = 0
END SUBROUTINE mpi_init_thread

SUBROUTINE mpi_finalize(ierr)
  INTEGER :: ierr
  ierr = 0
END SUBROUTINE mpi_finalize

SUBROUTINE mpi_comm_size(comm, csize, ierr)
  INTEGER :: comm, csize, ierr
  ierr = 0
  csize = 1
END SUBROUTINE mpi_comm_size

SUBROUTINE mpi_comm_rank(comm, rank, ierr)
  INTEGER :: comm, rank, ierr
  rank = 0
  ierr = 0
END SUBROUTINE mpi_comm_rank

SUBROUTINE mpi_comm_split(comm, color, key, newcomm, ierr)
  INTEGER :: comm, color, key, newcomm, ierr
  newcomm = comm
  ierr = 0
END SUBROUTINE mpi_comm_split

SUBROUTINE mpi_allreduce
  RETURN
END SUBROUTINE mpi_allreduce

SUBROUTINE mpi_buffer_detach
  RETURN
END SUBROUTINE mpi_buffer_detach

SUBROUTINE mpi_recv
  RETURN
END SUBROUTINE mpi_recv

SUBROUTINE mpi_buffer_attach
  RETURN
END SUBROUTINE mpi_buffer_attach

SUBROUTINE mpi_send
  RETURN
END SUBROUTINE mpi_send

SUBROUTINE mpi_barrier
  RETURN
END SUBROUTINE mpi_barrier

SUBROUTINE mpi_waitany
  RETURN
END SUBROUTINE mpi_waitany

SUBROUTINE mpi_bsend
  RETURN
END SUBROUTINE mpi_bsend

SUBROUTINE mpi_comm_free
  RETURN
END SUBROUTINE mpi_comm_free

SUBROUTINE mpi_waitall
  RETURN
END SUBROUTINE mpi_waitall

SUBROUTINE mpi_comm_group
  RETURN
END SUBROUTINE mpi_comm_group

SUBROUTINE mpi_group_incl
  RETURN
END SUBROUTINE mpi_group_incl

SUBROUTINE mpi_comm_create
  RETURN
END SUBROUTINE mpi_comm_create

SUBROUTINE mpi_irecv
  RETURN
END SUBROUTINE mpi_irecv

! Parpack 
SUBROUTINE pdseupd
  RETURN
END SUBROUTINE pdseupd

SUBROUTINE pdneupd
  RETURN
END SUBROUTINE pdneupd

SUBROUTINE pzneupd
  RETURN
END SUBROUTINE pzneupd

SUBROUTINE pdsaupd
  RETURN
END SUBROUTINE pdsaupd

SUBROUTINE pdnaupd
  RETURN
END SUBROUTINE pdnaupd

SUBROUTINE pznaupd
  RETURN
END SUBROUTINE pznaupd
