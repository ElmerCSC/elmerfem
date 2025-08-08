MODULE mo_mpi_handshake_stub


  PUBLIC MAX_GROUPNAME_LEN
  INTEGER, PARAMETER :: MAX_GROUPNAME_LEN = 64

CONTAINS

  SUBROUTINE mpi_handshake ( comm, group_names, group_comms )

    implicit none

    integer, intent(in) :: comm
    character(len=MAX_GROUPNAME_LEN), intent(in) :: group_names(:)
    integer, intent(inout) :: group_comms(SIZE(group_names))

  end subroutine mpi_handshake

END MODULE mo_mpi_handshake_stub
