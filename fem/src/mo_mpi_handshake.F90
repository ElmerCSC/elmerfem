MODULE mo_mpi_handshake


  PUBLIC MAX_GROUPNAME_LEN
  INTEGER, PARAMETER :: MAX_GROUPNAME_LEN = 64

CONTAINS

  SUBROUTINE mpi_handshake ( comm, group_names, group_comms )
    use, intrinsic :: iso_c_binding, only : c_ptr, c_char, c_null_char, c_loc

    implicit none

    interface
       SUBROUTINE mpi_handshake_c2f (n, group_names, group_comms, comm) &
          bind ( c, name='mpi_handshake_c2f')
         use, intrinsic :: iso_c_binding, only : c_ptr, c_int
         implicit none
         integer(c_int), intent(in), value :: n
         type (c_ptr) , intent(in) :: group_names(n)
         integer, intent(out) :: group_comms(n)
         integer, intent(in), value :: comm
       end subroutine mpi_handshake_c2f
    end interface

    integer, intent(in) :: comm
    character(len=MAX_GROUPNAME_LEN), intent(in) :: group_names(:)
    integer, intent(inout) :: group_comms(SIZE(group_names))

    CHARACTER (kind=c_char, len=MAX_GROUPNAME_LEN), TARGET :: group_names_cpy(SIZE(group_names))
    type( c_ptr ) :: group_name_ptr(SIZE(group_names))
    integer :: i
    DO i=1,SIZE(group_names)
       group_names_cpy(i) = TRIM(group_names(i)) // c_null_char
       group_name_ptr(i) = c_loc(group_names_cpy(i))
    END DO

    call mpi_handshake_c2f(SIZE(group_names), group_name_ptr, group_comms, comm)

  end subroutine mpi_handshake

END MODULE mo_mpi_handshake
