!*****************************************************************************************
!>
!  Kinds for Minpack

    module minpack_kinds

    use iso_fortran_env, only: real64

    implicit none

    private

    ! Use the same double precision type as rest of Elmer
    integer,parameter,public :: wp = SELECTED_REAL_KIND(12)

!*****************************************************************************************
    end module minpack_kinds
!*****************************************************************************************
