!
!  ElmerParam - A simple system for parametrized computing
! 
!  Copyright (C) 2006  CSC - IT Center for Science Ltd.
!
!  Authors: Erik Edelmann <Erik.Edelmann@csc.fi>
!           Peter Råback <Peter.Raback@csc.fi>
!  Address: CSC - IT Center for Science Ltd.
!           Keilaranta 14
!           02101 Espoo, Finland
!            
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU General Public License
!  as published by the Free Software Foundation; either version 2
!  of the License, or (at your option) any later version.
!  
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
! 
!  You should have received a copy of the GNU General Public License
!  along with this program (in file elmerparam/COPYING); if not, write to
!  the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
!  Boston, MA 02110-1301, USA.
!

!
! TODO: Use IOMSG instead of IOSTAT to report IO errors when more compilers
! support it (F2003 feature).
!
! TODO: Use get_command_argument() and command_argument_count() instead of
! getarg() and iargc when more compilers support them (F2003 features).
!

module parameters

    implicit none
    integer, parameter :: STDERR = 0
    integer, parameter :: STDIN = 5
    integer, parameter :: STDOUT = 6

    integer, parameter :: DP = selected_real_kind(15)

end module parameters


program elmerparam_cmdline

    use elmerparam
    use parameters

    implicit none
    integer :: nr, ni, nfun
    integer, allocatable :: xi(:)
    real(DP), allocatable :: xr(:), fun(:)
    character(256) :: tag, arg, infile, outfile
    integer :: inonopt
    integer :: iou, ios, i
    include "iargc_decl.inc"

    infile = "-"
    outfile = "-"
    tag = ""

    inonopt = 1  ! Index of non-option argument

    do i = 1, iargc()
        call getarg(i, arg)
        if (arg == "-h" .or. arg == "--help") then
            call write_help()
        else
            select case (inonopt)
            case(1)
                infile = arg
            case(2)
                outfile = arg
            case(3)
                tag = arg
            case default
                write (STDERR, *) 'Too many arguments'
                call write_help()
                STOP 1
            end select
            inonopt = inonopt + 1
        end if
    end do

    if (infile == "-") then
        iou = STDIN
    else
        iou = 10
        open (unit=iou, file=infile, action="read", status="old", iostat=ios)
        if (ios /= 0) then
            write (STDERR, *) "Can't open file ", trim(infile), &
                              "; error code ", ios
            STOP 1
        end if
    end if

    read (iou, *, iostat=ios) nr
    if (ios > 0) then
        write (STDERR, *) 'Error reading number of reals from file ', &
                          trim(infile), '; error code', ios
        STOP 1
    else if (ios < 0) then
        nr = 0
        ni = 0
        nfun = 1
        go to 100
    end if

    if (nr > 0) then
        allocate(xr(nr))
        read (iou, *, iostat=ios) xr
        if (ios > 0) then
            write (STDERR, *) 'Error reading array of reals from file ', &
                              trim(infile), '; error code', ios
            STOP 1
        else if (ios < 0) then
            write(STDERR,*)'Unexpected end of file while reading array of reals'
            STOP 1
        end if
    end if

    read (iou, *, iostat=ios) ni
    if (ios > 0) then
        write (STDERR, *) 'Error reading number of integers from file ', &
                          trim(infile), '; error code', ios
        STOP 1
    else if (ios < 0) then
        ni = 0
        nfun = 1
        go to 100
    end if

    if (ni > 0) then
        allocate(xi(ni))
        read (iou, *, iostat=ios) xi
        if (ios > 0) then
            write (STDERR, *) 'Error reading array of integers from file ', &
                              trim(infile), '; error code', ios
            STOP 1
        else if (ios < 0) then
            write (STDERR, *) 'Unexpected end of file while reading array of &
                              &integers'
            STOP 1
        end if
    end if

    read (iou, *, iostat=ios) nfun
    if (ios > 0) then
        write (STDERR, *) 'Error reading number of result values from file ', &
                          trim(infile), '; error code', ios
        STOP 1
    else if (ios < 0) then
        nfun = 1
        go to 100
    end if

100 close(iou)

    allocate(fun(nfun))
    if (.not.allocated(xr)) allocate(xr(0))
    if (.not.allocated(xi)) allocate(xi(0))

    fun = elmer_param(nfun, xr, xi, tag)

    if (outfile == "-") then
        iou = STDOUT
    else
        iou = 10
        open (unit=iou, file=outfile, action="write", status="unknown", &
              iostat=ios)
        if (ios /= 0) then
            write (STDERR, *) "Can't open file ", trim(outfile), &
                              " for writing; error code ", ios
            STOP 1
        end if
    end if

    write (iou, '(ES15.8)') fun

contains

    subroutine write_help()
        write (STDERR, '(A)') "Shell command interface to ElmerParam"
        write (STDERR, '(A)') "Usage:"
        write (STDERR, '(A)') "ElmerParam [-h|--help] [inputfile [outputfile &
                              &[tag]]]"
        write (STDERR, '(A)') ""
        write (STDERR, '(A)') "If inputfile is missing or -, stdin is used"
        write (STDERR, '(A)') "If outputfile is missing or -, stdout is used"
        write (STDERR, '(A)') "If tag is missing, '' is used"

    end subroutine write_help

end program elmerparam_cmdline
