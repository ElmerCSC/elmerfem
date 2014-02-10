program foobar

    use ElmerParam

    implicit none
    double precision :: xr(20)
    double precision :: x, y
    integer :: i, xi(3)
    character(512) :: line

    do i = 1, size(xr)/2
        xr(i) = i
        xr(i+size(xr)/2) = i/10.0d0
    end do

    y = elmer_param(xr)
    if (y /= 5.0d0) stop 1

    open(unit=10, file="foo")
    do i = 1, 7
        read (10, *) x, y
        if (abs(x - i/10.0d0) > 1.0e-6 .or. y /= i) stop 2
    end do
    close(10)

    open(unit=10, file="save.dat", status='old', action="read")
    ! I don't understand why it is necessary to read to a string first, and then
    ! read from the string (I get 'Fortran runtime error: End of file' if I try
    ! to read xi, xr and y from "save.dat" directly.  Compiler bug?  (gfortran
    ! 4.2.0 20060729)).
    read (10, '(A)') line
    read (line, *) xi, xr, y
    close(10)

    if (any(xi /= (/ 1, 2, 3 /))) stop 3
    if (maxval(xr(1:10) - (/ 1.0d0, 2.0d0, 3.0d0, 4.5d0, 3.1d0, 7.9d0, 7.0d0, &
                             8.0d0,  9.0d0, 10.0d0 /)) > 1.0e-4) stop 4
    if (y /= 5.0d0) stop 5

end program foobar

