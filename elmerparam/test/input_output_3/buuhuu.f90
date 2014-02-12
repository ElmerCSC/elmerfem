program buuhuu

    use elmerparam

    integer, parameter :: N = 10
    double precision :: x(N), y(N)
    integer :: i

    do i = 1, N
        x(i) = i
    end do

    y = elmer_param(N, x)

    if (any(y /= x)) stop 1

end program buuhuu
