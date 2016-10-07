program main

    use rkf45_mod
    implicit none

    integer, parameter :: dp = kind(1.0d0)

    integer, parameter :: simd_width = 4 ! For Edison
    real(kind=dp) :: y0(simd_width), t0, tf, dt0, eps
    real(kind=dp) :: t, dt, y(simd_width)
    integer :: flag

    y0 = 0.0d0
    t0 = 0.0d0
    tf = 5.0d0
    dt0 = 25.0d0
    eps = 1.0d-3

    print *, 'Integrating the equation y''(t) = 2*t from [', t0, ', ', tf, ']'

    print *, 'Initial values:'
    print *, 'y0 = ', y0
    print *, 't0 = ', t0
    print *, 'dt0 = ', dt0

    call rkf45(y0(:), t0, tf, dt0, eps, t, dt, y(:), flag)

    if (flag == 0) then
      print *, 'Success!'
      print *, 'Final values:'
      print *, 't = ', t
      print *, 'y = ', y
      print *, 'dt = ', dt
    else
      print *, 'ERROR: Integration failed!'
    end if

end program main
