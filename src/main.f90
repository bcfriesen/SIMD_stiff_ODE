program main

    use rkf45_mod
    implicit none

    integer, parameter :: dp = kind(1.0d0)

    integer, parameter :: simd_width = 4 ! For Edison
    real(kind=dp) :: y0(simd_width), t0, tf, dt0, eps
    real(kind=dp) :: t, dt, y(simd_width)
    integer :: flag

    y0 = 0.0d0
    ! Every ODE must use the same limits of integration.
    t0 = 0.0d0
    tf = 5.0d0
    ! Use the same initial guess for the time step.
    dt0 = 25.0d0
    ! Local truncation error tolerance.
    eps = 1.0d-3

    write (*, *), 'Integrating the equation y''(t) = 2*t from:'
    write (*, '(a3, es12.3e2, a3, es12.3e2, a3)') '[', t0, ', ', tf, ']'
    write (*, '(a15, es12.3e2)') 'LTE: ', eps

    write (*, *) 'Initial values:'
    write (*, '(a8, 4es12.3e2)') 'y0 = ', y0
    write (*, '(a8, es12.3e2)') 't0 = ', t0
    write (*, '(a8, es12.3e2)') 'dt0 = ', dt0

    call rkf45(y0(:), t0, tf, dt0, eps, t, dt, y(:), flag)

    if (flag == 0) then
      write (*, *) 'Success!'
      write (*, *) 'Final values:'
      write (*, '(a8, 4es12.3e2)') 'y = ', y
      write (*, '(a8, es12.3e2)')  't = ', t
      write (*, '(a8, es12.3e2)')  'dt = ', dt
    else
      print *, 'ERROR: Integration failed!'
    end if

end program main
