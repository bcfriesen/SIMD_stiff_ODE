program main

    use rkf45_mod
    implicit none

    integer, parameter :: dp = kind(1.0d0)

    integer, parameter :: dp_simd_width = 4 ! AVX2 (Ivy Bridge) has 256 bit-wide SIMD width = 4 doubles
    integer, parameter :: array_length = 8
    real(kind=dp) :: y0(array_length), t0, tf, dt0, eps
    real(kind=dp) :: t, dt, y(array_length)
    integer :: num_steps
    integer :: flag
    integer :: i

    ! Set different initial conditions.
    do i = 1, array_length
      y0(i) = real(i, dp)
    end do
    ! Every ODE must use the same limits of integration.
    t0 = 0.0d0
    tf = 5.0d0
    ! Use the same initial guess for the time step.
    dt0 = 25.0d0
    ! Local truncation error tolerance.
    eps = 1.0d-3

    write (*, *), 'Integrating the equation y''(t) = 2*t from:'
    write (*, '(a3, es12.3e2, a3, es12.3e2, a3)') '[', t0, ', ', tf, ']'
    write (*, '(a15, es18.6e2)') 'LTE: ', eps


    do i = 1, array_length, dp_simd_width
      call rkf45(y0(i:i+dp_simd_width-1), t0, tf, dt0, eps, t, dt, y(i:i+dp_simd_width), flag, num_steps)
    end do

    if (flag == 0) then
      write (*, *) 'Success!'
      write (*, '(a25, i8)') '# of time steps taken: ', num_steps
      write (*, *) 'Final values:'
      write (*, '(a8, 8es12.3e2)') 'y = ', y
      write (*, '(a8, es12.3e2)')  't = ', t
      write (*, '(a8, es12.3e2)')  'dt = ', dt
    else
      write (*, *) 'ERROR: Integration failed!'
    end if

end program main
