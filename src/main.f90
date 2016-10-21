program main

    use rkf45_mod
    use rkf45_scalar_mod
    implicit none

    integer, parameter :: dp = kind(1.0d0)

    integer, parameter :: dp_simd_width = 4 ! AVX2 (Ivy Bridge) has 256 bit-wide SIMD width = 4 doubles
    integer, parameter :: array_length = 131072
    real(kind=dp) :: y0(array_length), t0, tf, dt0, eps
    real(kind=dp) :: t, dt, y(array_length)
    real(kind=dp) :: t1, t2
    integer :: num_steps
    integer :: flag
    integer :: i

    ! Set different initial conditions.
    do i = 1, array_length
      y0(i) = real(i, kind=dp)
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

    ! First do the integrations in scalar.
    call cpu_time(t1)
    do i = 1, array_length
      call rkf45_scalar(y0(i), t0, tf, dt0, eps, t, dt, y(i), flag, num_steps)
    end do
    call cpu_time(t2)
    write (*, *) 'time in scalar RKF45 (sec): ', t2-t1
    open(unit=11, name='scalar_results.dat')
    do i = 1, array_length
      write (11, '(i8, 3es18.6e2)') i, y(i), t, dt
    end do
    close(11)

    ! Now do the same integrations in SIMD.
    call cpu_time(t1)
    do i = 1, array_length, dp_simd_width
      call rkf45(y0(i:i+dp_simd_width-1), t0, tf, dt0, eps, t, dt, y(i:i+dp_simd_width), flag, num_steps)
    end do
    call cpu_time(t2)
    write (*, *) 'time in SIMD RKF45 (sec): ', t2-t1

    if (flag == 0) then
      write (*, *) 'Success!'
      write (*, '(a25, i8)') '# of time steps taken: ', num_steps
      open(unit=11, name='simd_results.dat')
      do i = 1, array_length
        write (11, '(i8, 3es18.6e2)') i, y(i), t, dt
      end do
      close(11)
    else
      write (*, *) 'ERROR: Integration failed!'
    end if

end program main
