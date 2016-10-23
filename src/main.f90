program main

    use rkf45_mod
    use rkf45_scalar_mod
    implicit none

    integer, parameter :: dp = kind(1.0d0)

    integer, parameter :: dp_simd_width = 4 ! AVX2 (Ivy Bridge) has 256 bit-wide SIMD width = 4 doubles
    integer, parameter :: array_length = 2**18
    real(kind=dp) :: y0(array_length), t0, tf, dt0, eps
    real(kind=dp) :: t, dt, y(array_length)
    real(kind=dp) :: t1, t2
    integer :: num_steps
    integer :: flag
    integer :: i

    ! Set initial conditions. Just use the value of the loop index as the
    ! initial condition.
    do i = 1, array_length
      y0(i) = real(i, kind=dp)
    end do

    ! Every ODE must use the same limits of integration. This is the case for the
    ! Nyx problem, where the Strang splitting time step is the same in every
    ! cell during each grid sweep.
    t0 = 0.0d0
    tf = 5.0d0

    ! Use the same initial guess for the time step for every ODE.
    dt0 = 25.0d0

    ! Local truncation error tolerance.
    eps = 1.0d-3

    write (*, *), 'Integrating the equation y''(t) = 2*t from:'
    write (*, '(a3, es12.3e2, a3, es12.3e2, a3)') '[', t0, ', ', tf, ']'
    write (*, '(a45, es18.6e2)') 'local truncation error tolerance: ', eps

    ! First do the integrations using the scalar version of RKF45.
    call cpu_time(t1)
    do i = 1, array_length
      call rkf45_scalar(y0(i), t0, tf, dt0, eps, t, dt, y(i), flag, num_steps)
    end do
    call cpu_time(t2)
    write (*, *) 'time in scalar RKF45 (sec): ', t2-t1

    ! Now do the same sweep of integrations using the SIMD RKF45.
    call cpu_time(t1)
    do i = 1, array_length, dp_simd_width
      call rkf45(y0(i:i+dp_simd_width-1), t0, tf, dt0, eps, t, dt, y(i:i+dp_simd_width-1), flag, num_steps)
    end do
    call cpu_time(t2)
    write (*, *) 'time in SIMD RKF45 (sec): ', t2-t1

    if (flag == 0) then
      write (*, *) 'Success!'
      write (*, '(a25, i8)') '# of time steps taken: ', num_steps
    else
      write (*, *) 'ERROR: Integration failed!'
    end if

end program main
