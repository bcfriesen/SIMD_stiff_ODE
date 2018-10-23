program main

    use rkf45_simd_mod
    use rkf45_scalar_mod
    implicit none

    integer, parameter :: dp = selected_real_kind(15)

    ! The actual width of the array (analogous to # of cells to sweep
    ! through in Nyx).
    integer, parameter :: array_length = 2**19

    ! The width of the array to pass to the SIMD integrator. The
    ! integrator supports arrays of arbitrary width, but there is likely
    ! a sweet spot for the best with to give it. The minimum should be
    ! the SIMD width of the architecture (for HSW it's 256 bits -> 4
    ! doubles; for Knights it's 512 bits -> 8 doubles). Choosing a
    ! smaller number means the temporary arrays use less memory, but the
    ! code makes more function calls. The converse is true for larger
    ! array widths. Larger array widths may also incur memory access
    ! penalties due to spilling out of cache. So the sweet spot may be
    ! somewhere in the middle.
    integer :: rkf_array_width

    ! Initial and final values of the variable in the ODE. In contrast
    ! to most other integrators, I decided not to make this this one
    ! update the values in situ.
    real(kind=dp) :: y0(array_length), y(array_length)

    ! Initial and final times for the integration. In a Strang splitting
    ! grid sweep in Nyx, all cells have the same initial and final time,
    ! so these are scalars, not arrays.
    real(kind=dp) :: t0, tf

    ! Initial guess for the time step in the integration. In theory we
    ! could have a unique guess for each cell, but that would be an
    ! awful lot of work. So just make all cells use the same initial
    ! guess.
    real(kind=dp) :: dt0

    ! Tolerance for the local truncation error, which decides how small
    ! the time steps can be when integrating from t0 to tf.
    real(kind=dp) :: eps

    ! Final value of the time and last time step. t should be just tf,
    ! so there is no good reason to keep it, and it can probably go
    ! away. dt just shows the value of the last time step when arriving
    ! at tf.
    real(kind=dp) :: t, dt

    ! Used for wall clock timers.
    real(kind=dp) :: t1, t2

    ! Number of steps the integration took to get from t0 to tf.
    integer :: num_steps

    ! 0 means integration succeeded; < 0 means it failed.
    integer :: flag, err = 0

    integer :: i

    write (*,'(a20, i8)') 'Array length: ', array_length

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
      err = err + flag
    end do
    call cpu_time(t2)
    write (*, *) 'time in scalar RKF45 (sec): ', t2-t1
    write (*, *)

    ! Now do the same sweep of integrations using the SIMD RKF45, but
    ! repeat with varying array widths given to RKF.
    rkf_array_width = 2
    do while (rkf_array_width <= array_length)
      write (*,'(a40, i8)') 'Array width per call to RKF45: ', rkf_array_width
      call cpu_time(t1)
      do i = 1, array_length, rkf_array_width
        call rkf45_simd(y0(i:i+rkf_array_width-1), t0, tf, dt0, eps, t, dt, y(i:i+rkf_array_width-1), flag, num_steps)
        err = err + flag
      end do
      call cpu_time(t2)
      write (*, *) 'time in SIMD RKF45 (sec): ', t2-t1
      rkf_array_width = rkf_array_width*2
    end do

    if (err < 0) error stop 'ERROR: Integration failed!'

    write (*, *) 'Success!'
    write (*, '(a50, 2i8)') '# of time steps taken: ', num_steps

end program main
