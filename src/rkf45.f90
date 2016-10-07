module rkf45_mod
  implicit none

  contains
    subroutine rkf45 (y0, t0, tf, dt0, eps, t, dt, y, flag, num_steps)

        use rhs_mod
        implicit none

        integer, parameter :: dp = selected_real_kind(15)

        real(kind=dp), intent(in) :: y0(:), t0, tf, dt0, eps
        real(kind=dp), intent(out) :: t, dt, y(:)
        integer, intent(out) :: flag, num_steps

        real(kind=dp) :: k1(size(y0))
        real(kind=dp) :: k2(size(y0))
        real(kind=dp) :: k3(size(y0))
        real(kind=dp) :: k4(size(y0))
        real(kind=dp) :: k5(size(y0))
        real(kind=dp) :: k6(size(y0))

        real(kind=dp) :: y_rk4_approx(size(y0))
        real(kind=dp) :: y_rk5_approx(size(y0))
        real(kind=dp) :: s
        real(kind=dp) :: lte(size(y0))
        real(kind=dp) :: maxlte
        real(kind=dp) :: f(size(y0))

        integer, parameter :: maxiter = 10
        integer :: i

        num_steps = 0

        t = t0
        y = y0

        do while (t < tf)
          num_steps = num_steps + 1
          ! Initial guess at integration step size.
          dt = dt0
          ! Keep adjusting step size until we run out of iterations.
          do i = 1, maxiter

            ! A bunch of coefficients, taken from Hairer, Norsett, and Wanner (1993).

            call rhs(t, y, f)
            k1 = dt * f

            call rhs(t + 1.0d0/4.0d0*dt, y + 1.0d0/4.0d0*k1, f)
            k2 = dt * f

            call rhs(t + 3.0d0/8.0d0*dt, y + 3.0d0/32.0d0*k1 &
                                           + 9.0d0/32.0d0*k2, &
                     f)
            k3 = dt * f

            call rhs(t + 12.0d0/13.0d0*dt, y + 1932.0d0/2197.0d0*k1 &
                                             - 7200.0d0/2197.0d0*k2 &
                                             + 7296.0d0/2197.0d0*k3, &
                     f)
            k4 = dt * f

            call rhs(t + dt, y + 439.0d0/216.0d0*k1 &
                               - 8.0d0*k2 &
                               + 3680.0d0/513.0d0*k3 &
                               - 845.0d0/4104.0d0*k4, &
                     f)
            k5 = dt * f

            call rhs(t + 1.0d0/2.0d0*dt, y - 8.0d0/27.0d0*k1 &
                                           + 2.0d0*k2 &
                                           - 3544.0d0/2565.0d0*k3 &
                                           + 1859.0d0/4104.0d0*k4 &
                                           - 11.0d0/40.0d0*k5, &
                     f)
            k6 = dt * f

            ! 4th-order approximation to y(t+dt)
            y_rk4_approx = y + 25.0d0/216.0d0*k1 &
                              + 1408.0d0/2565.0d0*k3 &
                              + 2197.0d0/4101.0d0*k4 &
                              - 1.0d0/5.0d0*k5
            ! 5th-order approximation to y(t+dt)
            y_rk5_approx = y + 16.0d0/135.0d0*k1 &
                              + 6656.0d0/12825.0d0*k3 &
                              + 28561.0d0/56430.0d0*k4 &
                              - 9.0d0/50.0d0*k5 &
                              + 2.0d0/55.0d0*k6

            ! Local truncation error between 4th- and 5th-order approximations of y(t+dt)
            lte = abs(y_rk5_approx - y_rk4_approx)
            maxlte = maxval(lte)

            ! Scalar by which to multiply step size to get new step size to satisfy error tolerance
            s = (eps / (2.0d0 * maxlte))**(1.0d0/4.0d0)

            ! If the local truncation error satisfies the tolerance
            ! requirements, then we're done!
            if (maxlte < eps) then
                flag = 0
                exit
            end if

            ! Adjust step size
            dt = dt*s

            ! Make sure we don't go past the final value
            if (t + dt > tf) then
                dt = tf - t
                cycle
            end if

          end do

          if (maxlte > eps) flag = -1

          t = t + dt
          y = y_rk5_approx

        end do

    end subroutine rkf45

end module rkf45_mod
