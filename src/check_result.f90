module check_result_mod

  implicit none
  contains

    pure subroutine check_result(y)
      use iso_fortran_env, only: real64
      implicit none

      real(real64), intent(in), dimension(:) :: y
      integer :: i, flag
      real(real64) :: delta, correct
      real(real64), parameter :: eps = 1.0d-3

      flag = 0

      do i = 1, size(y)
        correct = 25.0d0 + real(i, real64)
        if (abs(y(i) - correct)/correct > eps) flag = flag - 1
      end do

      if (flag < 0) error stop 'ERROR: Wrong integration result'

    end subroutine check_result

end module check_result_mod
