module rhs_mod

  implicit none

  integer, parameter, private :: dp = selected_real_kind(15)

  contains

  subroutine rhs (t, y, f)
      implicit none

      real(kind=dp), intent(in) :: t, y(:)
      real(kind=dp), intent(out) :: f(:)

      f = 2.0d0*t
  end subroutine rhs

end module rhs_mod
