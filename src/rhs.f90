module rhs_mod

  implicit none

  integer, parameter, private :: dp = selected_real_kind(15)

  contains

  pure elemental real(kind=dp) function f (t, y)
      implicit none

      real(kind=dp), intent(in) :: t, y

      f = 2.0d0*t
  end function f

end module rhs_mod
