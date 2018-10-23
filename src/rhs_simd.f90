module rhs_simd_mod

  use iso_fortran_env, only: real64
  implicit none

  contains

  pure subroutine rhs_simd (t, y, f)
      implicit none

      real(real64), intent(in) :: t, y(:)
      real(real64), intent(out) :: f(:)

      f = 2.0d0*t
  end subroutine rhs_simd

end module rhs_simd_mod
