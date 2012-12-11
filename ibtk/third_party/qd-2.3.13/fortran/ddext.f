! module file describing the C++-to-Fortran interface found in f_dd.cpp.

module ddext
  implicit none

  interface
    pure subroutine f_dd_add(a, b, c)
      real*8, intent(in) :: a(2), b(2)
      real*8, intent(out) :: c(2)
    end subroutine

    pure subroutine f_dd_add_dd_d(a, b, c)
      real*8, intent(in) :: a(2), b
      real*8, intent(out) :: c(2)
    end subroutine

    pure subroutine f_dd_sub(a, b, c)
      real*8, intent(in) :: a(2), b(2)
      real*8, intent(out) :: c(2)
    end subroutine

    pure subroutine f_dd_sub_dd_d(a, b, c)
      real*8, intent(in) :: a(2), b
      real*8, intent(out) :: c(2)
    end subroutine

    pure subroutine f_dd_sub_d_dd(a, b, c)
      real*8, intent(in) :: a, b(2)
      real*8, intent(out) :: c(2)
    end subroutine

    pure subroutine f_dd_mul(a, b, c)
      real*8, intent(in) :: a(2), b(2)
      real*8, intent(out) :: c(2)
    end subroutine

    pure subroutine f_dd_mul_dd_d(a, b, c)
      real*8, intent(in) :: a(2), b
      real*8, intent(out) :: c(2)
    end subroutine

    pure subroutine f_dd_div(a, b, c)
      real*8, intent(in) :: a(2), b(2)
      real*8, intent(out) :: c(2)
    end subroutine

    pure subroutine f_dd_div_dd_d(a, b, c)
      real*8, intent(in) :: a(2), b
      real*8, intent(out) :: c(2)
    end subroutine

    pure subroutine f_dd_div_d_dd(a, b, c)
      real*8, intent(in) :: a, b(2)
      real*8, intent(out) :: c(2)
    end subroutine

    pure subroutine f_dd_sqrt(a, b)
      real*8, intent(in) :: a(2)
      real*8, intent(out) :: b(2)
    end subroutine

    pure subroutine f_dd_sqr(a, b)
      real*8, intent(in) :: a(2)
      real*8, intent(out) :: b(2)
    end subroutine

    pure subroutine f_dd_abs(a, b)
      real*8, intent(in) :: a(2)
      real*8, intent(out) :: b(2)
    end subroutine

    pure subroutine f_dd_npwr(a, n, b)
      real*8, intent(in) :: a(2)
      integer, intent(in) :: n
      real*8, intent(out) :: b(2)
    end subroutine

    pure subroutine f_dd_nroot(a, n, b)
      real*8, intent(in) :: a(2)
      integer, intent(in) :: n
      real*8, intent(out) :: b(2)
    end subroutine

    pure subroutine f_dd_nint(a, b)
      real*8, intent(in) :: a(2)
      real*8, intent(out) :: b(2)
    end subroutine

    pure subroutine f_dd_aint(a, b)
      real*8, intent(in) :: a(2)
      real*8, intent(out) :: b(2)
    end subroutine

    pure subroutine f_dd_floor(a, b)
      real*8, intent(in) :: a(2)
      real*8, intent(out) :: b(2)
    end subroutine

    pure subroutine f_dd_ceil(a, b)
      real*8, intent(in) :: a(2)
      real*8, intent(out) :: b(2)
    end subroutine

    pure subroutine f_dd_log(a, b)
      real*8, intent(in) :: a(2)
      real*8, intent(out) :: b(2)
    end subroutine

    pure subroutine f_dd_log10(a, b)
      real*8, intent(in) :: a(2)
      real*8, intent(out) :: b(2)
    end subroutine

    pure subroutine f_dd_exp(a, b)
      real*8, intent(in) :: a(2)
      real*8, intent(out) :: b(2)
    end subroutine

    pure subroutine f_dd_sin(a, b)
      real*8, intent(in) :: a(2)
      real*8, intent(out) :: b(2)
    end subroutine

    pure subroutine f_dd_cos(a, b)
      real*8, intent(in) :: a(2)
      real*8, intent(out) :: b(2)
    end subroutine

    pure subroutine f_dd_tan(a, b)
      real*8, intent(in) :: a(2)
      real*8, intent(out) :: b(2)
    end subroutine

    pure subroutine f_dd_asin(a, b)
      real*8, intent(in) :: a(2)
      real*8, intent(out) :: b(2)
    end subroutine

    pure subroutine f_dd_acos(a, b)
      real*8, intent(in) :: a(2)
      real*8, intent(out) :: b(2)
    end subroutine

    pure subroutine f_dd_atan(a, b)
      real*8, intent(in) :: a(2)
      real*8, intent(out) :: b(2)
    end subroutine

    pure subroutine f_dd_atan2(a, b, c)
      real*8, intent(in) :: a(2), b(2)
      real*8, intent(out) :: c(2)
    end subroutine

    pure subroutine f_dd_sinh(a, b)
      real*8, intent(in) :: a(2)
      real*8, intent(out) :: b(2)
    end subroutine

    pure subroutine f_dd_cosh(a, b)
      real*8, intent(in) :: a(2)
      real*8, intent(out) :: b(2)
    end subroutine

    pure subroutine f_dd_tanh(a, b)
      real*8, intent(in) :: a(2)
      real*8, intent(out) :: b(2)
    end subroutine

    pure subroutine f_dd_asinh(a, b)
      real*8, intent(in) :: a(2)
      real*8, intent(out) :: b(2)
    end subroutine

    pure subroutine f_dd_acosh(a, b)
      real*8, intent(in) :: a(2)
      real*8, intent(out) :: b(2)
    end subroutine

    pure subroutine f_dd_atanh(a, b)
      real*8, intent(in) :: a(2)
      real*8, intent(out) :: b(2)
    end subroutine

    pure subroutine f_dd_sincos(a, s, c)
      real*8, intent(in) :: a(2)
      real*8, intent(out) :: s(2), c(2)
    end subroutine

    pure subroutine f_dd_sincosh(a, s, c)
      real*8, intent(in) :: a(2)
      real*8, intent(out) :: s(2), c(2)
    end subroutine

    subroutine f_dd_swrite(a, prec, str, maxlen)
      real*8, intent(in) :: a(2)
      integer, intent(in) :: prec, maxlen
      character, intent(out) :: str(maxlen)
    end subroutine

    subroutine f_dd_rand(a)
      real*8, intent(out) :: a(2)
    end subroutine

    pure subroutine f_dd_comp(a, b, r)
      real*8, intent(in) :: a(2), b(2)
      integer, intent(out) :: r
    end subroutine

    pure subroutine f_dd_comp_dd_d(a, b, r)
      real*8, intent(in) :: a(2), b
      integer, intent(out) :: r
    end subroutine

    pure subroutine f_dd_comp_d_dd(a, b, r)
      real*8, intent(in) :: a, b(2)
      integer, intent(out) :: r
    end subroutine

    pure subroutine f_dd_pi(a)
      real*8, intent(out) :: a(2)
    end subroutine

    pure subroutine f_dd_nan(a)
      real*8, intent(out) :: a(2)
    end subroutine

  end interface
end
