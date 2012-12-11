! module file describing the C++-to-Fortran interface found in f_qd.cpp.

module qdext
  implicit none

  interface
    pure subroutine f_qd_add(a, b, c)
      real*8, intent(in) :: a(4), b(4)
      real*8, intent(out) :: c(4)
    end subroutine

    pure subroutine f_qd_add_qd_d(a, b, c)
      real*8, intent(in) :: a(4), b
      real*8, intent(out) :: c(4)
    end subroutine

    pure subroutine f_qd_add_d_qd(a, b, c)
      real*8, intent(in) :: a, b(4)
      real*8, intent(out) :: c(4)
    end subroutine

    pure subroutine f_qd_sub(a, b, c)
      real*8, intent(in) :: a(4), b(4)
      real*8, intent(out) :: c(4)
    end subroutine

    pure subroutine f_qd_sub_qd_d(a, b, c)
      real*8, intent(in) :: a(4), b
      real*8, intent(out) :: c(4)
    end subroutine

    pure subroutine f_qd_sub_d_qd(a, b, c)
      real*8, intent(in) :: a, b(4)
      real*8, intent(out) :: c(4)
    end subroutine

    pure subroutine f_qd_mul(a, b, c)
      real*8, intent(in) :: a(4), b(4)
      real*8, intent(out) :: c(4)
    end subroutine

    pure subroutine f_qd_mul_qd_d(a, b, c)
      real*8, intent(in) :: a(4), b
      real*8, intent(out) :: c(4)
    end subroutine

    pure subroutine f_qd_mul_d_qd(a, b, c)
      real*8, intent(in) :: a, b(4)
      real*8, intent(out) :: c(4)
    end subroutine

    pure subroutine f_qd_div(a, b, c)
      real*8, intent(in) :: a(4), b(4)
      real*8, intent(out) :: c(4)
    end subroutine

    pure subroutine f_qd_div_qd_d(a, b, c)
      real*8, intent(in) :: a(4), b
      real*8, intent(out) :: c(4)
    end subroutine

    pure subroutine f_qd_div_d_qd(a, b, c)
      real*8, intent(in) :: a, b(4)
      real*8, intent(out) :: c(4)
    end subroutine

    pure subroutine f_qd_sqrt(a, b)
      real*8, intent(in) :: a(4)
      real*8, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qd_sqr(a, b)
      real*8, intent(in) :: a(4)
      real*8, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qd_abs(a, b)
      real*8, intent(in) :: a(4)
      real*8, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qd_npwr(a, n, b)
      real*8, intent(in) :: a(4)
      integer, intent(in) :: n
      real*8, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qd_nroot(a, n, b)
      real*8, intent(in) :: a(4)
      integer, intent(in) :: n
      real*8, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qd_nint(a, b)
      real*8, intent(in) :: a(4)
      real*8, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qd_aint(a, b)
      real*8, intent(in) :: a(4)
      real*8, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qd_floor(a, b)
      real*8, intent(in) :: a(4)
      real*8, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qd_ceil(a, b)
      real*8, intent(in) :: a(4)
      real*8, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qd_log(a, b)
      real*8, intent(in) :: a(4)
      real*8, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qd_log10(a, b)
      real*8, intent(in) :: a(4)
      real*8, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qd_exp(a, b)
      real*8, intent(in) :: a(4)
      real*8, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qd_sin(a, b)
      real*8, intent(in) :: a(4)
      real*8, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qd_cos(a, b)
      real*8, intent(in) :: a(4)
      real*8, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qd_tan(a, b)
      real*8, intent(in) :: a(4)
      real*8, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qd_asin(a, b)
      real*8, intent(in) :: a(4)
      real*8, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qd_acos(a, b)
      real*8, intent(in) :: a(4)
      real*8, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qd_atan(a, b)
      real*8, intent(in) :: a(4)
      real*8, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qd_atan2(a, b, c)
      real*8, intent(in) :: a(4), b(4)
      real*8, intent(out) :: c(4)
    end subroutine

    pure subroutine f_qd_sinh(a, b)
      real*8, intent(in) :: a(4)
      real*8, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qd_cosh(a, b)
      real*8, intent(in) :: a(4)
      real*8, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qd_tanh(a, b)
      real*8, intent(in) :: a(4)
      real*8, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qd_asinh(a, b)
      real*8, intent(in) :: a(4)
      real*8, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qd_acosh(a, b)
      real*8, intent(in) :: a(4)
      real*8, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qd_atanh(a, b)
      real*8, intent(in) :: a(4)
      real*8, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qd_sincos(a, s, c)
      real*8, intent(in) :: a(4)
      real*8, intent(out) :: s(4), c(4)
    end subroutine

    pure subroutine f_qd_sincosh(a, s, c)
      real*8, intent(in) :: a(4)
      real*8, intent(out) :: s(4), c(4)
    end subroutine

    subroutine f_qd_swrite(a, prec, str, maxlen)
      real*8, intent(in) :: a(4)
      integer, intent(in) :: prec, maxlen
      character, intent(out) :: str(*)
    end subroutine

    subroutine f_qd_rand(a)
      real*8, intent(out) :: a(4)
    end subroutine

    pure subroutine f_qd_comp(a, b, r)
      real*8, intent(in) :: a(4), b(4)
      integer, intent(out) :: r
    end subroutine

    pure subroutine f_qd_comp_qd_d(a, b, r)
      real*8, intent(in) :: a(4), b
      integer, intent(out) :: r
    end subroutine

    pure subroutine f_qd_comp_d_qd(a, b, r)
      real*8, intent(in) :: a, b(4)
      integer, intent(out) :: r
    end subroutine

    pure subroutine f_qd_pi(a)
      real*8, intent(out) :: a(4)
    end subroutine

    pure subroutine f_qd_nan(a)
      real*8, intent(out) :: a(4)
    end subroutine

  end interface
end
