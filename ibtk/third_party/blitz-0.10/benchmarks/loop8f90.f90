
      SUBROUTINE loop8_F90(N, a, x, u)
      INTEGER i, N
      REAL*8 a(N), x(N), u

      x = u/a
      RETURN
      END


      SUBROUTINE loop8_F90Overhead(N, a, x, u)
      INTEGER i, N
      REAL*8 a(N), x(N), u

      RETURN
      END
