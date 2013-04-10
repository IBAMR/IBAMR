
      SUBROUTINE loop13_F90(N, y, x, a, b, c, d, u)
      INTEGER i, N
      REAL*8 y(N), x(N), a(N), b(N), c(N), d(N), u

      x = a+b+c+d; y = u+d
      RETURN
      END


      SUBROUTINE loop13_F90Overhead(N, y, x, a, b, c, d, u)
      INTEGER i, N
      REAL*8 y(N), x(N), a(N), b(N), c(N), d(N), u

      RETURN
      END
