
      SUBROUTINE floop16_F90(N, y, x, a, b, c, u)
      INTEGER i, N
      REAL*4 y(N), x(N), a(N), b(N), c(N), u

      x = a+b+c; y = x+c+u
      RETURN
      END


      SUBROUTINE floop16_F90Overhead(N, y, x, a, b, c, u)
      INTEGER i, N
      REAL*4 y(N), x(N), a(N), b(N), c(N), u

      RETURN
      END
