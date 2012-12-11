
      SUBROUTINE floop12_F90(N, y, x, a, b, c, d, u)
      INTEGER i, N
      REAL*4 y(N), x(N), a(N), b(N), c(N), d(N), u

      y = u+a; x = a+b+c+d
      RETURN
      END


      SUBROUTINE floop12_F90Overhead(N, y, x, a, b, c, d, u)
      INTEGER i, N
      REAL*4 y(N), x(N), a(N), b(N), c(N), d(N), u

      RETURN
      END
