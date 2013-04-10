
      SUBROUTINE floop23_F90(N, x, a, b, c, d, y)
      INTEGER i, N
      REAL*4 x(N), a(N), b(N), c(N), d(N), y(N)

      x = a*b + c*d; y = b+d
      RETURN
      END


      SUBROUTINE floop23_F90Overhead(N, x, a, b, c, d, y)
      INTEGER i, N
      REAL*4 x(N), a(N), b(N), c(N), d(N), y(N)

      RETURN
      END
