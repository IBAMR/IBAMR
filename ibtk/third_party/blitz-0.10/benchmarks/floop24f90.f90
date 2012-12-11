
      SUBROUTINE floop24_F90(N, x, a, b, c, d, y)
      INTEGER i, N
      REAL*4 x(N), a(N), b(N), c(N), d(N), y(N)

      x = a*c - b*c; y = a*d + b+c
      RETURN
      END


      SUBROUTINE floop24_F90Overhead(N, x, a, b, c, d, y)
      INTEGER i, N
      REAL*4 x(N), a(N), b(N), c(N), d(N), y(N)

      RETURN
      END
