
      SUBROUTINE floop22_F90(N, x, a, b, c, d)
      INTEGER i, N
      REAL*4 x(N), a(N), b(N), c(N), d(N)

      x = x + a*b + c*d
      RETURN
      END


      SUBROUTINE floop22_F90Overhead(N, x, a, b, c, d)
      INTEGER i, N
      REAL*4 x(N), a(N), b(N), c(N), d(N)

      RETURN
      END
