
      SUBROUTINE floop14_F90(N, y, x, a, b)
      INTEGER i, N
      REAL*4 y(N), x(N), a(N), b(N)

      x = a+b; y = a-b
      RETURN
      END


      SUBROUTINE floop14_F90Overhead(N, y, x, a, b)
      INTEGER i, N
      REAL*4 y(N), x(N), a(N), b(N)

      RETURN
      END
