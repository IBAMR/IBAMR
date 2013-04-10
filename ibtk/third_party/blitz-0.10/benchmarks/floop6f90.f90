
      SUBROUTINE floop6_F90(N, x, a, b)
      INTEGER i, N
      REAL*4 x(N), a(N), b(N)

      x = a*b
      RETURN
      END


      SUBROUTINE floop6_F90Overhead(N, x, a, b)
      INTEGER i, N
      REAL*4 x(N), a(N), b(N)

      RETURN
      END
