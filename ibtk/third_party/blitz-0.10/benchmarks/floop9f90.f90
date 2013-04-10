
      SUBROUTINE floop9_F90(N, a, x)
      INTEGER i, N
      REAL*4 a(N), x(N)

      x = x+a
      RETURN
      END


      SUBROUTINE floop9_F90Overhead(N, a, x)
      INTEGER i, N
      REAL*4 a(N), x(N)

      RETURN
      END
