
      SUBROUTINE floop8_F90(N, a, x, u)
      INTEGER i, N
      REAL*4 a(N), x(N), u

      x = u/a
      RETURN
      END


      SUBROUTINE floop8_F90Overhead(N, a, x, u)
      INTEGER i, N
      REAL*4 a(N), x(N), u

      RETURN
      END
