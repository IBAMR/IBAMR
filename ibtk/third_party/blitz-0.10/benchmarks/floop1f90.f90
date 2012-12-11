
      SUBROUTINE floop1_F90(N, x, y)
      INTEGER i, N
      REAL*4 x(N), y(N)

      x = sqrt(y)
      RETURN
      END


      SUBROUTINE floop1_F90Overhead(N, x, y)
      INTEGER i, N
      REAL*4 x(N), y(N)

      RETURN
      END
