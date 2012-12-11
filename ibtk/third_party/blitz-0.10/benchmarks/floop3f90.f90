
      SUBROUTINE floop3_F90(N, x, y, a)
      INTEGER i, N
      REAL*4 x(N), y(N), a

      y = y + a*x
      RETURN
      END


      SUBROUTINE floop3_F90Overhead(N, x, y, a)
      INTEGER i, N
      REAL*4 x(N), y(N), a

      RETURN
      END
