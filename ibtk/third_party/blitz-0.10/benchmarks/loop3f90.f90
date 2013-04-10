
      SUBROUTINE loop3_F90(N, x, y, a)
      INTEGER i, N
      REAL*8 x(N), y(N), a

      y = y + a*x
      RETURN
      END


      SUBROUTINE loop3_F90Overhead(N, x, y, a)
      INTEGER i, N
      REAL*8 x(N), y(N), a

      RETURN
      END
