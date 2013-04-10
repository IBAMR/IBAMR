
      SUBROUTINE loop2_F90(N, x, y, u)
      INTEGER i, N
      REAL*8 x(N), y(N), u

      x = y/u
      RETURN
      END


      SUBROUTINE loop2_F90Overhead(N, x, y, u)
      INTEGER i, N
      REAL*8 x(N), y(N), u

      RETURN
      END
