
      SUBROUTINE loop2_F77(N, x, y, u)
      INTEGER i, N
      REAL*8 x(N), y(N), u

      DO i=1,N
          x(i) = y(i)/u;
      END DO
      RETURN
      END


      SUBROUTINE loop2_F77Overhead(N, x, y, u)
      INTEGER i, N
      REAL*8 x(N), y(N), u
      RETURN
      END
