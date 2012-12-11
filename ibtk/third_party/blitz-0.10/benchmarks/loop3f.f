
      SUBROUTINE loop3_F77(N, x, y, a)
      INTEGER i, N
      REAL*8 x(N), y(N), a

      DO i=1,N
          y(i) = y(i) + a*x(i);
      END DO
      RETURN
      END


      SUBROUTINE loop3_F77Overhead(N, x, y, a)
      INTEGER i, N
      REAL*8 x(N), y(N), a
      RETURN
      END
