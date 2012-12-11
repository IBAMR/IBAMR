
      SUBROUTINE loop1_F77(N, x, y)
      INTEGER i, N
      REAL*8 x(N), y(N)

      DO i=1,N
          x(i) = sqrt(y(i));
      END DO
      RETURN
      END


      SUBROUTINE loop1_F77Overhead(N, x, y)
      INTEGER i, N
      REAL*8 x(N), y(N)
      RETURN
      END
