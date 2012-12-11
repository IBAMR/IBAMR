
      SUBROUTINE loop9_F77(N, a, x)
      INTEGER i, N
      REAL*8 a(N), x(N)

      DO i=1,N
          x(i) = x(i)+a(i);
      END DO
      RETURN
      END


      SUBROUTINE loop9_F77Overhead(N, a, x)
      INTEGER i, N
      REAL*8 a(N), x(N)
      RETURN
      END
