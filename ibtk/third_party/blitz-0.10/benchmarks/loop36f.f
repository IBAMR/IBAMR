
      SUBROUTINE loop36_F77(N, x, e)
      INTEGER i, N
      REAL*8 x(N), e(N)

      DO i=1,N
          x(i) = exp(e(i));
      END DO
      RETURN
      END


      SUBROUTINE loop36_F77Overhead(N, x, e)
      INTEGER i, N
      REAL*8 x(N), e(N)
      RETURN
      END
