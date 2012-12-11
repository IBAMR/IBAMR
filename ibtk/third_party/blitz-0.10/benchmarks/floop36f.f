
      SUBROUTINE floop36_F77(N, x, e)
      INTEGER i, N
      REAL*4 x(N), e(N)

      DO i=1,N
          x(i) = exp(e(i));
      END DO
      RETURN
      END


      SUBROUTINE floop36_F77Overhead(N, x, e)
      INTEGER i, N
      REAL*4 x(N), e(N)
      RETURN
      END
