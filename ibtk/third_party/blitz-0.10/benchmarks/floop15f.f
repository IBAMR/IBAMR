
      SUBROUTINE floop15_F77(N, x, a, b, c)
      INTEGER i, N
      REAL*4 x(N), a(N), b(N), c(N)

      DO i=1,N
          x(i) = c(i) + a(i)*b(i);
      END DO
      RETURN
      END


      SUBROUTINE floop15_F77Overhead(N, x, a, b, c)
      INTEGER i, N
      REAL*4 x(N), a(N), b(N), c(N)
      RETURN
      END
