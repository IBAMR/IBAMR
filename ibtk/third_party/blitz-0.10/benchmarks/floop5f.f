
      SUBROUTINE floop5_F77(N, x, a, b)
      INTEGER i, N
      REAL*4 x(N), a(N), b(N)

      DO i=1,N
          x(i) = a(i)+b(i);
      END DO
      RETURN
      END


      SUBROUTINE floop5_F77Overhead(N, x, a, b)
      INTEGER i, N
      REAL*4 x(N), a(N), b(N)
      RETURN
      END
